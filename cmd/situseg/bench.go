package main

// situseg simulate and situseg score: a simulated tissue written as a Xenium-shaped table of
// transcripts with the truth beside it, and a segmentation's transcripts scored against it.

import (
	"bufio"
	"compress/gzip"
	"context"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/robert-mcdermott/situseg"
	"github.com/robert-mcdermott/situseg/internal/parquet"
	"github.com/robert-mcdermott/situseg/internal/simulate"
)

func simulateCommand(args []string) error {
	fs := flag.NewFlagSet("simulate", flag.ExitOnError)
	preset := fs.String("preset", "small", "small (160 × 120 µm, about 300 cells) or standard (500 × 500 µm, about 3,000 cells)")
	seed := fs.Uint64("seed", 0, "the seed (default: the preset's)")
	width := fs.Float64("width-um", 0, "the tissue's width (µm; default: the preset's)")
	height := fs.Float64("height-um", 0, "the tissue's height (µm; default: the preset's)")
	fs.Usage = func() {
		fmt.Fprintln(fs.Output(), `Usage: situseg simulate [flags] transcripts.parquet|transcripts.csv.gz

Writes a simulated tissue's transcripts as Xenium's transcripts table (the columns situseg reads),
with the truth in more columns: each molecule's true cell (true_cell_id, empty for the
background), its origin and its true place. The same preset and seed give the same tissue on any
computer. Score a segmentation of it with situseg score.

Flags:`)
		fs.PrintDefaults()
	}
	fs.Parse(args)
	if fs.NArg() != 1 {
		fs.Usage()
		os.Exit(2)
	}
	p, ok := simulate.Presets[*preset]
	if !ok {
		return fmt.Errorf("no simulation preset %q (small or standard)", *preset)
	}
	if *seed != 0 {
		p.Seed = *seed
	}
	if *width > 0 {
		p.WidthUm = *width
	}
	if *height > 0 {
		p.HeightUm = *height
	}
	tis, err := simulate.Simulate(context.Background(), p)
	if err != nil {
		return err
	}
	path := fs.Arg(0)
	origins := []string{"cell", "background", "control probe", "control codeword"}
	cellOf := func(ids []string, c int32) string {
		if c < 0 {
			return "UNASSIGNED"
		}
		return ids[c]
	}
	n := len(tis.X)
	switch name := strings.ToLower(path); {
	case strings.HasSuffix(name, ".parquet"):
		f, err := os.Create(path)
		if err != nil {
			return err
		}
		bw := bufio.NewWriterSize(f, 1<<20)
		fields := []parquet.Field{{Name: "transcript_id", Type: parquet.Int64, Annotation: "uint"}, {Name: "cell_id", Type: parquet.ByteArray, Annotation: "string"},
			{Name: "overlaps_nucleus", Type: parquet.Int32, Annotation: "uint8"}, {Name: "feature_name", Type: parquet.ByteArray, Annotation: "string"},
			{Name: "x_location", Type: parquet.Float}, {Name: "y_location", Type: parquet.Float}, {Name: "z_location", Type: parquet.Float}, {Name: "qv", Type: parquet.Float},
			{Name: "true_cell_id", Type: parquet.ByteArray, Annotation: "string"}, {Name: "origin", Type: parquet.ByteArray, Annotation: "string"},
			{Name: "true_x_location", Type: parquet.Float}, {Name: "true_y_location", Type: parquet.Float}}
		w, err := parquet.NewWriter(bw, fields, "situseg "+version())
		if err != nil {
			f.Close()
			return err
		}
		ids, nuc := &parquet.Values{Len: n, Int64: make([]int64, n)}, &parquet.Values{Len: n, Int32: make([]int32, n)}
		cells, feats := &parquet.Values{Len: n, Bytes: make([][]byte, n)}, &parquet.Values{Len: n, Bytes: make([][]byte, n)}
		truth, origin := &parquet.Values{Len: n, Bytes: make([][]byte, n)}, &parquet.Values{Len: n, Bytes: make([][]byte, n)}
		for i := range n {
			ids.Int64[i] = int64(tis.ID[i])
			if tis.Nuclear[i] {
				nuc.Int32[i] = 1
			}
			cells.Bytes[i] = []byte(cellOf(tis.CellIDs, tis.Cell[i]))
			feats.Bytes[i] = []byte(tis.Features[tis.Feature[i]])
			truth.Bytes[i] = []byte{}
			if c := tis.TrueCell[i]; c >= 0 {
				truth.Bytes[i] = []byte(tis.TrueIDs[c])
			}
			origin.Bytes[i] = []byte(origins[tis.Origin[i]])
		}
		cols := []*parquet.Values{ids, cells, nuc, feats, {Len: n, Float32: tis.X}, {Len: n, Float32: tis.Y}, {Len: n, Float32: tis.Z}, {Len: n, Float32: tis.QV},
			truth, origin, {Len: n, Float32: tis.TrueX}, {Len: n, Float32: tis.TrueY}}
		if err := w.WriteRowGroup(cols); err != nil {
			f.Close()
			return err
		}
		if err := w.Close(); err != nil {
			f.Close()
			return err
		}
		if err := bw.Flush(); err != nil {
			f.Close()
			return err
		}
		if err := f.Close(); err != nil {
			return err
		}
	case strings.HasSuffix(name, ".csv.gz"), strings.HasSuffix(name, ".csv"):
		f, err := os.Create(path)
		if err != nil {
			return err
		}
		bw := bufio.NewWriterSize(f, 1<<20)
		var w io.Writer = bw
		var zw *gzip.Writer
		if strings.HasSuffix(name, ".gz") {
			zw = gzip.NewWriter(bw)
			w = zw
		}
		b := []byte("transcript_id,cell_id,overlaps_nucleus,feature_name,x_location,y_location,z_location,qv,true_cell_id,origin,true_x_location,true_y_location\n")
		for i := range n {
			truth := ""
			if c := tis.TrueCell[i]; c >= 0 {
				truth = tis.TrueIDs[c]
			}
			nuclear := 0
			if tis.Nuclear[i] {
				nuclear = 1
			}
			b = fmt.Appendf(b, "%d,%s,%d,%s,%g,%g,%g,%g,%s,%s,%g,%g\n", tis.ID[i], cellOf(tis.CellIDs, tis.Cell[i]), nuclear, tis.Features[tis.Feature[i]],
				tis.X[i], tis.Y[i], tis.Z[i], tis.QV[i], truth, origins[tis.Origin[i]], tis.TrueX[i], tis.TrueY[i])
			if len(b) > 1<<16 {
				if _, err := w.Write(b); err != nil {
					f.Close()
					return err
				}
				b = b[:0]
			}
		}
		_, err = w.Write(b)
		if err == nil && zw != nil {
			err = zw.Close()
		}
		if err == nil {
			err = bw.Flush()
		}
		if cerr := f.Close(); err == nil {
			err = cerr
		}
		if err != nil {
			return err
		}
	default:
		return errors.New("the output ends in .parquet, .csv or .csv.gz")
	}
	fmt.Fprintf(os.Stderr, "Wrote %s: %d transcripts, %d true cells, %d instrument cells.\n", path, n, len(tis.TrueIDs), len(tis.CellIDs))
	return nil
}

func scoreCommand(args []string) error {
	fs := flag.NewFlagSet("score", flag.ExitOnError)
	asJSON := fs.Bool("json", false, "print the scores as JSON")
	fs.Usage = func() {
		fmt.Fprintln(fs.Output(), `Usage: situseg score [flags] truth.parquet|truth.csv.gz segmentation ...

Scores segmentations of a simulated tissue (situseg simulate) against its truth, and the
instrument's cells (the table's cell_id) for comparison: the share of the gene transcripts of
quality 20 or more that the truth puts in a cell and the segmentation puts in the same cell (cells
matched one to one where each is the other's best overlap), and the matched cells' median
Jaccard index. A segmentation is situseg's output folder, or a table of transcript_id and cell_id
(CSV or Parquet).

Flags:`)
		fs.PrintDefaults()
	}
	fs.Parse(args)
	if fs.NArg() < 2 {
		fs.Usage()
		os.Exit(2)
	}
	truthCols := situseg.XeniumColumns()
	inst, err := situseg.ReadTable(fs.Arg(0), truthCols)
	if err != nil {
		return err
	}
	truthCols.Cell, truthCols.NoCell = "true_cell_id", []string{"", "UNASSIGNED"}
	truth, err := situseg.ReadTable(fs.Arg(0), truthCols)
	if err != nil {
		return err
	}
	if inst.ID == nil || inst.QV == nil {
		return errors.New("the truth needs its transcript_id and qv columns")
	}
	counted := func(i int) bool { return !truth.Control[truth.Feature[i]] && truth.QV[i] >= 20 }
	scores := map[string]simulate.Score{"instrument": simulate.Compare(truth.Cell, inst.Cell, len(truth.CellIDs), len(inst.CellIDs), counted)}
	order := []string{"instrument"}
	for _, path := range fs.Args()[1:] {
		file := path
		if info, err := os.Stat(path); err == nil && info.IsDir() {
			file = filepath.Join(path, situseg.FileTranscripts)
			if _, err := os.Stat(file); err != nil {
				file = filepath.Join(path, situseg.FileTranscriptsPq)
			}
		}
		header, err := situseg.TableHeader(file)
		if err != nil {
			return err
		}
		has := map[string]bool{}
		for _, h := range header {
			has[h] = true
		}
		cols := situseg.Columns{X: "x", Y: "y", Feature: "feature_name", Cell: "cell_id", ID: "transcript_id", NoCell: []string{"", "UNASSIGNED"}}
		if has["x_location"] {
			cols.X, cols.Y = "x_location", "y_location"
		}
		seg, err := situseg.ReadTable(file, cols)
		if err != nil {
			return err
		}
		if seg.ID == nil {
			return fmt.Errorf("%s has no transcript_id column", file)
		}
		byID := make(map[uint64]int32, seg.Rows())
		for i, id := range seg.ID {
			byID[id] = seg.Cell[i]
		}
		pred := make([]int32, truth.Rows())
		for i, id := range inst.ID {
			c, ok := byID[id]
			if !ok {
				c = -1
			}
			pred[i] = c
		}
		scores[path] = simulate.Compare(truth.Cell, pred, len(truth.CellIDs), len(seg.CellIDs), counted)
		order = append(order, path)
	}
	if *asJSON {
		raw, _ := json.MarshalIndent(scores, "", " ")
		fmt.Println(string(raw))
		return nil
	}
	for _, name := range order {
		s := scores[name]
		fmt.Printf("%-30s %5.1f%% in their true cell, median Jaccard %.3f, %d of %d true cells matched, %s in a cell (the truth %s)\n",
			name, 100*s.InTrueCell, s.MedianJaccard, s.MatchedCells, s.TrueCells, percent(s.Assigned), percent(s.TrueAssigned))
	}
	return nil
}

func percent(x float64) string { return strconv.FormatFloat(100*x, 'f', 1, 64) + "%" }
