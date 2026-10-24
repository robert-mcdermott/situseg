package situseg

// Reading tables: every format (Parquet, CSV, TSV, gzipped), each instrument's columns, the cell
// column's no-cell values and groups, the nucleus column's kinds, the controls; and the inputs a
// table gives.

import (
	"bytes"
	"compress/gzip"
	"context"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/robert-mcdermott/situseg/internal/parquet"
)

func writeText(t *testing.T, name, text string) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	data := []byte(text)
	if strings.HasSuffix(name, ".gz") {
		var buf bytes.Buffer
		zw := gzip.NewWriter(&buf)
		zw.Write(data)
		zw.Close()
		data = buf.Bytes()
	}
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestReadXeniumCSV(t *testing.T) {
	text := "\uFEFFtranscript_id,cell_id,overlaps_nucleus,feature_name,x_location,y_location,z_location,qv\n" +
		"1,c1,1,GENEA,1.5,2.5,3,40\n" +
		"2,UNASSIGNED,0,NegControlProbe_X,10,20,1,30\n" +
		"3,c2,0,GENEB,4,5,6,12.5\n" +
		"4,c1,0,\"GENE,C\",7,8,9,39\n"
	for _, name := range []string{"transcripts.csv", "transcripts.csv.gz"} {
		path := writeText(t, name, text)
		header, err := TableHeader(path)
		if err != nil {
			t.Fatal(err)
		}
		cols, inst, ok := DetectColumns(header)
		if !ok || inst != "xenium" {
			t.Fatalf("%s: detected %q (%v) from %v", name, inst, ok, header)
		}
		tab, err := ReadTable(path, cols)
		if err != nil {
			t.Fatal(err)
		}
		if tab.Rows() != 4 || len(tab.CellIDs) != 2 || len(tab.Features) != 4 {
			t.Fatalf("%s: %d rows, cells %v, features %v", name, tab.Rows(), tab.CellIDs, tab.Features)
		}
		if tab.Cell[1] != -1 || tab.CellIDs[tab.Cell[3]] != "c1" || !tab.Nuclear[0] || tab.Nuclear[2] {
			t.Errorf("%s: cells %v, nuclear %v", name, tab.Cell, tab.Nuclear)
		}
		if !tab.Control[tab.Feature[1]] || tab.Control[tab.Feature[0]] || tab.Features[tab.Feature[3]] != "GENE,C" {
			t.Errorf("%s: controls %v of %v", name, tab.Control, tab.Features)
		}
		if tab.X[0] != 1.5 || tab.Z[2] != 6 || tab.QV[2] != 12.5 || tab.ID[3] != 4 {
			t.Errorf("%s: x %v z %v qv %v ids %v", name, tab.X, tab.Z, tab.QV, tab.ID)
		}
	}
}

// CosMx numbers cells within each field of view, in pixels; MERSCOPE has no nucleus column and
// numbers its z planes.
func TestReadCosMxAndMERSCOPE(t *testing.T) {
	cosmx := "fov,cell_ID,x_global_px,y_global_px,x_local_px,y_local_px,z,target,CellComp\n" +
		"1,5,100,200,1,2,0,GeneA,Nuclear\n" +
		"2,5,300,400,1,2,1,GeneB,Cytoplasm\n" +
		"2,0,10,20,1,2,1,NegPrb1,0\n"
	path := writeText(t, "tx_file.csv", cosmx)
	header, _ := TableHeader(path)
	cols, inst, ok := DetectColumns(header)
	if !ok || inst != "cosmx" {
		t.Fatalf("detected %q from %v", inst, header)
	}
	tab, err := ReadTable(path, cols)
	if err != nil {
		t.Fatal(err)
	}
	if len(tab.CellIDs) != 2 || tab.CellIDs[0] != "1_5" || tab.CellIDs[1] != "2_5" || tab.Cell[2] != -1 {
		t.Errorf("cells %v %v", tab.CellIDs, tab.Cell)
	}
	if !tab.Nuclear[0] || tab.Nuclear[1] || !tab.Control[tab.Feature[2]] {
		t.Errorf("nuclear %v, controls %v", tab.Nuclear, tab.Control)
	}
	if math.Abs(float64(tab.X[0])-12.028) > 1e-4 || tab.QV != nil {
		t.Errorf("x %v (scaled), qv %v", tab.X, tab.QV)
	}
	merscope := ",barcode_id,global_x,global_y,global_z,x,y,fov,gene,transcript_id,cell_id\n" +
		"0,1,10.5,20.5,0,1,2,0,GeneA,ENST1,123\n" +
		"1,2,11,21,2,1,2,0,Blank-1,ENST2,-1\n"
	path = writeText(t, "detected_transcripts.csv", merscope)
	header, _ = TableHeader(path)
	cols, inst, ok = DetectColumns(header)
	if !ok || inst != "merscope" {
		t.Fatalf("detected %q from %v", inst, header)
	}
	tab, err = ReadTable(path, cols)
	if err != nil {
		t.Fatal(err)
	}
	if tab.Nuclear != nil || tab.Cell[1] != -1 || tab.CellIDs[0] != "123" || tab.Z[1] != 3 || !tab.Control[tab.Feature[1]] {
		t.Errorf("nuclear %v, cells %v %v, z %v, controls %v", tab.Nuclear, tab.Cell, tab.CellIDs, tab.Z, tab.Control)
	}
}

func TestReadTSVAndErrors(t *testing.T) {
	path := writeText(t, "t.tsv", "x\ty\tgene\tcell\n1\t2\tA\tk1\n3\t4\tB\t\n")
	tab, err := ReadTable(path, Columns{X: "x", Y: "y", Feature: "gene", Cell: "cell", NoCell: []string{""}})
	if err != nil {
		t.Fatal(err)
	}
	if tab.Rows() != 2 || tab.Cell[1] != -1 || tab.Z != nil || tab.QV != nil {
		t.Errorf("%+v", tab)
	}
	if _, err := ReadTable(path, Columns{X: "x", Y: "missing", Feature: "gene", Cell: "cell"}); err == nil || !strings.Contains(err.Error(), "missing") {
		t.Errorf("a missing column: %v", err)
	}
	bad := writeText(t, "bad.csv", "x,y,gene,cell\n1,two,A,k\n")
	if _, err := ReadTable(bad, Columns{X: "x", Y: "y", Feature: "gene", Cell: "cell"}); err == nil {
		t.Error("a coordinate that is not a number was read")
	}
	none := writeText(t, "none.csv", "x,y,gene,cell\n1,2,A,UNASSIGNED\n")
	if _, err := ReadTable(none, XeniumColumnsFor("x", "y", "gene", "cell")); err == nil {
		t.Error("a table with no cells was read")
	}
}

// XeniumColumnsFor is Xenium's columns with other names, for the tests.
func XeniumColumnsFor(x, y, gene, cell string) Columns {
	c := XeniumColumns()
	c.X, c.Y, c.Z, c.Feature, c.Cell, c.QV, c.Nuclear, c.ID = x, y, "", gene, cell, "", "", ""
	return c
}

// A Parquet table gives what the same CSV table gives.
func TestReadParquetLikeCSV(t *testing.T) {
	tab := tableOf(simulated(t, "small", 60, 40))
	n := tab.Rows()
	var buf bytes.Buffer
	fields := []parquet.Field{{Name: "transcript_id", Type: parquet.Int64, Annotation: "uint"}, {Name: "cell_id", Type: parquet.ByteArray, Annotation: "string"},
		{Name: "overlaps_nucleus", Type: parquet.Int32, Annotation: "uint8"}, {Name: "feature_name", Type: parquet.ByteArray, Annotation: "string"},
		{Name: "x_location", Type: parquet.Float}, {Name: "y_location", Type: parquet.Float}, {Name: "z_location", Type: parquet.Float}, {Name: "qv", Type: parquet.Float}}
	w, err := parquet.NewWriter(&buf, fields, "situseg test")
	if err != nil {
		t.Fatal(err)
	}
	ids, nuc := &parquet.Values{Len: n, Int64: make([]int64, n)}, &parquet.Values{Len: n, Int32: make([]int32, n)}
	cells, feats := &parquet.Values{Len: n, Bytes: make([][]byte, n)}, &parquet.Values{Len: n, Bytes: make([][]byte, n)}
	var csv strings.Builder
	csv.WriteString("transcript_id,cell_id,overlaps_nucleus,feature_name,x_location,y_location,z_location,qv\n")
	for i := range n {
		ids.Int64[i] = int64(tab.ID[i])
		cell := "UNASSIGNED"
		if c := tab.Cell[i]; c >= 0 {
			cell = tab.CellIDs[c]
		}
		cells.Bytes[i], feats.Bytes[i] = []byte(cell), []byte(tab.Features[tab.Feature[i]])
		if tab.Nuclear[i] {
			nuc.Int32[i] = 1
		}
		csv.WriteString(strings.Join([]string{formatU(tab.ID[i]), cell, map[bool]string{true: "1", false: "0"}[tab.Nuclear[i]], tab.Features[tab.Feature[i]],
			formatF(tab.X[i]), formatF(tab.Y[i]), formatF(tab.Z[i]), formatF(tab.QV[i])}, ",") + "\n")
	}
	if err := w.WriteRowGroup([]*parquet.Values{ids, cells, nuc, feats, {Len: n, Float32: tab.X}, {Len: n, Float32: tab.Y}, {Len: n, Float32: tab.Z}, {Len: n, Float32: tab.QV}}); err != nil {
		t.Fatal(err)
	}
	w.Close()
	dir := t.TempDir()
	pq, cv := filepath.Join(dir, "transcripts.parquet"), filepath.Join(dir, "transcripts.csv")
	os.WriteFile(pq, buf.Bytes(), 0o644)
	os.WriteFile(cv, []byte(csv.String()), 0o644)
	a, err := ReadTable(pq, XeniumColumns())
	if err != nil {
		t.Fatal(err)
	}
	b, err := ReadTable(cv, XeniumColumns())
	if err != nil {
		t.Fatal(err)
	}
	for i := range n {
		if a.X[i] != b.X[i] || a.Y[i] != b.Y[i] || a.Z[i] != b.Z[i] || a.QV[i] != b.QV[i] || a.ID[i] != b.ID[i] ||
			a.Features[a.Feature[i]] != b.Features[b.Feature[i]] || a.Nuclear[i] != b.Nuclear[i] || (a.Cell[i] < 0) != (b.Cell[i] < 0) ||
			a.Cell[i] >= 0 && a.CellIDs[a.Cell[i]] != b.CellIDs[b.Cell[i]] {
			t.Fatalf("row %d differs between Parquet and CSV", i)
		}
		if a.Cell[i] >= 0 && a.CellIDs[a.Cell[i]] != tab.CellIDs[tab.Cell[i]] {
			t.Fatalf("row %d's cell %s, written %s", i, a.CellIDs[a.Cell[i]], tab.CellIDs[tab.Cell[i]])
		}
	}
}

func formatU(v uint64) string { return strconv.FormatUint(v, 10) }

// formatF writes a float32 so that it reads back to the same value.
func formatF(v float32) string { return strconv.FormatFloat(float64(v), 'g', -1, 32) }

// The inputs: the start raster gives every cell its transcripts' pixels and the pixels near them,
// the seed is the centroid of its nuclear transcripts, the priors are its nuclear transcripts.
func TestInputsFromTable(t *testing.T) {
	// Two cells side by side, 6 µm apart; a background transcript far away; a control.
	tab := &Table{Features: []string{"A", "B", "NegControlProbe_1"}, Control: []bool{false, false, true}, CellIDs: []string{"left", "right"}}
	add := func(x, y float32, f uint32, cell int32, nuclear bool, qv float32) {
		tab.X, tab.Y, tab.Feature, tab.Cell = append(tab.X, x), append(tab.Y, y), append(tab.Feature, f), append(tab.Cell, cell)
		tab.Nuclear, tab.QV = append(tab.Nuclear, nuclear), append(tab.QV, qv)
	}
	for i := range 5 {
		add(10+float32(i)*0.5, 10, 0, 0, i < 2, 40)
		add(16+float32(i)*0.5, 10, 1, 1, i >= 3, 40)
	}
	add(13, 30, 0, -1, false, 40) // the background
	add(12, 12, 2, -1, false, 40) // a control
	add(11, 11, 0, 0, false, 10)  // low quality
	in, idx, err := InputsFromTable(context.Background(), tab, StartOptions{MinQV: 20, PixelUm: 0.5, RadiusUm: 2})
	if err != nil {
		t.Fatal(err)
	}
	if in.Cells != 2 || in.Genes != 2 || len(idx.Rows) != 11 || idx.Controls != 1 || idx.LowQuality != 1 || len(in.ControlX) != 1 {
		t.Fatalf("cells %d genes %d modelled %d controls %d low %d", in.Cells, in.Genes, len(idx.Rows), idx.Controls, idx.LowQuality)
	}
	// Seeds: the nuclear transcripts' centroids.
	if math.Abs(in.SeedX[0]-10.25) > 1e-6 || math.Abs(in.SeedX[1]-17.75) > 1e-6 {
		t.Errorf("seeds %v", in.SeedX)
	}
	// The start: each cell's transcripts' pixels, and the pixels within 2 µm; the midline between
	// them; nothing near the background transcript.
	at := func(x, y float64) uint32 { return in.Start.AtMicrons(x-in.X0, y-in.Y0) }
	if at(10.2, 10.2) != 1 || at(17.2, 10.2) != 2 || at(8.4, 10.2) != 1 || at(13.4, 10.2) != 1 || at(14.6, 10.2) != 2 || at(13.1, 30.1) != 0 || at(10.2, 13.6) != 0 {
		t.Errorf("start labels at (10.2) %d, (17.2) %d, (8.4) %d, (13.4) %d, (14.6) %d, background %d, above %d", at(10.2, 10.2), at(17.2, 10.2), at(8.4, 10.2), at(13.4, 10.2), at(14.6, 10.2), at(13.1, 30.1), at(10.2, 13.6))
	}
	// Priors: nuclear transcripts only.
	priors := 0
	for i, q := range in.Transcripts.Prior {
		r := idx.Rows[i]
		if (q != 0) != tab.Nuclear[r] || q != 0 && int32(q)-1 != tab.Cell[r] {
			t.Errorf("row %d: prior %d", r, q)
		}
		if q != 0 {
			priors++
		}
	}
	if priors != 4 || idx.WithPrior != 4 || !in.Nucleus[0] || !in.Nucleus[1] {
		t.Errorf("%d priors, nuclei %v", priors, in.Nucleus)
	}
	// Without the nucleus column, every assigned transcript's cell is its prior, and no cell is
	// known to have a nucleus or not.
	tab.Nuclear = nil
	in, idx, err = InputsFromTable(context.Background(), tab, StartOptions{MinQV: 20})
	if err != nil {
		t.Fatal(err)
	}
	if idx.WithPrior != 10 || in.Nucleus != nil {
		t.Errorf("%d priors without the nucleus column, nuclei %v", idx.WithPrior, in.Nucleus)
	}
}
