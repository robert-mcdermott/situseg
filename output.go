package situseg

// A segmentation's files, in formats other tools read:
//
//	cells.csv.gz            each cell: its id, centroid and area (µm), transcripts, genes,
//	                        whether it was refined, its boundary certainty and its flags
//	transcripts.csv.gz      each row of the table: its cell ("" none) and the posterior
//	                        probability of it (empty: assigned by place), with its place and
//	                        feature and the cell the table gave it
//	cell_feature_matrix/    the counts in Matrix Market form, as 10x Genomics writes them
//	                        (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz): scanpy's
//	                        read_10x_mtx and Seurat's Read10X read it
//	cell_boundaries.geojson each cell's outline (µm), a GeoJSON polygon
//	cell_boundaries_layers.geojson  each layer's outlines, with more than one layer
//	diagnostics.json        convergence, the background and the cells flagged
//	cell_labels.tif         (WriteOptions.Labels) the label image, 32-bit, deflated
//
// With WriteOptions.Parquet, cells, transcripts and the boundaries (as Xenium's
// cell_boundaries.parquet: cell_id, vertex_x, vertex_y) are Parquet instead.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strconv"

	"github.com/robert-mcdermott/situseg/internal/parquet"
)

// WriteOptions say which files a segmentation writes.
type WriteOptions struct {
	// Parquet writes the cells, transcripts and boundaries as Parquet instead of CSV and GeoJSON.
	Parquet bool
	// Labels writes the label image as a TIFF.
	Labels bool
	// CreatedBy names the writer in the files that record one.
	CreatedBy string
}

// The files' names.
const (
	FileCells          = "cells.csv.gz"
	FileTranscripts    = "transcripts.csv.gz"
	FileCellsPq        = "cells.parquet"
	FileTranscriptsPq  = "transcripts.parquet"
	DirMatrix          = "cell_feature_matrix"
	FileBoundaries     = "cell_boundaries.geojson"
	FileBoundariesPq   = "cell_boundaries.parquet"
	FileLayers         = "cell_boundaries_layers.geojson"
	FileDiagnostics    = "diagnostics.json"
	FileLabels         = "cell_labels.tif"
	defaultCreatedByID = "situseg"
)

// WriteFiles writes a segmentation of a table to dir.
func (seg *Segmentation) WriteFiles(dir string, t *Table, o WriteOptions) error {
	if o.CreatedBy == "" {
		o.CreatedBy = defaultCreatedByID
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	writes := []func() error{
		func() error { return seg.writeMatrix(filepath.Join(dir, DirMatrix)) },
		func() error {
			raw, err := json.MarshalIndent(seg.Diagnostics, "", " ")
			if err != nil {
				return err
			}
			return os.WriteFile(filepath.Join(dir, FileDiagnostics), append(raw, '\n'), 0o644)
		},
	}
	if o.Parquet {
		writes = append(writes,
			func() error {
				return create(filepath.Join(dir, FileCellsPq), func(w io.Writer) error { return seg.cellsParquet(w, o.CreatedBy) })
			},
			func() error {
				return create(filepath.Join(dir, FileTranscriptsPq), func(w io.Writer) error { return seg.transcriptsParquet(w, t, o.CreatedBy) })
			},
			func() error {
				return create(filepath.Join(dir, FileBoundariesPq), func(w io.Writer) error { return seg.boundariesParquet(w, o.CreatedBy) })
			})
	} else {
		writes = append(writes,
			func() error { return createGzip(filepath.Join(dir, FileCells), gzip.DefaultCompression, seg.cellsCSV) },
			func() error {
				return createGzip(filepath.Join(dir, FileTranscripts), gzip.BestSpeed, func(w io.Writer) error { return seg.transcriptsCSV(w, t) })
			},
			func() error {
				return create(filepath.Join(dir, FileBoundaries), func(w io.Writer) error {
					return writeGeoJSON(w, seg.Cells, []geoGroup{{outlines: seg.Outlines}})
				})
			})
	}
	if len(seg.Layers) > 1 {
		writes = append(writes, func() error {
			return create(filepath.Join(dir, FileLayers), func(w io.Writer) error {
				var groups []geoGroup
				for z, l := range seg.Layers {
					groups = append(groups, geoGroup{outlines: l.Outlines, props: map[string]any{"layer": z, "z0": l.Z0, "z1": l.Z1}})
				}
				return writeGeoJSON(w, seg.Cells, groups)
			})
		})
	}
	if o.Labels {
		writes = append(writes, func() error {
			return create(filepath.Join(dir, FileLabels), func(w io.Writer) error { return writeLabelTIFF(w, seg.Labels) })
		})
	}
	for _, write := range writes {
		if err := write(); err != nil {
			if f := Classify(err); f != nil {
				return f
			}
			return err
		}
	}
	return nil
}

// create writes a file through a buffer.
func create(path string, write func(io.Writer) error) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	bw := bufio.NewWriterSize(f, 1<<20)
	if err := write(bw); err != nil {
		f.Close()
		return fmt.Errorf("%s: %w", filepath.Base(path), err)
	}
	if err := bw.Flush(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// createGzip writes a gzipped file.
func createGzip(path string, level int, write func(io.Writer) error) error {
	return create(path, func(w io.Writer) error {
		zw, err := gzip.NewWriterLevel(w, level)
		if err != nil {
			return err
		}
		if err := write(zw); err != nil {
			return err
		}
		return zw.Close()
	})
}

// csvField quotes a text field when it needs it.
func csvField(b []byte, s string) []byte {
	need := false
	for i := 0; i < len(s); i++ {
		if c := s[i]; c == ',' || c == '"' || c == '\n' || c == '\r' {
			need = true
			break
		}
	}
	if !need {
		return append(b, s...)
	}
	b = append(b, '"')
	for i := 0; i < len(s); i++ {
		if s[i] == '"' {
			b = append(b, '"')
		}
		b = append(b, s[i])
	}
	return append(b, '"')
}

// appendFloat writes a number, empty for NaN.
func appendFloat(b []byte, v float64, prec int) []byte {
	if math.IsNaN(v) {
		return b
	}
	return strconv.AppendFloat(b, v, 'f', prec, 64)
}

func (seg *Segmentation) nucleusOf(c int) string {
	if seg.HasNucleus == nil {
		return ""
	}
	return strconv.FormatBool(seg.HasNucleus[c])
}

// flags are the cell's limitations: few transcripts, and uncertain boundaries.
func (seg *Segmentation) flags(c int) (low, uncertain bool) {
	cell := &seg.Cells[c]
	return cell.Refined && cell.Transcripts < lowTranscripts, cell.Refined && !(cell.Certainty >= uncertainShare)
}

func (seg *Segmentation) cellsCSV(w io.Writer) error {
	b := []byte("cell_id,x_centroid,y_centroid,area_um2,transcripts,genes,assigned,refined,boundary_certainty,has_nucleus,low_transcripts,uncertain_boundary\n")
	for c := range seg.Cells {
		cell := &seg.Cells[c]
		b = csvField(b, cell.ID)
		b = append(b, ',')
		b = appendFloat(b, cell.X, 3)
		b = append(b, ',')
		b = appendFloat(b, cell.Y, 3)
		b = append(b, ',')
		b = appendFloat(b, cell.AreaUm2, 3)
		b = fmt.Appendf(b, ",%d,%d,%d,%v,", cell.Transcripts, cell.Genes, cell.Assigned, cell.Refined)
		b = appendFloat(b, cell.Certainty, 4)
		low, uncertain := seg.flags(c)
		b = fmt.Appendf(b, ",%s,%v,%v\n", seg.nucleusOf(c), low, uncertain)
		if len(b) > 1<<16 {
			if _, err := w.Write(b); err != nil {
				return err
			}
			b = b[:0]
		}
	}
	_, err := w.Write(b)
	return err
}

func (seg *Segmentation) transcriptsCSV(w io.Writer, t *Table) error {
	b := []byte("row")
	if t.ID != nil {
		b = append(b, ",transcript_id"...)
	}
	b = append(b, ",feature_name,x,y"...)
	if t.Z != nil {
		b = append(b, ",z"...)
	}
	if t.QV != nil {
		b = append(b, ",qv"...)
	}
	b = append(b, ",cell_id,probability,prior_cell_id\n"...)
	for r := range seg.Assignment {
		b = strconv.AppendInt(b, int64(r), 10)
		if t.ID != nil {
			b = append(b, ',')
			b = strconv.AppendUint(b, t.ID[r], 10)
		}
		b = append(b, ',')
		b = csvField(b, t.Features[t.Feature[r]])
		b = append(b, ',')
		b = appendFloat(b, float64(t.X[r]), 3)
		b = append(b, ',')
		b = appendFloat(b, float64(t.Y[r]), 3)
		if t.Z != nil {
			b = append(b, ',')
			b = appendFloat(b, float64(t.Z[r]), 3)
		}
		if t.QV != nil {
			b = append(b, ',')
			b = appendFloat(b, float64(t.QV[r]), 2)
		}
		b = append(b, ',')
		if c := seg.Assignment[r]; c >= 0 {
			b = csvField(b, seg.Cells[c].ID)
		}
		b = append(b, ',')
		b = appendFloat(b, float64(seg.Probability[r]), 4)
		b = append(b, ',')
		if c := t.Cell[r]; c >= 0 {
			b = csvField(b, t.CellIDs[c])
		}
		b = append(b, '\n')
		if len(b) > 1<<16 {
			if _, err := w.Write(b); err != nil {
				return err
			}
			b = b[:0]
		}
	}
	_, err := w.Write(b)
	return err
}

// writeMatrix writes the counts as 10x Genomics' Matrix Market folder: genes by cells.
func (seg *Segmentation) writeMatrix(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	k := seg.Counts
	err := createGzip(filepath.Join(dir, "features.tsv.gz"), gzip.DefaultCompression, func(w io.Writer) error {
		var b []byte
		for _, g := range k.Genes {
			b = fmt.Appendf(b, "%s\t%s\tGene Expression\n", g, g)
		}
		_, err := w.Write(b)
		return err
	})
	if err != nil {
		return err
	}
	err = createGzip(filepath.Join(dir, "barcodes.tsv.gz"), gzip.DefaultCompression, func(w io.Writer) error {
		var b []byte
		for _, c := range seg.Cells {
			b = append(append(b, c.ID...), '\n')
		}
		_, err := w.Write(b)
		return err
	})
	if err != nil {
		return err
	}
	return createGzip(filepath.Join(dir, "matrix.mtx.gz"), gzip.DefaultCompression, func(w io.Writer) error {
		b := fmt.Appendf(nil, "%%%%MatrixMarket matrix coordinate integer general\n%%metadata_json: {\"software_version\": \"situseg\", \"format_version\": 2}\n%d %d %d\n",
			len(k.Genes), len(seg.Cells), len(k.Gene))
		for c := range seg.Cells {
			for i := k.Offsets[c]; i < k.Offsets[c+1]; i++ {
				b = fmt.Appendf(b, "%d %d %d\n", k.Gene[i]+1, c+1, k.Count[i])
			}
			if len(b) > 1<<16 {
				if _, err := w.Write(b); err != nil {
					return err
				}
				b = b[:0]
			}
		}
		_, err := w.Write(b)
		return err
	})
}

// ring is an outline's closed ring of [x, y] pairs (µm), with the first vertex repeated.
func appendRing(b []byte, o *Outline) []byte {
	b = append(b, '[')
	for i := 0; i <= len(o.X); i++ {
		j := i % len(o.X)
		if i > 0 {
			b = append(b, ',')
		}
		b = append(b, '[')
		b = strconv.AppendFloat(b, float64(o.X[j]), 'f', 3, 32)
		b = append(b, ',')
		b = strconv.AppendFloat(b, float64(o.Y[j]), 'f', 3, 32)
		b = append(b, ']')
	}
	return append(b, ']')
}

// geoGroup is a set of outlines written as GeoJSON features, with properties of their own.
type geoGroup struct {
	outlines []Outline
	props    map[string]any
}

// writeGeoJSON writes outlines as a FeatureCollection of polygons, one per cell with an outline in
// each group.
func writeGeoJSON(w io.Writer, cells []Cell, groups []geoGroup) error {
	b := []byte(`{"type":"FeatureCollection","features":[`)
	first := true
	for _, g := range groups {
		var props []byte
		for _, k := range []string{"layer", "z0", "z1"} {
			if v, ok := g.props[k]; ok {
				raw, _ := json.Marshal(v)
				props = fmt.Appendf(props, `,%q:%s`, k, raw)
			}
		}
		for c := range cells {
			o := &g.outlines[c]
			if len(o.X) < 3 {
				continue
			}
			if !first {
				b = append(b, ',')
			}
			first = false
			id, _ := json.Marshal(cells[c].ID)
			b = fmt.Appendf(b, `{"type":"Feature","id":%s,"properties":{"cell_id":%s%s},"geometry":{"type":"Polygon","coordinates":[`, id, id, props)
			b = appendRing(b, o)
			b = append(b, "]}}\n"...)
			if len(b) > 1<<16 {
				if _, err := w.Write(b); err != nil {
					return err
				}
				b = b[:0]
			}
		}
	}
	b = append(b, "]}\n"...)
	_, err := w.Write(b)
	return err
}

func (seg *Segmentation) cellsParquet(w io.Writer, createdBy string) error {
	n := len(seg.Cells)
	ids := &parquet.Values{Len: n, Bytes: make([][]byte, n)}
	xs, ys, area, cert := &parquet.Values{Len: n, Float64: make([]float64, n)}, &parquet.Values{Len: n, Float64: make([]float64, n)},
		&parquet.Values{Len: n, Float64: make([]float64, n)}, &parquet.Values{Len: n, Float64: make([]float64, n)}
	tr, genes, assigned := &parquet.Values{Len: n, Int64: make([]int64, n)}, &parquet.Values{Len: n, Int64: make([]int64, n)}, &parquet.Values{Len: n, Int64: make([]int64, n)}
	refined, low, uncertain := &parquet.Values{Len: n, Bool: make([]bool, n)}, &parquet.Values{Len: n, Bool: make([]bool, n)}, &parquet.Values{Len: n, Bool: make([]bool, n)}
	for c, cell := range seg.Cells {
		ids.Bytes[c] = []byte(cell.ID)
		xs.Float64[c], ys.Float64[c], area.Float64[c], cert.Float64[c] = cell.X, cell.Y, cell.AreaUm2, cell.Certainty
		tr.Int64[c], genes.Int64[c], assigned.Int64[c] = int64(cell.Transcripts), int64(cell.Genes), int64(cell.Assigned)
		refined.Bool[c] = cell.Refined
		low.Bool[c], uncertain.Bool[c] = seg.flags(c)
	}
	fields := []parquet.Field{{Name: "cell_id", Type: parquet.ByteArray, Annotation: "string"}, {Name: "x_centroid", Type: parquet.Double},
		{Name: "y_centroid", Type: parquet.Double}, {Name: "area_um2", Type: parquet.Double}, {Name: "transcripts", Type: parquet.Int64},
		{Name: "genes", Type: parquet.Int64}, {Name: "assigned", Type: parquet.Int64}, {Name: "refined", Type: parquet.Boolean},
		{Name: "boundary_certainty", Type: parquet.Double}, {Name: "low_transcripts", Type: parquet.Boolean}, {Name: "uncertain_boundary", Type: parquet.Boolean}}
	cols := []*parquet.Values{ids, xs, ys, area, tr, genes, assigned, refined, cert, low, uncertain}
	if seg.HasNucleus != nil {
		nuc := &parquet.Values{Len: n, Bool: append([]bool(nil), seg.HasNucleus...)}
		fields = append(fields, parquet.Field{Name: "has_nucleus", Type: parquet.Boolean})
		cols = append(cols, nuc)
	}
	return writeParquet(w, fields, cols, n, createdBy)
}

func (seg *Segmentation) transcriptsParquet(w io.Writer, t *Table, createdBy string) error {
	n := len(seg.Assignment)
	fields := []parquet.Field{{Name: "row", Type: parquet.Int64}}
	if t.ID != nil {
		fields = append(fields, parquet.Field{Name: "transcript_id", Type: parquet.Int64, Annotation: "uint"})
	}
	fields = append(fields, parquet.Field{Name: "feature_name", Type: parquet.ByteArray, Annotation: "string"},
		parquet.Field{Name: "x", Type: parquet.Float}, parquet.Field{Name: "y", Type: parquet.Float})
	if t.Z != nil {
		fields = append(fields, parquet.Field{Name: "z", Type: parquet.Float})
	}
	if t.QV != nil {
		fields = append(fields, parquet.Field{Name: "qv", Type: parquet.Float})
	}
	fields = append(fields, parquet.Field{Name: "cell_id", Type: parquet.ByteArray, Annotation: "string"},
		parquet.Field{Name: "probability", Type: parquet.Float}, parquet.Field{Name: "prior_cell_id", Type: parquet.ByteArray, Annotation: "string"})
	pw, err := parquet.NewWriter(w, fields, createdBy)
	if err != nil {
		return err
	}
	const group = 1 << 20
	for lo := 0; lo < n; lo += group {
		hi := min(n, lo+group)
		m := hi - lo
		rows := &parquet.Values{Len: m, Int64: make([]int64, m)}
		feat, cell, prior := &parquet.Values{Len: m, Bytes: make([][]byte, m)}, &parquet.Values{Len: m, Bytes: make([][]byte, m)}, &parquet.Values{Len: m, Bytes: make([][]byte, m)}
		cols := []*parquet.Values{rows}
		if t.ID != nil {
			ids := &parquet.Values{Len: m, Int64: make([]int64, m)}
			for i := range m {
				ids.Int64[i] = int64(t.ID[lo+i])
			}
			cols = append(cols, ids)
		}
		features := make([][]byte, len(t.Features))
		for f, name := range t.Features {
			features[f] = []byte(name)
		}
		cellIDs := make([][]byte, len(seg.Cells))
		for c := range seg.Cells {
			cellIDs[c] = []byte(seg.Cells[c].ID)
		}
		empty := []byte{}
		for i := range m {
			r := lo + i
			rows.Int64[i] = int64(r)
			feat.Bytes[i] = features[t.Feature[r]]
			cell.Bytes[i], prior.Bytes[i] = empty, empty
			if c := seg.Assignment[r]; c >= 0 {
				cell.Bytes[i] = cellIDs[c]
			}
			if c := t.Cell[r]; c >= 0 {
				prior.Bytes[i] = []byte(t.CellIDs[c])
			}
		}
		cols = append(cols, feat, &parquet.Values{Len: m, Float32: t.X[lo:hi]}, &parquet.Values{Len: m, Float32: t.Y[lo:hi]})
		if t.Z != nil {
			cols = append(cols, &parquet.Values{Len: m, Float32: t.Z[lo:hi]})
		}
		if t.QV != nil {
			cols = append(cols, &parquet.Values{Len: m, Float32: t.QV[lo:hi]})
		}
		cols = append(cols, cell, &parquet.Values{Len: m, Float32: seg.Probability[lo:hi]}, prior)
		if err := pw.WriteRowGroup(cols); err != nil {
			return err
		}
	}
	return pw.Close()
}

func (seg *Segmentation) boundariesParquet(w io.Writer, createdBy string) error {
	var ids [][]byte
	var xs, ys []float32
	for c := range seg.Cells {
		o := &seg.Outlines[c]
		if len(o.X) < 3 {
			continue
		}
		id := []byte(seg.Cells[c].ID)
		for i := 0; i <= len(o.X); i++ {
			j := i % len(o.X)
			ids, xs, ys = append(ids, id), append(xs, o.X[j]), append(ys, o.Y[j])
		}
	}
	n := len(ids)
	fields := []parquet.Field{{Name: "cell_id", Type: parquet.ByteArray, Annotation: "string"}, {Name: "vertex_x", Type: parquet.Float}, {Name: "vertex_y", Type: parquet.Float}}
	return writeParquet(w, fields, []*parquet.Values{{Len: n, Bytes: ids}, {Len: n, Float32: xs}, {Len: n, Float32: ys}}, n, createdBy)
}

// writeParquet writes columns as one Parquet row group.
func writeParquet(w io.Writer, fields []parquet.Field, cols []*parquet.Values, n int, createdBy string) error {
	pw, err := parquet.NewWriter(w, fields, createdBy)
	if err != nil {
		return err
	}
	if n > 0 {
		if err := pw.WriteRowGroup(cols); err != nil {
			return err
		}
	}
	return pw.Close()
}

// writeLabelTIFF writes a label raster as a baseline TIFF of 32-bit unsigned samples, in strips of
// 64 rows, each deflated (compression 8), with the pixel size as its resolution (pixels per cm).
func writeLabelTIFF(w io.Writer, r *Raster) error {
	const rowsPerStrip = 64
	strips := (r.Height + rowsPerStrip - 1) / rowsPerStrip
	var data []byte
	offsets, counts := make([]uint32, strips), make([]uint32, strips)
	row := make([]uint32, r.Width)
	raw := make([]byte, 4*r.Width)
	const header = 8
	for s := range strips {
		var buf bytes.Buffer
		zw := zlib.NewWriter(&buf)
		for y := s * rowsPerStrip; y < min(r.Height, (s+1)*rowsPerStrip); y++ {
			r.Row(y, row)
			for x, v := range row {
				binary.LittleEndian.PutUint32(raw[4*x:], v)
			}
			zw.Write(raw)
		}
		zw.Close()
		if uint64(header+len(data)+buf.Len()) > math.MaxUint32-1<<20 {
			return fmt.Errorf("the label image compresses to more than 4 GB, more than a TIFF holds")
		}
		offsets[s], counts[s] = uint32(header+len(data)), uint32(buf.Len())
		data = append(data, buf.Bytes()...)
	}
	if len(data)%2 == 1 {
		data = append(data, 0)
	}
	// The IFD, after the strips: its entries, then the arrays and rationals they point to.
	ifd := header + len(data)
	type entry struct {
		tag, typ uint16
		count    uint32
		value    []byte // inline when it fits 4 bytes, else stored after the IFD
	}
	le := binary.LittleEndian
	u16 := func(v uint16) []byte { return le.AppendUint16(nil, v) }
	u32s := func(v []uint32) []byte {
		var b []byte
		for _, x := range v {
			b = le.AppendUint32(b, x)
		}
		return b
	}
	perCm := uint32(math.Round(1e4 / r.PixelSize * 1000)) // pixels per cm, over 1000
	rational := u32s([]uint32{perCm, 1000})
	entries := []entry{
		{256, 4, 1, u32s([]uint32{uint32(r.Width)})}, // ImageWidth
		{257, 4, 1, u32s([]uint32{uint32(r.Height)})},
		{258, 3, 1, u16(32)}, // BitsPerSample
		{259, 3, 1, u16(8)},  // Compression: Adobe deflate
		{262, 3, 1, u16(1)},  // PhotometricInterpretation: black is zero
		{273, 4, uint32(strips), u32s(offsets)},
		{277, 3, 1, u16(1)}, // SamplesPerPixel
		{278, 4, 1, u32s([]uint32{rowsPerStrip})},
		{279, 4, uint32(strips), u32s(counts)},
		{282, 5, 1, rational}, // XResolution
		{283, 5, 1, rational},
		{296, 3, 1, u16(3)}, // ResolutionUnit: cm
		{339, 3, 1, u16(1)}, // SampleFormat: unsigned
	}
	extra := ifd + 2 + 12*len(entries) + 4
	var tail []byte
	var dir []byte
	dir = le.AppendUint16(dir, uint16(len(entries)))
	for _, e := range entries {
		dir = le.AppendUint16(dir, e.tag)
		dir = le.AppendUint16(dir, e.typ)
		dir = le.AppendUint32(dir, e.count)
		if len(e.value) <= 4 {
			v := append(append([]byte(nil), e.value...), make([]byte, 4-len(e.value))...)
			dir = append(dir, v...)
			continue
		}
		dir = le.AppendUint32(dir, uint32(extra+len(tail)))
		tail = append(tail, e.value...)
	}
	dir = le.AppendUint32(dir, 0) // no next IFD
	head := []byte{'I', 'I', 42, 0}
	head = le.AppendUint32(head, uint32(ifd))
	for _, part := range [][]byte{head, data, dir, tail} {
		if _, err := w.Write(part); err != nil {
			return err
		}
	}
	return nil
}
