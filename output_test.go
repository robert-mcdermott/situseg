package situseg

// A segmentation's files read back: every table, the matrix, the outlines and the label image
// agree with the segmentation they were written from; and the diagnostics' convergence measure.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"encoding/binary"
	"encoding/csv"
	"encoding/json"
	"io"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/robert-mcdermott/situseg/internal/parquet"
)

func readGzipCSV(t *testing.T, path string) [][]string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	zr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	rows, err := csv.NewReader(zr).ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	return rows
}

func readGzipLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	zr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	var lines []string
	sc := bufio.NewScanner(zr)
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	return lines
}

func TestFilesAreTheSegmentation(t *testing.T) {
	tab := tableOf(simulated(t, "small", 80, 60))
	spec := shortSpec(60, 40)
	spec.Options.Stages[0].Layers = 2
	seg, _ := segmentOf(t, tab, spec, 0)
	dir := t.TempDir()
	if err := seg.WriteFiles(dir, tab, WriteOptions{Labels: true}); err != nil {
		t.Fatal(err)
	}
	// The cells.
	cells := readGzipCSV(t, filepath.Join(dir, FileCells))
	if len(cells) != len(seg.Cells)+1 || cells[0][0] != "cell_id" {
		t.Fatalf("%d cell rows for %d cells", len(cells)-1, len(seg.Cells))
	}
	for c, row := range cells[1:] {
		if row[0] != seg.Cells[c].ID || row[4] != strconv.Itoa(seg.Cells[c].Transcripts) {
			t.Fatalf("cell row %d: %v", c, row)
		}
	}
	// The transcripts: every row, its cell and probability.
	tr := readGzipCSV(t, filepath.Join(dir, FileTranscripts))
	if len(tr) != tab.Rows()+1 {
		t.Fatalf("%d transcript rows for %d", len(tr)-1, tab.Rows())
	}
	head := map[string]int{}
	for i, h := range tr[0] {
		head[h] = i
	}
	for r, row := range tr[1:] {
		want := ""
		if c := seg.Assignment[r]; c >= 0 {
			want = seg.Cells[c].ID
		}
		if row[head["cell_id"]] != want || row[head["row"]] != strconv.Itoa(r) {
			t.Fatalf("transcript row %d: %v, want cell %q", r, row, want)
		}
		if p := seg.Probability[r]; math.IsNaN(float64(p)) != (row[head["probability"]] == "") {
			t.Fatalf("transcript row %d: probability %q for %g", r, row[head["probability"]], p)
		}
	}
	// The matrix: genes by cells, every count.
	mtx := readGzipLines(t, filepath.Join(dir, DirMatrix, "matrix.mtx.gz"))
	features := readGzipLines(t, filepath.Join(dir, DirMatrix, "features.tsv.gz"))
	barcodes := readGzipLines(t, filepath.Join(dir, DirMatrix, "barcodes.tsv.gz"))
	if !strings.HasPrefix(mtx[0], "%%MatrixMarket matrix coordinate integer general") || len(barcodes) != len(seg.Cells) || len(features) != len(seg.Counts.Genes) {
		t.Fatalf("matrix header %q, %d barcodes, %d features", mtx[0], len(barcodes), len(features))
	}
	k := 1
	for k < len(mtx) && strings.HasPrefix(mtx[k], "%") {
		k++
	}
	dims := strings.Fields(mtx[k])
	if dims[0] != strconv.Itoa(len(features)) || dims[1] != strconv.Itoa(len(barcodes)) || dims[2] != strconv.Itoa(len(seg.Counts.Gene)) || len(mtx)-k-1 != len(seg.Counts.Gene) {
		t.Fatalf("matrix dimensions %v, %d entries", dims, len(mtx)-k-1)
	}
	total := 0
	for _, line := range mtx[k+1:] {
		f := strings.Fields(line)
		n, _ := strconv.Atoi(f[2])
		total += n
	}
	want := 0
	for _, c := range seg.Cells {
		want += c.Transcripts
	}
	if total != want {
		t.Fatalf("the matrix counts %d transcripts, the cells %d", total, want)
	}
	// The outlines: a polygon for each cell with one, each ring closed.
	for _, name := range []string{FileBoundaries, FileLayers} {
		raw, err := os.ReadFile(filepath.Join(dir, name))
		if err != nil {
			t.Fatal(err)
		}
		var fc struct {
			Type     string `json:"type"`
			Features []struct {
				ID         string         `json:"id"`
				Properties map[string]any `json:"properties"`
				Geometry   struct {
					Type        string         `json:"type"`
					Coordinates [][][2]float64 `json:"coordinates"`
				} `json:"geometry"`
			} `json:"features"`
		}
		if err := json.Unmarshal(raw, &fc); err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if fc.Type != "FeatureCollection" || len(fc.Features) == 0 {
			t.Fatalf("%s: %s of %d features", name, fc.Type, len(fc.Features))
		}
		for _, f := range fc.Features {
			ring := f.Geometry.Coordinates[0]
			if f.Geometry.Type != "Polygon" || len(ring) < 4 || ring[0] != ring[len(ring)-1] {
				t.Fatalf("%s: feature %s's ring %v", name, f.ID, ring)
			}
		}
		if name == FileBoundaries {
			n := 0
			for _, o := range seg.Outlines {
				if len(o.X) >= 3 {
					n++
				}
			}
			if len(fc.Features) != n {
				t.Fatalf("%d polygons for %d outlines", len(fc.Features), n)
			}
		} else if _, ok := fc.Features[0].Properties["layer"]; !ok {
			t.Fatalf("a layer's feature without its layer: %v", fc.Features[0].Properties)
		}
	}
	// The label image: its size and its first strip's rows.
	raw, err := os.ReadFile(filepath.Join(dir, FileLabels))
	if err != nil {
		t.Fatal(err)
	}
	checkLabelTIFF(t, raw, seg.Labels)
	// The diagnostics.
	var diag Diagnostics
	b, _ := os.ReadFile(filepath.Join(dir, FileDiagnostics))
	if err := json.Unmarshal(b, &diag); err != nil || diag.Summary.Refined == 0 {
		t.Fatalf("diagnostics %+v, %v", diag.Summary, err)
	}
	// As Parquet.
	pdir := t.TempDir()
	if err := seg.WriteFiles(pdir, tab, WriteOptions{Parquet: true}); err != nil {
		t.Fatal(err)
	}
	for name, rows := range map[string]int{FileCellsPq: len(seg.Cells), FileTranscriptsPq: tab.Rows()} {
		f, err := os.Open(filepath.Join(pdir, name))
		if err != nil {
			t.Fatal(err)
		}
		info, _ := f.Stat()
		p, err := parquet.Open(f, info.Size())
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if int(p.NumRows) != rows {
			t.Fatalf("%s: %d rows, want %d", name, p.NumRows, rows)
		}
		if name == FileTranscriptsPq {
			v, err := p.ReadAll(p.Index("cell_id"))
			if err != nil {
				t.Fatal(err)
			}
			for r, id := range v.Bytes {
				want := ""
				if c := seg.Assignment[r]; c >= 0 {
					want = seg.Cells[c].ID
				}
				if string(id) != want {
					t.Fatalf("%s row %d: cell %q, want %q", name, r, id, want)
				}
			}
		}
		f.Close()
	}
}

// checkLabelTIFF reads a label TIFF's size and its first strip.
func checkLabelTIFF(t *testing.T, raw []byte, r *Raster) {
	t.Helper()
	le := binary.LittleEndian
	if string(raw[:4]) != "II*\x00" {
		t.Fatalf("not a little-endian TIFF: % x", raw[:4])
	}
	ifd := le.Uint32(raw[4:])
	n := int(le.Uint16(raw[ifd:]))
	tags := map[uint16]uint32{}
	for i := range n {
		e := raw[int(ifd)+2+12*i:]
		tag, typ := le.Uint16(e), le.Uint16(e[2:])
		v := le.Uint32(e[8:])
		if typ == 3 {
			v = uint32(le.Uint16(e[8:]))
		}
		tags[tag] = v
	}
	if int(tags[256]) != r.Width || int(tags[257]) != r.Height || tags[258] != 32 || tags[259] != 8 {
		t.Fatalf("TIFF tags %v for %d×%d", tags, r.Width, r.Height)
	}
	strips := (r.Height + 63) / 64
	first, count := tags[273], tags[279]
	if strips > 1 { // arrays: the first entries
		first, count = le.Uint32(raw[tags[273]:]), le.Uint32(raw[tags[279]:])
	}
	zr, err := zlib.NewReader(bytes.NewReader(raw[first : first+count]))
	if err != nil {
		t.Fatal(err)
	}
	strip, _ := io.ReadAll(zr)
	row := make([]uint32, r.Width)
	for y := 0; y < min(64, r.Height); y++ {
		r.Row(y, row)
		for x, v := range row {
			if le.Uint32(strip[4*(y*r.Width+x):]) != v {
				t.Fatalf("TIFF pixel (%d, %d) differs", x, y)
			}
		}
	}
}

func TestGewekeDiagnostic(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	noise, trend, flat := make([]float64, 500), make([]float64, 500), make([]float64, 500)
	for i := range noise {
		noise[i] = r.NormFloat64()
		trend[i] = float64(i)*0.01 + r.NormFloat64()
	}
	if z, ok := gewekeZ(noise); math.Abs(z) > 2 || !ok {
		t.Errorf("noise: z %g (%v)", z, ok)
	}
	if z, ok := gewekeZ(trend); z < 5 || !ok {
		t.Errorf("a trend: z %g (%v)", z, ok)
	}
	if z, ok := gewekeZ(flat); z != 0 || !ok {
		t.Errorf("a constant: z %g (%v)", z, ok)
	}
	// Too short to judge: under minGeweke samples.
	if z, ok := gewekeZ(noise[:minGeweke-1]); z != 0 || ok {
		t.Errorf("too short: z %g (%v)", z, ok)
	}
	if _, ok := gewekeZ(noise[:minGeweke]); !ok {
		t.Errorf("%d samples not judged", minGeweke)
	}
}
