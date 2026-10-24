package situseg

// Reading a table of transcripts from Parquet or CSV (comma- or tab-separated, gzipped or not),
// with the columns an instrument names: Xenium's transcripts.parquet or transcripts.csv.gz, CosMx's
// tx_file.csv, MERSCOPE's detected_transcripts.csv, or any table with a place, a feature and a cell
// for each molecule.

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/robert-mcdermott/situseg/internal/parquet"
)

// Columns say which of a table's columns hold what, and how to read them.
type Columns struct {
	// X, Y, Feature and Cell are required; Z, QV, Nuclear and ID may be "" (the table has none).
	X, Y, Z string
	Feature string
	QV      string
	Cell    string
	Nuclear string
	ID      string
	// CellGroup, when set, is a column whose value makes cell ids unique with the cell's (CosMx
	// numbers cells within each field of view): a cell's id is "<group>_<cell>".
	CellGroup string
	// NoCell are the cell column's values that mean no cell.
	NoCell []string
	// NuclearValues are a text nucleus column's values that mean the molecule lies over its cell's
	// nucleus (compared without case); in a numeric column any value but 0 does, in a boolean one
	// true.
	NuclearValues []string
	// Controls is a regular expression matching the features that are negative controls.
	Controls string
	// Scale and ZScale turn the coordinates into µm (0: 1).
	Scale, ZScale float64
}

// XeniumColumns are 10x Genomics Xenium's transcripts.parquet's (and transcripts.csv.gz's).
func XeniumColumns() Columns {
	return Columns{X: "x_location", Y: "y_location", Z: "z_location", Feature: "feature_name", QV: "qv",
		Cell: "cell_id", Nuclear: "overlaps_nucleus", ID: "transcript_id",
		NoCell: []string{"UNASSIGNED", "-1", ""}, NuclearValues: []string{"1", "true"},
		Controls: `^(NegControlProbe_|NegControlCodeword_|UnassignedCodeword_|BLANK_|DeprecatedCodeword_|Intergenic_Region|GenomicControl)`}
}

// CosMxColumns are NanoString CosMx's tx_file.csv's: places in pixels of 0.12028 µm, cells
// numbered within each field of view (0: none), and the compartment each molecule lies in.
func CosMxColumns() Columns {
	return Columns{X: "x_global_px", Y: "y_global_px", Z: "z", Feature: "target", Cell: "cell_ID", CellGroup: "fov",
		Nuclear: "CellComp", NoCell: []string{"0", ""}, NuclearValues: []string{"Nuclear"},
		Controls: `^(NegPrb|Negative|SystemControl|FalseCode)`, Scale: 0.12028}
}

// MERSCOPEColumns are Vizgen MERSCOPE's detected_transcripts.csv's: depths as z-plane numbers, 1.5
// µm apart, and no nucleus column.
func MERSCOPEColumns() Columns {
	return Columns{X: "global_x", Y: "global_y", Z: "global_z", Feature: "gene", Cell: "cell_id",
		NoCell: []string{"-1", ""}, Controls: `^Blank-`, ZScale: 1.5}
}

// tableBuilder gathers a table's rows.
type tableBuilder struct {
	c        Columns
	t        *Table
	features map[string]uint32
	cells    map[string]int32
	noCell   map[string]bool
	nuclear  map[string]bool
	controls *regexp.Regexp
	scale    float32
	zscale   float32
	key      []byte
}

func newTableBuilder(c Columns) (*tableBuilder, error) {
	if c.X == "" || c.Y == "" || c.Feature == "" || c.Cell == "" {
		return nil, errors.New("situseg: the columns of x, y, the feature and the cell are needed")
	}
	b := &tableBuilder{c: c, t: &Table{}, features: map[string]uint32{}, cells: map[string]int32{},
		noCell: map[string]bool{}, nuclear: map[string]bool{}, scale: 1, zscale: 1}
	for _, v := range c.NoCell {
		b.noCell[v] = true
	}
	for _, v := range c.NuclearValues {
		b.nuclear[strings.ToLower(v)] = true
	}
	if c.Controls != "" {
		re, err := regexp.Compile(c.Controls)
		if err != nil {
			return nil, fmt.Errorf("situseg: the controls' pattern: %w", err)
		}
		b.controls = re
	}
	if c.Scale != 0 {
		b.scale = float32(c.Scale)
	}
	if c.ZScale != 0 {
		b.zscale = float32(c.ZScale)
	}
	if c.Z != "" {
		b.t.Z = []float32{}
	}
	if c.QV != "" {
		b.t.QV = []float32{}
	}
	if c.Nuclear != "" {
		b.t.Nuclear = []bool{}
	}
	if c.ID != "" {
		b.t.ID = []uint64{}
	}
	return b, nil
}

func (b *tableBuilder) feature(name []byte) uint32 {
	if f, ok := b.features[string(name)]; ok {
		return f
	}
	f := uint32(len(b.t.Features))
	s := string(name)
	b.features[s] = f
	b.t.Features = append(b.t.Features, s)
	b.t.Control = append(b.t.Control, b.controls != nil && b.controls.MatchString(s))
	return f
}

// cell interns a cell id (in its group, when the table has one); −1 for none.
func (b *tableBuilder) cell(group, id []byte) int32 {
	if b.noCell[string(id)] {
		return -1
	}
	key := id
	if b.c.CellGroup != "" {
		b.key = append(append(append(b.key[:0], group...), '_'), id...)
		key = b.key
	}
	if c, ok := b.cells[string(key)]; ok {
		return c
	}
	c := int32(len(b.t.CellIDs))
	s := string(key)
	b.cells[s] = c
	b.t.CellIDs = append(b.t.CellIDs, s)
	return c
}

func (b *tableBuilder) nuclearText(v []byte) bool { return b.nuclear[strings.ToLower(string(v))] }

// finish checks what was read.
func (b *tableBuilder) finish(name string) (*Table, error) {
	t := b.t
	if t.Rows() == 0 {
		return nil, fmt.Errorf("%s has no rows", name)
	}
	if len(t.CellIDs) == 0 {
		return nil, fail(ClassPrior, remedyPrior, "%s: no row is in a cell (column %q)", name, b.c.Cell)
	}
	return t, t.Check()
}

// ReadTable reads a table of transcripts: Parquet (.parquet) or CSV (.csv, comma-separated, or
// .tsv and .txt, tab-separated; each may be gzipped, .gz).
func ReadTable(path string, c Columns) (*Table, error) {
	name := strings.ToLower(filepath.Base(path))
	if strings.HasSuffix(name, ".parquet") {
		return readParquetTable(path, c)
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var r io.Reader = bufio.NewReaderSize(f, 1<<20)
	if strings.HasSuffix(name, ".gz") {
		zr, err := gzip.NewReader(r)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		defer zr.Close()
		r = bufio.NewReaderSize(zr, 1<<20)
		name = strings.TrimSuffix(name, ".gz")
	}
	sep := ','
	if strings.HasSuffix(name, ".tsv") || strings.HasSuffix(name, ".txt") {
		sep = '\t'
	}
	return readCSVTable(r, sep, path, c)
}

func readCSVTable(r io.Reader, sep rune, path string, c Columns) (*Table, error) {
	b, err := newTableBuilder(c)
	if err != nil {
		return nil, err
	}
	cr := csv.NewReader(r)
	cr.Comma = sep
	cr.ReuseRecord = true
	cr.FieldsPerRecord = -1
	head, err := cr.Read()
	if err != nil {
		return nil, fmt.Errorf("%s: the header: %w", path, err)
	}
	index := map[string]int{}
	for i, h := range head {
		if i == 0 {
			h = strings.TrimPrefix(h, "\uFEFF")
		}
		index[strings.TrimSpace(h)] = i
	}
	col := func(name string, need bool) (int, error) {
		if name == "" {
			return -1, nil
		}
		i, ok := index[name]
		if !ok {
			if need {
				return -1, fail(ClassInput, remedyInput, "%s has no column %q (it has %s)", path, name, strings.Join(head, ", "))
			}
			return -1, nil
		}
		return i, nil
	}
	var ix, iy, iz, ifeat, iqv, icell, inuc, iid, igroup int
	for _, s := range []struct {
		name string
		need bool
		dst  *int
	}{{c.X, true, &ix}, {c.Y, true, &iy}, {c.Feature, true, &ifeat}, {c.Cell, true, &icell}, {c.Z, false, &iz},
		{c.QV, false, &iqv}, {c.Nuclear, false, &inuc}, {c.ID, false, &iid}, {c.CellGroup, c.CellGroup != "", &igroup}} {
		if *s.dst, err = col(s.name, s.need); err != nil {
			return nil, err
		}
	}
	t := b.t
	if iz < 0 {
		t.Z = nil
	}
	if iqv < 0 {
		t.QV = nil
	}
	if inuc < 0 {
		t.Nuclear = nil
	}
	if iid < 0 {
		t.ID = nil
	}
	num := func(rec []string, i int, line int, what string) (float32, error) {
		s := strings.TrimSpace(rec[i])
		if s == "" || strings.EqualFold(s, "nan") || strings.EqualFold(s, "na") {
			return float32(math.NaN()), nil
		}
		v, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return 0, fail(ClassInput, remedyInput, "%s line %d: %s %q is not a number", path, line, what, s)
		}
		return float32(v), nil
	}
	maxIndex := max(ix, iy, iz, ifeat, iqv, icell, inuc, iid, igroup)
	for line := 2; ; line++ {
		rec, err := cr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		if len(rec) <= maxIndex {
			return nil, fail(ClassInput, remedyInput, "%s line %d has %d fields, the header %d", path, line, len(rec), len(head))
		}
		x, err := num(rec, ix, line, c.X)
		if err != nil {
			return nil, err
		}
		y, err := num(rec, iy, line, c.Y)
		if err != nil {
			return nil, err
		}
		t.X, t.Y = append(t.X, x*b.scale), append(t.Y, y*b.scale)
		if iz >= 0 {
			z, err := num(rec, iz, line, c.Z)
			if err != nil {
				return nil, err
			}
			t.Z = append(t.Z, z*b.zscale)
		}
		if iqv >= 0 {
			q, err := num(rec, iqv, line, c.QV)
			if err != nil {
				return nil, err
			}
			t.QV = append(t.QV, q)
		}
		t.Feature = append(t.Feature, b.feature([]byte(rec[ifeat])))
		var group []byte
		if igroup >= 0 {
			group = []byte(strings.TrimSpace(rec[igroup]))
		}
		t.Cell = append(t.Cell, b.cell(group, []byte(strings.TrimSpace(rec[icell]))))
		if inuc >= 0 {
			v := strings.TrimSpace(rec[inuc])
			nuc := b.nuclearText([]byte(v))
			if f, err := strconv.ParseFloat(v, 64); err == nil {
				nuc = f != 0
			}
			t.Nuclear = append(t.Nuclear, nuc)
		}
		if t.ID != nil {
			id, err := strconv.ParseUint(strings.TrimSpace(rec[iid]), 10, 64)
			if err != nil {
				t.ID = nil // not numeric: rows are known by their number instead
			} else {
				t.ID = append(t.ID, id)
			}
		}
	}
	return b.finish(path)
}

func readParquetTable(path string, c Columns) (*Table, error) {
	b, err := newTableBuilder(c)
	if err != nil {
		return nil, err
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	st, err := f.Stat()
	if err != nil {
		return nil, err
	}
	p, err := parquet.Open(f, st.Size())
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	col := func(name string, need bool) (int, error) {
		if name == "" {
			return -1, nil
		}
		i := p.Index(name)
		if i < 0 && need {
			var names []string
			for _, c := range p.Columns {
				names = append(names, c.Name)
			}
			return -1, fail(ClassInput, remedyInput, "%s has no column %q (it has %s)", path, name, strings.Join(names, ", "))
		}
		return i, nil
	}
	var ix, iy, iz, ifeat, iqv, icell, inuc, iid, igroup int
	for _, s := range []struct {
		name string
		need bool
		dst  *int
	}{{c.X, true, &ix}, {c.Y, true, &iy}, {c.Feature, true, &ifeat}, {c.Cell, true, &icell}, {c.Z, false, &iz},
		{c.QV, false, &iqv}, {c.Nuclear, false, &inuc}, {c.ID, false, &iid}, {c.CellGroup, c.CellGroup != "", &igroup}} {
		if *s.dst, err = col(s.name, s.need); err != nil {
			return nil, err
		}
	}
	t := b.t
	if iz < 0 {
		t.Z = nil
	}
	if iqv < 0 {
		t.QV = nil
	}
	if inuc < 0 {
		t.Nuclear = nil
	}
	if iid < 0 {
		t.ID = nil
	}
	n := int(p.NumRows)
	t.X, t.Y, t.Feature, t.Cell = make([]float32, 0, n), make([]float32, 0, n), make([]uint32, 0, n), make([]int32, 0, n)
	floats := func(v *parquet.Values, scale float32, dst []float32) ([]float32, error) {
		null := float32(math.NaN())
		for i := 0; i < v.Len; i++ {
			var x float32
			switch {
			case v.Valid != nil && !v.Valid[i]:
				dst = append(dst, null)
				continue
			case v.Float32 != nil:
				x = v.Float32[i]
			case v.Float64 != nil:
				x = float32(v.Float64[i])
			case v.Int32 != nil:
				x = float32(v.Int32[i])
			case v.Int64 != nil:
				x = float32(v.Int64[i])
			default:
				return nil, fmt.Errorf("%s: a coordinate or quality column is not numeric", path)
			}
			dst = append(dst, x*scale)
		}
		return dst, nil
	}
	text := func(v *parquet.Values, i int, buf []byte) []byte {
		switch {
		case v.Valid != nil && !v.Valid[i]:
			return buf[:0]
		case v.Bytes != nil:
			return v.Bytes[i]
		case v.Int64 != nil:
			return strconv.AppendInt(buf[:0], v.Int64[i], 10)
		case v.Int32 != nil:
			return strconv.AppendInt(buf[:0], int64(v.Int32[i]), 10)
		case v.Bool != nil:
			return strconv.AppendBool(buf[:0], v.Bool[i])
		}
		return buf[:0]
	}
	var buf, gbuf []byte
	for g := range p.RowGroups {
		read := func(i int) (*parquet.Values, error) {
			if i < 0 {
				return nil, nil
			}
			return p.ReadColumn(g, i)
		}
		var cols [9]*parquet.Values
		for k, i := range []int{ix, iy, iz, ifeat, iqv, icell, inuc, iid, igroup} {
			if cols[k], err = read(i); err != nil {
				return nil, fmt.Errorf("%s: %w", path, err)
			}
		}
		vx, vy, vz, vf, vq, vc, vn, vid, vg := cols[0], cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
		if t.X, err = floats(vx, b.scale, t.X); err != nil {
			return nil, err
		}
		if t.Y, err = floats(vy, b.scale, t.Y); err != nil {
			return nil, err
		}
		if vz != nil {
			if t.Z, err = floats(vz, b.zscale, t.Z); err != nil {
				return nil, err
			}
		}
		if vq != nil {
			if t.QV, err = floats(vq, 1, t.QV); err != nil {
				return nil, err
			}
		}
		if vf.Bytes == nil {
			return nil, fail(ClassInput, remedyInput, "%s: the feature column %q is not text", path, c.Feature)
		}
		for i := 0; i < vf.Len; i++ {
			t.Feature = append(t.Feature, b.feature(text(vf, i, buf)))
		}
		for i := 0; i < vc.Len; i++ {
			buf = text(vc, i, buf)
			var group []byte
			if vg != nil {
				gbuf = append(gbuf[:0], text(vg, i, nil)...)
				group = gbuf
			}
			t.Cell = append(t.Cell, b.cell(group, buf))
		}
		if vn != nil {
			for i := 0; i < vn.Len; i++ {
				var nuc bool
				switch {
				case vn.Valid != nil && !vn.Valid[i]:
				case vn.Bool != nil:
					nuc = vn.Bool[i]
				case vn.Int32 != nil:
					nuc = vn.Int32[i] != 0
				case vn.Int64 != nil:
					nuc = vn.Int64[i] != 0
				case vn.Float32 != nil:
					nuc = vn.Float32[i] != 0
				case vn.Float64 != nil:
					nuc = vn.Float64[i] != 0
				case vn.Bytes != nil:
					nuc = b.nuclearText(vn.Bytes[i])
				}
				t.Nuclear = append(t.Nuclear, nuc)
			}
		}
		if vid != nil && t.ID != nil {
			switch {
			case vid.Int64 != nil:
				for _, id := range vid.Int64 {
					t.ID = append(t.ID, uint64(id))
				}
			case vid.Int32 != nil:
				for _, id := range vid.Int32 {
					t.ID = append(t.ID, uint64(uint32(id)))
				}
			default:
				t.ID = nil
			}
		}
	}
	return b.finish(path)
}

// TableHeader is a table's column names: a Parquet file's schema, a CSV file's first line.
func TableHeader(path string) ([]string, error) {
	name := strings.ToLower(filepath.Base(path))
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	if strings.HasSuffix(name, ".parquet") {
		st, err := f.Stat()
		if err != nil {
			return nil, err
		}
		p, err := parquet.Open(f, st.Size())
		if err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		var names []string
		for _, c := range p.Columns {
			names = append(names, c.Name)
		}
		return names, nil
	}
	var r io.Reader = f
	if strings.HasSuffix(name, ".gz") {
		zr, err := gzip.NewReader(f)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		defer zr.Close()
		r = zr
		name = strings.TrimSuffix(name, ".gz")
	}
	cr := csv.NewReader(bufio.NewReader(r))
	if strings.HasSuffix(name, ".tsv") || strings.HasSuffix(name, ".txt") {
		cr.Comma = '\t'
	}
	cr.FieldsPerRecord = -1
	head, err := cr.Read()
	if err != nil {
		return nil, fmt.Errorf("%s: the header: %w", path, err)
	}
	for i := range head {
		head[i] = strings.TrimSpace(strings.TrimPrefix(head[i], "\uFEFF"))
	}
	return head, nil
}

// DetectColumns names the instrument whose columns a header has ("xenium", "cosmx" or
// "merscope") and answers its columns; ok is false for a header of none of them.
func DetectColumns(header []string) (c Columns, instrument string, ok bool) {
	has := map[string]bool{}
	for _, h := range header {
		has[h] = true
	}
	switch {
	case has["x_location"] && has["y_location"] && has["feature_name"] && has["cell_id"]:
		return XeniumColumns(), "xenium", true
	case has["x_global_px"] && has["y_global_px"] && has["target"] && has["cell_ID"]:
		return CosMxColumns(), "cosmx", true
	case has["global_x"] && has["global_y"] && has["gene"] && has["cell_id"]:
		return MERSCOPEColumns(), "merscope", true
	}
	return Columns{}, "", false
}
