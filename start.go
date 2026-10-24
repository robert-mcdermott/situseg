package situseg

// A segmentation's inputs from a table of transcripts, as imaging-based spatial transcriptomics
// instruments write them: each molecule's place, its gene (or negative control), its quality, and
// the cell the instrument's own segmentation gave it, with whether it lies over that cell's
// nucleus. The cells to start from are the instrument's, drawn from their transcripts alone:
//   - the start raster: every pixel holding transcripts the instrument gave a cell is that cell's
//     (the cell most of them were given), and every other pixel the cell of the nearest such pixel
//     within RadiusUm (an exact Euclidean expansion, expand.go): each cell's region about its
//     transcripts, split between neighbors along the line midway between them;
//   - each cell's seed: the centroid of its transcripts over its nucleus, else of all of them;
//   - each modeled transcript's prior cell: its cell when it lies over the cell's nucleus (or,
//     when the table does not say, its cell whatever it lies over);
//   - the negative controls' places, which measure the false-decode rate.

import (
	"context"
	"errors"
	"fmt"
	"math"
	"sort"
)

// Table is a table of transcripts: one row per decoded molecule.
type Table struct {
	// X, Y and Z are each row's place (µm); Z is nil for a table without depths.
	X, Y, Z []float32
	// Feature is each row's feature (a gene or a negative control), an index into Features;
	// Control marks the features that are negative controls (probes or codewords that match no
	// gene: false decodes, not a gene's molecules).
	Feature  []uint32
	Features []string
	Control  []bool
	// QV is each row's quality score (Phred-scaled); nil for a table without them.
	QV []float32
	// Cell is the cell the instrument gave each row, an index into CellIDs, or −1 for none.
	Cell    []int32
	CellIDs []string
	// Nuclear marks the rows over their cell's nucleus; nil for a table that does not say.
	Nuclear []bool
	// ID is each row's transcript id; nil for a table without them.
	ID []uint64
}

// Rows is the table's number of rows.
func (t *Table) Rows() int { return len(t.X) }

// Check reports a table whose columns disagree.
func (t *Table) Check() error {
	n := len(t.X)
	switch {
	case len(t.Y) != n || len(t.Feature) != n || len(t.Cell) != n:
		return errors.New("situseg: the table's columns differ in length")
	case t.Z != nil && len(t.Z) != n, t.QV != nil && len(t.QV) != n, t.Nuclear != nil && len(t.Nuclear) != n, t.ID != nil && len(t.ID) != n:
		return errors.New("situseg: the table's optional columns differ in length from its rows")
	case len(t.Control) != len(t.Features):
		return errors.New("situseg: the table's control flags differ in number from its features")
	case n == 0:
		return errors.New("situseg: the table has no rows")
	}
	for i, f := range t.Feature {
		if int(f) >= len(t.Features) {
			return fmt.Errorf("situseg: row %d is of feature %d, of %d", i, f, len(t.Features))
		}
		if c := t.Cell[i]; c < -1 || int(c) >= len(t.CellIDs) {
			return fmt.Errorf("situseg: row %d is of cell %d, of %d", i, c, len(t.CellIDs))
		}
	}
	return nil
}

// StartOptions say how a table's inputs are made.
type StartOptions struct {
	// MinQV is the quality the modeled transcripts reach (when the table has quality scores).
	MinQV float64
	// PixelUm is the start raster's pixel side (µm), also the segmentation's (0: 0.5).
	PixelUm float64
	// RadiusUm is how far about its transcripts a cell's start region reaches (µm; 0: 2.5).
	RadiusUm float64
	// Workers share the expansion (0: 1).
	Workers int
}

// DefaultStartOptions are the defaults: transcripts of quality 20 or more, a raster of 0.5 µm
// pixels, and each cell's start region within 2.5 µm of its transcripts.
func DefaultStartOptions() StartOptions { return StartOptions{MinQV: 20, PixelUm: 0.5, RadiusUm: 2.5} }

// Index relates a segmentation's inputs to its table: the genes' names (the modeled transcripts'
// Gene values index them), each modeled transcript's row, and each cell's id.
type Index struct {
	Genes   []string
	Rows    []int
	CellIDs []string
	// Outside counts the gene transcripts beyond the start raster (none, for a table's own
	// inputs); Controls the rows of negative controls; LowQuality the gene transcripts below the
	// quality; WithPrior the modeled transcripts with a prior cell.
	Outside    int
	Controls   int
	LowQuality int
	WithPrior  int
	// Total is the table's rows; RowGroup and MaskBlock what a reader decodes at once (rows, and
	// a mask's pixels; 0 for a table in memory).
	Total     int
	RowGroup  int
	MaskBlock int
}

// InputsFromTable makes a segmentation's inputs from a table, its cells started from those the
// instrument gave the transcripts.
func InputsFromTable(ctx context.Context, t *Table, o StartOptions) (*Inputs, *Index, error) {
	if err := t.Check(); err != nil {
		return nil, nil, fail(ClassInput, remedyInput, "%s", plain(err))
	}
	if o.PixelUm == 0 {
		o.PixelUm = 0.5
	}
	if o.RadiusUm == 0 {
		o.RadiusUm = 2.5
	}
	if !(o.PixelUm > 0) || math.IsInf(o.PixelUm, 0) || !(o.RadiusUm >= 0) || math.IsInf(o.RadiusUm, 0) {
		return nil, nil, errors.New("situseg: the pixel size is positive and the start radius at least 0")
	}
	n := len(t.CellIDs)
	if n == 0 {
		return nil, nil, fail(ClassPrior, remedyPrior, "no transcript was given a cell: the table names no cells to start from")
	}
	// The extent: every finite place, with a margin of the radius and a pixel.
	minX, minY, maxX, maxY := math.Inf(1), math.Inf(1), math.Inf(-1), math.Inf(-1)
	finite := func(i int) bool {
		x, y := float64(t.X[i]), float64(t.Y[i])
		return !math.IsNaN(x) && !math.IsInf(x, 0) && !math.IsNaN(y) && !math.IsInf(y, 0)
	}
	for i := range t.X {
		if finite(i) {
			x, y := float64(t.X[i]), float64(t.Y[i])
			minX, minY, maxX, maxY = math.Min(minX, x), math.Min(minY, y), math.Max(maxX, x), math.Max(maxY, y)
		}
	}
	if math.IsInf(minX, 0) {
		return nil, nil, fail(ClassInput, remedyInput, "no transcript has a place")
	}
	ps := o.PixelUm
	pad := o.RadiusUm + ps
	x0, y0 := mul(math.Floor((minX-pad)/ps), ps), mul(math.Floor((minY-pad)/ps), ps)
	w, h := int(math.Ceil((maxX+pad-x0)/ps)), int(math.Ceil((maxY+pad-y0)/ps))
	if float64(w)*float64(h) > 1<<40 {
		return nil, nil, fail(ClassInput, remedyInput, "the transcripts span %.0f × %.0f µm, too large for %g µm pixels", maxX-minX, maxY-minY, ps)
	}
	pixel := func(i int) (int, int) {
		return int(math.Floor((float64(t.X[i]) - x0) / ps)), int(math.Floor((float64(t.Y[i]) - y0) / ps))
	}
	// The rows the instrument gave a cell, by pixel row, for the seed image.
	rowStart := make([]int32, h+1)
	for i, c := range t.Cell {
		if c >= 0 && finite(i) {
			_, py := pixel(i)
			rowStart[py+1]++
		}
	}
	for y := 1; y <= h; y++ {
		rowStart[y] += rowStart[y-1]
	}
	byRow := make([]int32, rowStart[h])
	next := append([]int32(nil), rowStart[:h]...)
	for i, c := range t.Cell {
		if c >= 0 && finite(i) {
			_, py := pixel(i)
			byRow[next[py]] = int32(i)
			next[py]++
		}
	}
	src := &seedSource{t: t, w: w, h: h, rowStart: rowStart, byRow: byRow, pixel: pixel}
	start, _, err := Expand(ctx, src, ExpandOptions{Distance: o.RadiusUm / ps, Workers: max(1, o.Workers)})
	if err != nil {
		return nil, nil, err
	}
	start.PixelSize = ps
	stats := start.Stats(n)

	in := &Inputs{X0: x0, Y0: y0, WidthUm: mul(float64(w), ps), HeightUm: mul(float64(h), ps), Cells: n, Start: start,
		SeedX: make([]float64, n), SeedY: make([]float64, n), Box: make([]Region, n),
		Origin: fmt.Sprintf("the cells the table's transcripts were given: each the pixels within %g µm of its transcripts (%g µm pixels), seeded at the centroid of its transcripts over its nucleus (else of all of them)", o.RadiusUm, ps)}
	for c := range in.Box {
		if s := stats[c]; s.Pixels > 0 {
			in.Box[c] = Region{x0 + mul(float64(s.MinX), ps), y0 + mul(float64(s.MinY), ps), x0 + mul(float64(s.MaxX+1), ps), y0 + mul(float64(s.MaxY+1), ps)}
		}
	}
	// Seeds: the centroid of each cell's nuclear transcripts, else of all of them.
	type sums struct {
		n      int
		sx, sy float64
	}
	all, nuclear := make([]sums, n), make([]sums, n)
	for i, c := range t.Cell {
		if c < 0 || !finite(i) {
			continue
		}
		x, y := float64(t.X[i]), float64(t.Y[i])
		s := &all[c]
		s.n, s.sx, s.sy = s.n+1, s.sx+x, s.sy+y
		if t.Nuclear != nil && t.Nuclear[i] {
			s := &nuclear[c]
			s.n, s.sx, s.sy = s.n+1, s.sx+x, s.sy+y
		}
	}
	if t.Nuclear != nil {
		in.Nucleus = make([]bool, n)
	}
	for c := 0; c < n; c++ {
		s := all[c]
		if nuclear[c].n > 0 {
			s = nuclear[c]
			in.Nucleus[c] = true
		}
		if s.n > 0 {
			in.SeedX[c], in.SeedY[c] = s.sx/float64(s.n), s.sy/float64(s.n)
		}
	}
	// The modeled transcripts: genes of the quality or more, in the table's order.
	idx := &Index{CellIDs: t.CellIDs, Total: t.Rows()}
	in.MinQV = o.MinQV
	geneOf := make([]int32, len(t.Features))
	for f := range geneOf {
		geneOf[f] = -1
	}
	controlKinds := map[uint32]bool{}
	zMin, zMax := math.Inf(1), math.Inf(-1)
	for i, f := range t.Feature {
		good := t.QV == nil || float64(t.QV[i]) >= o.MinQV
		if t.Control[f] {
			idx.Controls++
			controlKinds[f] = true
			if good && finite(i) {
				in.ControlX, in.ControlY = append(in.ControlX, t.X[i]), append(in.ControlY, t.Y[i])
			}
			continue
		}
		if !good {
			idx.LowQuality++
			continue
		}
		if !finite(i) {
			idx.Outside++
			continue
		}
		g := geneOf[f]
		if g < 0 {
			g = int32(len(idx.Genes))
			geneOf[f] = g
			idx.Genes = append(idx.Genes, t.Features[f])
		}
		prior := uint32(0)
		if c := t.Cell[i]; c >= 0 && (t.Nuclear == nil || t.Nuclear[i]) {
			prior = uint32(c) + 1
			idx.WithPrior++
		}
		z := float32(0)
		if t.Z != nil {
			z = t.Z[i]
		}
		zMin, zMax = math.Min(zMin, float64(z)), math.Max(zMax, float64(z))
		tr := &in.Transcripts
		tr.X, tr.Y, tr.Z = append(tr.X, t.X[i]), append(tr.Y, t.Y[i]), append(tr.Z, z)
		tr.Gene, tr.Prior = append(tr.Gene, uint32(g)), append(tr.Prior, prior)
		idx.Rows = append(idx.Rows, i)
	}
	if len(idx.Rows) == 0 {
		return nil, nil, fail(ClassInput, "Lower the quality threshold (--min-qv), or check which features are taken for negative controls.", "the table has no gene transcripts of quality %g or more", o.MinQV)
	}
	in.Genes = len(idx.Genes)
	in.ControlKinds = len(controlKinds)
	in.Z0, in.DepthUm = zMin, zMax-zMin
	return in, idx, nil
}

// seedSource is the seed image of a table's cells: each pixel holding transcripts the instrument
// gave a cell is the cell most of them were given (ties to the cell listed first), every other
// pixel none.
type seedSource struct {
	t        *Table
	w, h     int
	rowStart []int32
	byRow    []int32
	pixel    func(int) (int, int)
}

func (s *seedSource) Size() (int, int) { return s.w, s.h }

func (s *seedSource) Rows(y0, y1 int, dst []uint32) error {
	clear(dst[:(y1-y0)*s.w])
	type vote struct {
		x    int
		cell int32
	}
	var votes []vote
	for y := y0; y < y1; y++ {
		votes = votes[:0]
		for _, i := range s.byRow[s.rowStart[y]:s.rowStart[y+1]] {
			x, _ := s.pixel(int(i))
			votes = append(votes, vote{x, s.t.Cell[i]})
		}
		sort.Slice(votes, func(a, b int) bool {
			if votes[a].x != votes[b].x {
				return votes[a].x < votes[b].x
			}
			return votes[a].cell < votes[b].cell
		})
		row := dst[(y-y0)*s.w:][:s.w]
		for a := 0; a < len(votes); {
			b := a
			best, bestN := int32(-1), 0
			for b < len(votes) && votes[b].x == votes[a].x {
				c := b
				for c < len(votes) && votes[c].x == votes[a].x && votes[c].cell == votes[b].cell {
					c++
				}
				if c-b > bestN {
					best, bestN = votes[b].cell, c-b
				}
				b = c
			}
			if x := votes[a].x; x >= 0 && x < s.w {
				row[x] = uint32(best) + 1
			}
			a = b
		}
	}
	return nil
}
