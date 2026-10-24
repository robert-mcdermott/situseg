package situseg

// A segmentation of a table, stitched from its runs' results. Its cells are those it started from,
// one for one, with their ids:
//   - the label raster at the start raster's resolution: in each refined region, a refined cell's
//     pixels are those its 2D footprint holds (each lattice column's dominant label, holes
//     filled), and every other pixel keeps the start's, as do the pixels beyond the regions;
//   - the outlines traced from the raster, and with more than one layer each layer's;
//   - each row's cell: for a modeled transcript in a refined cell, the cell it lay in most often
//     over the summarized sweeps, with that share as its probability; for every other (a cell held
//     where it was, controls, low quality, beyond the regions), the cell holding the pixel it was
//     observed in;
//   - the counts of gene transcripts (of quality 20 or more when the table has quality scores),
//     by those assignments;
//   - the diagnostics: each region's traces and convergence, the cells flagged for what limits
//     them (few transcripts, no nucleus, uncertain boundaries), and those left out.

import (
	"math"
	"sort"
)

// The method as a segmentation records it.
const (
	Method        = "situseg transcript-informed segmentation"
	MethodVersion = "1"
	Citation      = "Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR, Newell EW. Cell simulation as cell segmentation. Nat Methods 2025;22(6):1331–1342. doi:10.1038/s41592-025-02697-0"
)

// countsMinQV is the quality of the transcripts the counts are of, for a table with quality
// scores (Xenium's Q20, a 1 in 100 chance the gene was misread).
const countsMinQV = 20

// Cell is one cell of a segmentation.
type Cell struct {
	// ID is the cell's id in the table it started from.
	ID string
	// X, Y are its centroid (µm) and AreaUm2 its area, from its pixels (0 for a cell with none).
	X, Y    float64
	AreaUm2 float64
	// Transcripts counts its gene transcripts (of the counts' quality), Genes the genes among
	// them, and Assigned every row assigned to it.
	Transcripts, Genes, Assigned int
	// Refined says whether a run refined it (else it was held, or lay beyond the regions);
	// Certainty is its boundary certainty, the mean share of the summarized sweeps its pixels were
	// its (NaN when not refined).
	Refined   bool
	Certainty float64
}

// Counts is a cell-by-gene matrix in compressed sparse rows: cell c's genes are
// Gene[Offsets[c]:Offsets[c+1]], with their counts.
type Counts struct {
	Genes   []string
	Offsets []uint32
	Gene    []uint32
	Count   []uint32
}

// Layer is one z layer's outlines: those of the cells refined in it (µm), and the layer's depths.
type Layer struct {
	Z0, Z1   float32
	Outlines []Outline
}

// Segmentation is a stitched segmentation of a table.
type Segmentation struct {
	Cells []Cell
	// Labels is the label raster: label k is Cells[k−1], 0 none; pixel (0, 0)'s corner is at
	// (X0, Y0) µm.
	Labels *Raster
	X0, Y0 float64
	// Assignment is each table row's cell (an index into Cells, −1 none); Probability its
	// posterior probability (NaN for a row assigned by place).
	Assignment  []int32
	Probability []float32
	// Counts are the cells' gene transcripts of CountsMinQV or more (all, for a table without
	// quality scores).
	Counts      Counts
	CountsMinQV float64
	// Outlines are the cells' outlines (µm); Layers each layer's, when the finest stage has more
	// than one.
	Outlines []Outline
	Layers   []Layer
	// HasNucleus says whether each cell has a nucleus (nil when the table does not say).
	HasNucleus  []bool
	Diagnostics *Diagnostics
	Recorded    Recorded
}

// Stitch makes a segmentation of a table from its runs' results (RunArea's).
func Stitch(t *Table, in *Inputs, idx *Index, results []*Result, opt Options, rec Recorded) (*Segmentation, error) {
	n := in.Cells
	img := in.Start
	ps := img.PixelSize
	// The cells each run refines (free, or left out), which decides each (−1: none), and each
	// run's footprint.
	refinedBy := make([]int, n)
	for c := range refinedBy {
		refinedBy[c] = -1
	}
	footprints := make([][]uint32, len(results))
	freeIn := make([]map[int]bool, len(results))
	for r, res := range results {
		footprints[r] = Footprint(res)
		free := map[int]bool{}
		for k, c := range res.Kept {
			if !res.Frozen[k] {
				free[c] = true
			}
		}
		for _, c := range res.LeftOut {
			free[c] = true
		}
		for c := range free {
			if res.decidesCell(in, c) {
				refinedBy[c] = r
			}
		}
		freeIn[r] = free
	}
	// The raster: the start's, each region's lattice (each patch's square) redrawn.
	raster := &Raster{Width: img.Width, Height: img.Height, PixelSize: ps, Rows: make([][]LabelRun, img.Height)}
	row := make([]uint32, img.Width)
	for y := 0; y < img.Height; y++ {
		img.Row(y, row)
		cyUm := in.Y0 + mul(float64(y), ps) + ps/2
		for r, res := range results {
			l := res.State.Lat
			cy := int(math.Floor((cyUm - l.Y0) / l.VoxelUm))
			if cy < 0 || cy >= l.H {
				continue
			}
			fp := footprints[r]
			// The pixels whose centers lie in the lattice.
			x0 := max(0, int(math.Floor((l.X0-in.X0)/ps)))
			x1 := min(len(row), int(math.Ceil((l.X0+mul(float64(l.W), l.VoxelUm)-in.X0)/ps))+1)
			for x := x0; x < x1; x++ {
				cxUm := in.X0 + mul(float64(x), ps) + ps/2
				cx := int(math.Floor((cxUm - l.X0) / l.VoxelUm))
				if cx < 0 || cx >= l.W || res.patched() && !res.Own.Contains(cxUm, cyUm) {
					continue
				}
				if f := fp[cy*l.W+cx]; f > 0 && !res.Frozen[f-1] {
					row[x] = uint32(res.Kept[f-1] + 1)
				} else if c := row[x]; c > 0 && freeIn[r][int(c)-1] {
					row[x] = 0 // a refined cell's pixel its footprint gave up
				}
			}
		}
		raster.Rows[y] = EncodeRow(row)
	}
	stats := raster.Stats(n)
	at := func(i int) uint32 { // the cell holding row i's pixel
		x, y := float64(t.X[i]), float64(t.Y[i])
		if math.IsNaN(x) || math.IsNaN(y) || math.IsInf(x, 0) || math.IsInf(y, 0) {
			return 0
		}
		return raster.AtMicrons(x-in.X0, y-in.Y0)
	}

	// Where each modeled transcript was refined: its result and index there.
	type place struct{ result, index int32 }
	modelledAt := make([]place, len(in.Transcripts.Gene))
	for i := range modelledAt {
		modelledAt[i] = place{-1, -1}
	}
	for r, res := range results {
		for i, j := range res.Index {
			if res.decidesTranscript(in, i) {
				modelledAt[j] = place{int32(r), int32(i)}
			}
		}
	}
	rowOf := make([]int32, idx.Total)
	for r := range rowOf {
		rowOf[r] = -1
	}
	for i, r := range idx.Rows {
		rowOf[r] = int32(i)
	}
	// Assignments, every row of the table; counts of gene transcripts by gene.
	seg := &Segmentation{Labels: raster, X0: in.X0, Y0: in.Y0, Assignment: make([]int32, idx.Total), Probability: make([]float32, idx.Total)}
	if t.QV != nil {
		seg.CountsMinQV = countsMinQV
	}
	nan := float32(math.NaN())
	genes := append([]string(nil), idx.Genes...)
	geneIndex := make(map[string]uint32, len(genes))
	for g, name := range genes {
		geneIndex[name] = uint32(g)
	}
	featureGene := make([]int32, len(t.Features)) // each feature's gene in the counts (−1: not yet)
	for f := range featureGene {
		featureGene[f] = -1
		if g, ok := geneIndex[t.Features[f]]; ok && !t.Control[f] {
			featureGene[f] = int32(g)
		}
	}
	// Each cell's genes and counts, (gene, count) pairs: a few dozen genes a cell, so a cell's list
	// is searched rather than hashed, and the memory is the pairs' own.
	perCell := make([][]uint32, n)
	for r := 0; r < idx.Total; r++ {
		c := uint32(0)
		p := nan
		assigned := false
		if j := rowOf[r]; j >= 0 {
			if pl := modelledAt[j]; pl.result >= 0 {
				res := results[pl.result]
				f := res.Assigned[pl.index]
				switch {
				case f > 0 && !res.Frozen[f-1]:
					c, p, assigned = uint32(res.Kept[f-1]+1), res.Probability[pl.index], true
				case f == 0:
					// The background, most often: so it stays, unless it lies on a held cell's
					// pixel (assigned by place, below).
					if pix := at(r); pix == 0 || freeIn[pl.result][int(pix)-1] {
						c, p, assigned = 0, res.Probability[pl.index], true
					}
				}
			}
		}
		if !assigned {
			c = at(r)
		}
		seg.Assignment[r], seg.Probability[r] = int32(c)-1, p
		f := t.Feature[r]
		if c == 0 || t.Control[f] || t.QV != nil && float64(t.QV[r]) < countsMinQV {
			continue
		}
		g := featureGene[f]
		if g < 0 {
			g = int32(len(genes))
			featureGene[f] = g
			genes = append(genes, t.Features[f])
		}
		list := perCell[c-1]
		found := false
		for k := 0; k < len(list); k += 2 {
			if list[k] == uint32(g) {
				list[k+1]++
				found = true
				break
			}
		}
		if !found {
			perCell[c-1] = append(list, uint32(g), 1)
		}
	}

	// Cells and counts.
	cells := make([]Cell, n)
	px2 := mul(ps, ps)
	for c := range cells {
		cell := &cells[c]
		cell.ID = idx.CellIDs[c]
		cell.Certainty = math.NaN()
		if s := stats[c]; s.Pixels > 0 {
			cell.AreaUm2 = mul(float64(s.Pixels), px2)
			cell.X, cell.Y = mul(s.SumX/float64(s.Pixels), ps)+in.X0, mul(s.SumY/float64(s.Pixels), ps)+in.Y0
		}
	}
	for _, c := range seg.Assignment {
		if c >= 0 {
			cells[c].Assigned++
		}
	}
	total := 0
	for _, list := range perCell {
		total += len(list) / 2
	}
	counts := Counts{Genes: genes, Offsets: make([]uint32, n+1), Gene: make([]uint32, 0, total), Count: make([]uint32, 0, total)}
	order := []int{}
	for c, list := range perCell {
		order = order[:0]
		for k := 0; k < len(list); k += 2 {
			order = append(order, k)
		}
		sort.Slice(order, func(i, j int) bool { return list[order[i]] < list[order[j]] })
		for _, k := range order {
			counts.Gene = append(counts.Gene, list[k])
			counts.Count = append(counts.Count, list[k+1])
			cells[c].Transcripts += int(list[k+1])
			cells[c].Genes++
		}
		counts.Offsets[c+1] = uint32(len(counts.Gene))
		perCell[c] = nil
	}
	seg.Counts = counts

	// Outlines: the raster's, simplified to half a voxel; and each layer's in the refined regions.
	finest := opt.Stages[len(opt.Stages)-1]
	outlines := make([]Outline, n)
	tolerance := math.Max(outlineTolerance, outlineTolerance*finest.VoxelUm/ps)
	for c := range outlines {
		outlines[c] = raster.TraceOutline(uint32(c+1), stats[c], tolerance)
		shift(&outlines[c], in.X0, in.Y0)
	}
	seg.Outlines = outlines
	if finest.Layers > 1 {
		layers := make([]Layer, finest.Layers)
		for z := range layers {
			layers[z].Outlines = make([]Outline, n)
		}
		for _, res := range results {
			l := res.State.Lat
			for z := 0; z < l.L; z++ {
				layers[z].Z0, layers[z].Z1 = float32(l.Z0+mul(float64(z), l.LayerUm)), float32(l.Z0+mul(float64(z+1), l.LayerUm))
				lr := &Raster{Width: l.W, Height: l.H, PixelSize: l.VoxelUm, Rows: make([][]LabelRun, l.H)}
				plane := res.Labels[z*l.W*l.H : (z+1)*l.W*l.H]
				for y := 0; y < l.H; y++ {
					lr.Rows[y] = EncodeRow(plane[y*l.W : (y+1)*l.W])
				}
				ls := lr.Stats(len(res.Kept))
				for k, c := range res.Kept {
					if res.Frozen[k] || res.patched() && !res.decidesCell(in, c) {
						continue
					}
					o := lr.TraceOutline(uint32(k+1), ls[k], outlineTolerance)
					shift(&o, l.X0, l.Y0)
					layers[z].Outlines[c] = o
				}
			}
		}
		seg.Layers = layers
	}

	// Each refined cell's boundary certainty (NaN for the cells held, and those outside the regions).
	for _, res := range results {
		for k, x := range boundaryCertainty(res) {
			if c := res.Kept[k]; !res.Frozen[k] && res.decidesCell(in, c) {
				cells[c].Certainty = x
			}
		}
	}
	for c := range cells {
		cells[c].Refined = refinedBy[c] >= 0
	}
	seg.Cells = cells

	// Diagnostics, and what was recorded.
	diag := diagnostics(results, cells, in, idx, refinedBy)
	if len(rec.Second) == len(results) && len(results) > 0 {
		var same, all int
		for r, res := range results {
			sr, ar := agreement(res, rec.Second[r], in)
			x := float64(sr) / float64(max(ar, 1))
			diag.Regions[r].Summary.SeedAgreement = &x
			same += sr
			all += ar
		}
		x := float64(same) / float64(max(all, 1))
		diag.Summary.SeedAgreement = &x
	} else {
		rec.SecondSeed = 0
	}
	rec.Regions, rec.StartStages, rec.Patches = nil, nil, nil
	for _, res := range results {
		rec.Regions = append(rec.Regions, res.Region)
		rec.StartStages = append(rec.StartStages, res.StartStage+1)
		if res.patched() {
			rec.Patches = append(rec.Patches, res.Own)
		}
	}
	if len(results) > 0 {
		rec.Model = results[0].Model.H
	}
	rec.Sampler, rec.Stages, rec.BlockUm, rec.Seed, rec.Samples = opt.Sampler, opt.Stages, opt.BlockUm, opt.Seed, opt.Samples
	rec.MinQV, rec.CountsMinQuality = in.MinQV, seg.CountsMinQV
	if rec.Start == "" {
		rec.Start = in.Origin
	}
	if rec.Area.Mode == "" {
		rec.Area.Mode = AreaWhole
	}
	seg.Diagnostics, seg.Recorded, seg.HasNucleus = diag, rec, in.Nucleus
	return seg, nil
}
