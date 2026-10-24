package situseg

// On simulated tissues (internal/simulate), started from the simulated instrument's cells as its
// transcripts carry them: the segmentation assigns more transcripts to their true cells, and draws
// cells that match the true ones better, than the instrument does; a seed gives the same result on
// any number of threads and on every processor; a region refines its cells and keeps the rest;
// patches stitch into one segmentation.

import (
	"context"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"hash"
	"math"
	"runtime"
	"sync"
	"testing"

	"github.com/robert-mcdermott/situseg/internal/simulate"
)

var tissues sync.Map // preset (and size) → *simulate.Tissue, made once

// simulated makes (once) a simulated tissue of a preset, optionally cut to a size (µm).
func simulated(t testing.TB, preset string, w, h float64) *simulate.Tissue {
	t.Helper()
	key := [3]any{preset, w, h}
	if v, ok := tissues.Load(key); ok {
		return v.(*simulate.Tissue)
	}
	p := simulate.Presets[preset]
	if w > 0 {
		p.WidthUm, p.HeightUm = w, h
	}
	tis, err := simulate.Simulate(context.Background(), p)
	if err != nil {
		t.Fatal(err)
	}
	tissues.Store(key, tis)
	return tis
}

// tableOf is a simulated tissue's table, as an instrument writes it.
func tableOf(tis *simulate.Tissue) *Table {
	return &Table{X: tis.X, Y: tis.Y, Z: tis.Z, Feature: tis.Feature, Features: tis.Features, Control: tis.Control,
		QV: tis.QV, Cell: tis.Cell, CellIDs: tis.CellIDs, Nuclear: tis.Nuclear, ID: tis.ID}
}

// scoreOf scores an assignment of a tissue's rows against the truth, on its gene transcripts of
// quality 20 or more.
func scoreOf(tis *simulate.Tissue, assignment []int32, cells int) simulate.Score {
	return simulate.Compare(tis.TrueCell, assignment, len(tis.TrueIDs), cells, func(i int) bool {
		return !tis.Control[tis.Feature[i]] && tis.QV[i] >= 20
	})
}

// segmentOf segments a table with a spec on a number of threads.
func segmentOf(t testing.TB, tab *Table, spec Spec, workers int) (*Segmentation, *Inputs) {
	t.Helper()
	in, idx, err := InputsFromTable(context.Background(), tab, DefaultStartOptions())
	if err != nil {
		t.Fatal(err)
	}
	seg, err := Segment(context.Background(), tab, in, idx, spec, Job{Workers: workers})
	if err != nil {
		t.Fatal(err)
	}
	return seg, in
}

// shortSpec is Balanced with fewer sweeps, for the tests.
func shortSpec(sweeps, samples int) Spec {
	s := DefaultSpec()
	s.Options.Stages = []Stage{{VoxelUm: 1, Layers: 1, Sweeps: sweeps}}
	s.Options.Samples = samples
	return s
}

// The segmentation beats its start: more transcripts in their true cell, cells closer to the true
// ones. On the small tissue (300 sweeps, the last 250 summarized).
func TestSegmentationBeatsTheInstrument(t *testing.T) {
	if testing.Short() {
		t.Skip("short mode")
	}
	tis := simulated(t, "small", 0, 0)
	tab := tableOf(tis)
	start := scoreOf(tis, tis.Cell, len(tis.CellIDs))
	seg, _ := segmentOf(t, tab, shortSpec(300, 250), 0)
	got := scoreOf(tis, seg.Assignment, len(seg.Cells))
	t.Logf("instrument: %.1f%% in their true cell, median Jaccard %.3f; segmentation %.1f%%, %.3f", 100*start.InTrueCell, start.MedianJaccard, 100*got.InTrueCell, got.MedianJaccard)
	if got.InTrueCell < start.InTrueCell+0.02 {
		t.Errorf("%.1f%% of transcripts in their true cell, the instrument's %.1f%%: not 2 points better", 100*got.InTrueCell, 100*start.InTrueCell)
	}
	if got.MedianJaccard < start.MedianJaccard+0.03 {
		t.Errorf("median Jaccard %.3f, the instrument's %.3f", got.MedianJaccard, start.MedianJaccard)
	}
	if got.MatchedCells < start.MatchedCells-start.MatchedCells/50 {
		t.Errorf("%d cells matched, the instrument %d", got.MatchedCells, start.MatchedCells)
	}
}

// shortRun is a run on a part of the small tissue, its options shortened.
func shortRun(t *testing.T, workers int) *Result {
	t.Helper()
	tab := tableOf(simulated(t, "small", 80, 60))
	in, _, err := InputsFromTable(context.Background(), tab, DefaultStartOptions())
	if err != nil {
		t.Fatal(err)
	}
	opt := DefaultOptions()
	opt.Stages = []Stage{{2, 1, 10}, {1, 2, 30}}
	opt.Samples, opt.Workers, opt.BlockUm = 10, workers, 16
	res, err := Run(context.Background(), in, Region{}, opt, nil)
	if err != nil {
		t.Fatal(err)
	}
	return res
}

// digestOf hashes everything a run makes.
func digestOf(res *Result) string {
	h := sha256.New()
	put := func(h hash.Hash, v ...float64) {
		for _, x := range v {
			binary.Write(h, binary.LittleEndian, math.Float64bits(x))
		}
	}
	for _, v := range [][]uint32{res.Labels, res.Assigned, res.State.Labels} {
		binary.Write(h, binary.LittleEndian, v)
	}
	binary.Write(h, binary.LittleEndian, res.Probability)
	put(h, res.State.posX...)
	put(h, res.State.posY...)
	m := res.Model
	put(h, m.Rate...)
	put(h, m.Shape...)
	put(h, m.RateParam...)
	put(h, m.Bg...)
	put(h, m.Ambient...)
	put(h, m.False)
	put(h, m.Pi...)
	put(h, m.VolMean...)
	for _, tr := range res.Trace {
		put(h, float64(tr.Moved), float64(tr.Repositioned), tr.LogLikelihood)
	}
	return hex.EncodeToString(h.Sum(nil))[:16]
}

// A seed gives the same run, byte for byte, on 1, 4 and 16 threads, and on every processor (the
// pinned digest, made on arm64, is checked on amd64 by CI).
func TestRunsAreTheSameOnAnyThreadsAndProcessor(t *testing.T) {
	one := shortRun(t, 1)
	want := digestOf(one)
	for _, workers := range []int{4, 16} {
		if got := digestOf(shortRun(t, workers)); got != want {
			t.Errorf("%d threads: digest %s, 1 thread %s", workers, got, want)
		}
	}
	const pinned = "3185a011f3c50265"
	if want != pinned {
		t.Errorf("the run's digest is %s, not %s (on %s)", want, pinned, runtime.GOARCH)
	}
}

// A region refines the cells whose seeds lie in it and keeps every other cell where it started.
func TestARegionKeepsTheRest(t *testing.T) {
	if testing.Short() {
		t.Skip("short mode")
	}
	tis := simulated(t, "small", 0, 0)
	tab := tableOf(tis)
	spec := shortSpec(120, 100)
	region := Region{X0: 50, Y0: 40, X1: 110, Y1: 80}
	spec.Area = Area{Mode: AreaRegions, Regions: []Region{region}}
	seg, in := segmentOf(t, tab, spec, 0)
	refined, held := 0, 0
	for c, cell := range seg.Cells {
		inside := region.Contains(in.SeedX[c], in.SeedY[c])
		if cell.Refined != inside {
			t.Fatalf("cell %s (seed %.1f, %.1f) refined %v", cell.ID, in.SeedX[c], in.SeedY[c], cell.Refined)
		}
		if inside {
			refined++
		} else {
			held++
		}
	}
	if refined == 0 || held == 0 {
		t.Fatalf("%d cells refined, %d kept", refined, held)
	}
	// Far from the region (beyond its margin), every pixel is the start's.
	ps := in.Start.PixelSize
	far := region.grown(marginUm + 2)
	row, start := make([]uint32, in.Start.Width), make([]uint32, in.Start.Width)
	for y := 0; y < in.Start.Height; y++ {
		seg.Labels.Row(y, row)
		in.Start.Row(y, start)
		for x := range row {
			cx, cy := in.X0+(float64(x)+0.5)*ps, in.Y0+(float64(y)+0.5)*ps
			if !far.Contains(cx, cy) && row[x] != start[x] {
				t.Fatalf("pixel (%d, %d) far from the region changed from %d to %d", x, y, start[x], row[x])
			}
		}
	}
	t.Logf("%d cells refined, %d kept", refined, held)
}

// Patches: a tissue cut into patches is stitched into one segmentation, close to the whole run's.
func TestPatchesAreStitched(t *testing.T) {
	if testing.Short() {
		t.Skip("short mode")
	}
	tis := simulated(t, "small", 0, 0)
	tab := tableOf(tis)
	whole, _ := segmentOf(t, tab, shortSpec(150, 120), 0)
	spec := shortSpec(150, 120)
	spec.Area.PatchTranscripts = 10_000
	patched, _ := segmentOf(t, tab, spec, 0)
	if patched.Diagnostics.Summary.Patches < 4 {
		t.Fatalf("%d patches", patched.Diagnostics.Summary.Patches)
	}
	if len(patched.Cells) != len(whole.Cells) {
		t.Fatalf("%d cells patched, %d whole", len(patched.Cells), len(whole.Cells))
	}
	sw, sp := scoreOf(tis, whole.Assignment, len(whole.Cells)), scoreOf(tis, patched.Assignment, len(patched.Cells))
	t.Logf("%d patches: %.1f%% in their true cell (whole %.1f%%), median Jaccard %.3f (%.3f)", patched.Diagnostics.Summary.Patches,
		100*sp.InTrueCell, 100*sw.InTrueCell, sp.MedianJaccard, sw.MedianJaccard)
	if sp.InTrueCell < sw.InTrueCell-0.015 || sp.MedianJaccard < sw.MedianJaccard-0.03 {
		t.Errorf("patched: %.1f%% and %.3f, whole %.1f%% and %.3f", 100*sp.InTrueCell, sp.MedianJaccard, 100*sw.InTrueCell, sw.MedianJaccard)
	}
}

// A seed check reports the two seeds' agreement, high on a tissue the model fits.
func TestSeedCheck(t *testing.T) {
	if testing.Short() {
		t.Skip("short mode")
	}
	tab := tableOf(simulated(t, "small", 80, 60))
	spec := shortSpec(120, 100)
	spec.SeedCheck = true
	seg, _ := segmentOf(t, tab, spec, 0)
	a := seg.Diagnostics.Summary.SeedAgreement
	if a == nil || *a < 0.8 || *a > 1 {
		t.Fatalf("seed agreement %v", a)
	}
	t.Logf("two seeds agree on %.1f%% of the transcripts", 100**a)
}
