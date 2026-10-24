package situseg

// The presets on a simulated tissue (SITUSEG_PRESETS=small or standard): each scored against the
// truth, with its time and peak heap; the numbers the presets' Measured texts quote. Each puts
// more transcripts in their true cells than the instrument; Conservative draws more compact
// outlines than Balanced, Flexible less compact ones; Low memory needs less memory; the ROI
// pilot's two seeds mostly agree.

import (
	"context"
	"os"
	"runtime"
	"sync/atomic"
	"testing"
	"time"

	"github.com/robert-mcdermott/situseg/internal/simulate"
)

// peakHeap samples the heap while fn runs and answers its largest.
func peakHeap(fn func()) uint64 {
	var peak atomic.Uint64
	done := make(chan struct{})
	go func() {
		var m runtime.MemStats
		for {
			runtime.ReadMemStats(&m)
			if m.HeapAlloc > peak.Load() {
				peak.Store(m.HeapAlloc)
			}
			select {
			case <-done:
				return
			case <-time.After(250 * time.Millisecond):
			}
		}
	}()
	fn()
	close(done)
	return peak.Load()
}

// compactness is the cells' mean perimeter, in pixel edges of their outlines' rings, over the least
// for their area (2⌈2√n⌉ for n pixels): 1 for the most compact.
func compactness(seg *Segmentation) float64 {
	var sum float64
	n := 0
	ps := seg.Labels.PixelSize
	for c := range seg.Outlines {
		o := &seg.Outlines[c]
		if len(o.X) < 3 || o.Pixels < 4 {
			continue
		}
		var perim float64
		for i := range o.X {
			j := (i + 1) % len(o.X)
			perim += (absf(float64(o.X[j]-o.X[i])) + absf(float64(o.Y[j]-o.Y[i]))) / ps
		}
		sum += perim / float64(minPerimeter(o.Pixels))
		n++
	}
	return sum / float64(max(n, 1))
}

func absf(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}

func TestPresetsOnSimulatedTissue(t *testing.T) {
	tissue := os.Getenv("SITUSEG_PRESETS")
	if tissue == "" {
		t.Skip("set SITUSEG_PRESETS=small or standard to measure the presets")
	}
	tis := simulated(t, tissue, 0, 0)
	tab := tableOf(tis)
	in, idx, err := InputsFromTable(context.Background(), tab, DefaultStartOptions())
	if err != nil {
		t.Fatal(err)
	}
	start := scoreOf(tis, tis.Cell, len(tis.CellIDs))
	t.Logf("instrument: %.1f%% in their true cell, median Jaccard %.3f", 100*start.InTrueCell, start.MedianJaccard)
	w, h := in.WidthUm, in.HeightUm
	pilot := Region{in.X0 + 0.3*w, in.Y0 + 0.3*h, in.X0 + 0.7*w, in.Y0 + 0.7*h}
	only := os.Getenv("SITUSEG_PRESETS_ONLY")
	for _, p := range Presets() {
		if p.ID == PresetExpert || only != "" && p.ID != only {
			continue
		}
		spec, _ := SpecOf(p.ID)
		if p.Area == AreaRegions {
			spec.Area.Regions = []Region{pilot}
		}
		var seg *Segmentation
		began := time.Now()
		heap := peakHeap(func() {
			seg, err = Segment(context.Background(), tab, in, idx, spec, Job{})
		})
		if err != nil {
			t.Fatal(err)
		}
		s := scoreOf(tis, seg.Assignment, len(seg.Cells))
		if p.Area == AreaRegions {
			// The pilot's region alone: its transcripts, by where they were observed.
			s = simulate.Compare(tis.TrueCell, seg.Assignment, len(tis.TrueIDs), len(seg.Cells), func(i int) bool {
				return !tis.Control[tis.Feature[i]] && tis.QV[i] >= 20 && pilot.Contains(float64(tis.X[i]), float64(tis.Y[i]))
			})
			whole, _ := segmentOf(t, tab, DefaultSpec(), 0)
			w := simulate.Compare(tis.TrueCell, whole.Assignment, len(tis.TrueIDs), len(whole.Cells), func(i int) bool {
				return !tis.Control[tis.Feature[i]] && tis.QV[i] >= 20 && pilot.Contains(float64(tis.X[i]), float64(tis.Y[i]))
			})
			t.Logf("%s: Balanced on the whole tissue, in the pilot's region: %.1f%% in their true cell, median Jaccard %.3f", p.ID, 100*w.InTrueCell, w.MedianJaccard)
		}
		d := seg.Diagnostics.Summary
		agree := 0.0
		if d.SeedAgreement != nil {
			agree = *d.SeedAgreement
		}
		t.Logf("%s: %.1f%% in their true cell, median Jaccard %.3f, outlines %.2f times the least; %.1f s, peak heap %d MB; Geweke z %.1f; seeds agree %.1f%%",
			p.ID, 100*s.InTrueCell, s.MedianJaccard, compactness(seg), time.Since(began).Seconds(), heap>>20, d.GewekeZ, 100*agree)
		if p.Area != AreaRegions && s.InTrueCell <= start.InTrueCell {
			t.Errorf("%s: %.1f%% in their true cell, no more than the instrument's %.1f%%", p.ID, 100*s.InTrueCell, 100*start.InTrueCell)
		}
	}
}
