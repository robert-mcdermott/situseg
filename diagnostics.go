package situseg

// A segmentation's diagnostics: each region's convergence (its trace, the split-half agreement of
// the summarized sweeps, Geweke's diagnostic), the background's shares, the cells flagged for what
// limits them (few transcripts, no nucleus, uncertain boundaries) and those left out; and the 2D
// footprint of a run's summarized labelling.

import (
	"math"
)

// Thresholds of the cells' limitations (this package's choices, docs/method.md).
const (
	// A cell with fewer gene transcripts than this has boundaries its prior draws more than its
	// transcripts.
	lowTranscripts = 10
	// A cell whose voxels' labels are, on average, its own in fewer of the summarized sweeps than
	// this has uncertain boundaries.
	uncertainShare = 0.75
	// outlineTolerance simplifies outlines to within half a voxel (in pixels of the lattice); at
	// the image's resolution, its pixels are finer than the voxels, and the outline is simplified to
	// half a voxel too, so that the voxels' steps are not drawn.
	outlineTolerance = 0.5
)

// Recorded are a segmentation's parameters as it records them: every option after defaults (the
// volume prior's mean resolved), what it covered and started from.
type Recorded struct {
	Model   Hyper    `json:"model"`
	Sampler Settings `json:"sampler"`
	Stages  []Stage  `json:"stages"`
	BlockUm float64  `json:"block_um"`
	Seed    uint64   `json:"seed"`
	Samples int      `json:"samples"`
	// Area is what was asked for; Regions the rectangles refined. The caller gives the area and
	// the preset; the rest is filled from the options, the inputs and the results.
	Area    Area     `json:"area"`
	Regions []Region `json:"regions"`
	// Patches are each patch's square, when a region was refined in patches (Regions then holds
	// each patch's region of refined cells).
	Patches []Region `json:"patches,omitempty"`
	// MinQV is the quality the modelled transcripts reach; the counts are of gene transcripts of
	// CountsMinQuality or more (every gene transcript, for a table without quality scores).
	MinQV            float64 `json:"min_qv"`
	CountsMinQuality float64 `json:"counts_min_quality"`
	// Start says what the cells started from; StartStage the stage each region started at (from 1).
	Start       string `json:"start"`
	StartStages []int  `json:"start_stages"`
	// Preset is the preset the options came from ("" for none), and the presets' version.
	Preset         string `json:"preset,omitempty"`
	PresetsVersion string `json:"presets_version,omitempty"`
	// SecondSeed is the seed of a second run checking the first (0: none), Second its results.
	SecondSeed uint64    `json:"second_seed,omitempty"`
	Second     []*Result `json:"-"`
}

// Limitations are the cells a segmentation flags, by id.
type Limitations struct {
	LowTranscripts []string `json:"lowTranscripts"`
	NoNucleus      []string `json:"noNucleus"`
	Uncertain      []string `json:"uncertainBoundaries"`
	LeftOut        []string `json:"leftOut"`
	Thresholds     struct {
		LowTranscripts int     `json:"lowTranscripts"`
		UncertainShare float64 `json:"uncertainShare"`
	} `json:"thresholds"`
}

// DiagnosticsSummary is a region's, or a segmentation's, convergence and flagged cells.
type DiagnosticsSummary struct {
	Sweeps  int `json:"sweeps"`
	Samples int `json:"samples"`
	// SplitHalf is the share of modelled transcripts whose cell is the same over both halves of
	// the samples; MeanProbability their assignments' mean probability; LowProbability the share
	// below 0.5.
	SplitHalf       float64 `json:"splitHalfAgreement"`
	MeanProbability float64 `json:"meanProbability"`
	LowProbability  float64 `json:"lowProbabilityShare"`
	// LogLikelihoodDrift is the log-likelihood's mean change from the first half of the samples to
	// the second, per modelled transcript (near 0 when it has settled).
	LogLikelihoodDrift float64 `json:"logLikelihoodDriftPerTranscript"`
	// GewekeZ is Geweke's convergence diagnostic of the log-likelihood over the summarized sweeps:
	// the difference between its mean over their first tenth and their last half, in standard
	// errors (batch means). Within ±2 it has settled. Over several regions, the one farthest from 0.
	// GewekeJudged says whether the sweeps summarized were enough to judge (100 or more; in every
	// region).
	GewekeZ      float64 `json:"gewekeZ"`
	GewekeJudged bool    `json:"gewekeJudged"`
	// MovedLast is the voxels moved per sweep over the last tenth of the sweeps.
	MovedLast float64 `json:"voxelsMovedPerSweep"`
	// The cells refined and held, the cells flagged, by limitation, and left out.
	Refined        int `json:"cellsRefined"`
	Held           int `json:"cellsHeld"`
	LowTranscripts int `json:"cellsLowTranscripts"`
	NoNucleus      int `json:"cellsNoNucleus"`
	Uncertain      int `json:"cellsUncertainBoundaries"`
	LeftOut        int `json:"cellsLeftOut"`
	// Regions are the regions refined; Patches the patches they were refined in, if any.
	Regions int `json:"regions"`
	Patches int `json:"patches,omitempty"`
	// Transcripts are the modelled transcripts (in the regions' lattices).
	Transcripts int `json:"transcripts"`
	// The background (docs/method.md): the shares of the modelled transcripts its rates account
	// for, false decodes (measured by the NegativeControls in the regions) and ambient RNA (less
	// surely: where cells cover the tissue, the ambient RNA of genes they express cannot be told
	// from their own); InBackground the share the summary leaves in no cell.
	FalseDecodeShare float64 `json:"falseDecodeShare"`
	AmbientShare     float64 `json:"ambientShare"`
	NegativeControls int     `json:"negativeControls"`
	InBackground     float64 `json:"inBackgroundShare"`
	// SeedAgreement is the share of the transcripts a run with a second seed puts in the same cell,
	// of those either run puts in a refined cell (absent without a second run).
	SeedAgreement *float64 `json:"seedAgreement,omitempty"`
}

// RegionDiagnostics are one region's, or one patch's (Patch its square).
type RegionDiagnostics struct {
	Region  Region             `json:"region"`
	Patch   *Region            `json:"patch,omitempty"`
	Summary DiagnosticsSummary `json:"summary"`
	Trace   []SweepTrace       `json:"trace"`
}

// Diagnostics is diagnostics.json.
type Diagnostics struct {
	Summary     DiagnosticsSummary  `json:"summary"`
	Limitations Limitations         `json:"limitations"`
	Regions     []RegionDiagnostics `json:"regions"`
}

// patched: a patch, which decides only its square; a region decides its whole lattice.
func (r *Result) patched() bool { return r.Own != r.Region }

// decidesCell: the run decides the cell it refines (a region every one, a patch those whose seeds
// lie in its square), so the stitched segmentation takes the cell's certainty and flags from it.
func (r *Result) decidesCell(in *Inputs, c int) bool {
	return !r.patched() || r.Own.Contains(in.SeedX[c], in.SeedY[c])
}

// decidesTranscript: the run's transcript i lies where the stitched segmentation takes its cell
// from the run (a region's every one, a patch's in its square).
func (r *Result) decidesTranscript(in *Inputs, i int) bool {
	if !r.patched() {
		return true
	}
	j := r.Index[i]
	return r.Own.Contains(float64(in.Transcripts.X[j]), float64(in.Transcripts.Y[j]))
}

// shift moves an outline by an origin (µm).
func shift(o *Outline, x0, y0 float64) {
	if x0 == 0 && y0 == 0 {
		return
	}
	for i := range o.X {
		o.X[i] += float32(x0)
		o.Y[i] += float32(y0)
	}
}

// regionSummary summarizes one region's convergence and flags its refined cells.
func regionSummary(res *Result, cells []Cell, in *Inputs, lim *Limitations, refinedBy []int) DiagnosticsSummary {
	var s DiagnosticsSummary
	s.Sweeps, s.SplitHalf, s.Samples, s.Regions = len(res.Trace), res.SplitHalf, res.Samples, 1
	if res.patched() {
		s.Patches = 1
	}
	samples := min(res.Samples, len(res.Trace))
	// The transcripts it decides (a region every one; a patch those in its square); its background's
	// shares are of all it modelled, as its rates were.
	var sum float64
	low, decided, background := 0, 0, 0
	for i, p := range res.Probability {
		if !res.decidesTranscript(in, i) {
			continue
		}
		decided++
		sum += float64(p)
		if p < 0.5 {
			low++
		}
		if res.Assigned[i] == 0 {
			background++
		}
	}
	s.Transcripts = decided
	n := float64(max(1, decided))
	s.MeanProbability, s.LowProbability, s.InBackground = sum/n, float64(low)/n, float64(background)/n
	all := float64(max(1, len(res.Probability)))
	s.FalseDecodeShare, s.AmbientShare, s.NegativeControls = res.FalseDecodes/all, res.Ambient/all, res.Controls
	if res.patched() {
		// The negative controls in its square.
		s.NegativeControls = 0
		for i := range in.ControlX {
			if res.Own.Contains(float64(in.ControlX[i]), float64(in.ControlY[i])) {
				s.NegativeControls++
			}
		}
	}
	if samples >= 2 {
		window := res.Trace[len(res.Trace)-samples:]
		var first, second float64
		h := samples / 2
		ll := make([]float64, samples)
		for i, t := range window {
			ll[i] = t.LogLikelihood
			if i < h {
				first += t.LogLikelihood / float64(h)
			} else {
				second += t.LogLikelihood / float64(samples-h)
			}
		}
		s.LogLikelihoodDrift = (second - first) / all
		s.GewekeZ, s.GewekeJudged = gewekeZ(ll)
	}
	tenth := max(1, len(res.Trace)/10)
	for _, t := range res.Trace[len(res.Trace)-tenth:] {
		s.MovedLast += float64(t.Moved) / float64(tenth)
	}
	// The refined cells' limitations.
	certainty := boundaryCertainty(res)
	for k, c := range res.Kept {
		if res.Frozen[k] {
			// Held, unless another run refines it.
			if refinedBy == nil || refinedBy[c] < 0 {
				s.Held++
			}
			continue
		}
		if !res.decidesCell(in, c) {
			continue // a patch's cell another patch decides
		}
		s.Refined++
		id := cells[c].ID
		if cells[c].Transcripts < lowTranscripts {
			lim.LowTranscripts = append(lim.LowTranscripts, id)
			s.LowTranscripts++
		}
		if in.Nucleus != nil && !in.Nucleus[c] {
			lim.NoNucleus = append(lim.NoNucleus, id)
			s.NoNucleus++
		}
		if !(certainty[k] >= uncertainShare) {
			lim.Uncertain = append(lim.Uncertain, id)
			s.Uncertain++
		}
	}
	for _, c := range res.LeftOut {
		if res.decidesCell(in, c) {
			s.LeftOut++
		}
	}
	return s
}

// boundaryCertainty is each of a region's cells' mean share of the summarized samples its voxels
// were its (its voxels by their most frequent label; NaN for a cell with none).
func boundaryCertainty(res *Result) []float64 {
	shares := make([]float64, len(res.Kept))
	voxels := make([]int, len(res.Kept))
	for v, c := range res.Labels {
		if c > 0 {
			shares[c-1] += float64(res.LabelShare[v])
			voxels[c-1]++
		}
	}
	for k := range shares {
		shares[k] /= float64(voxels[k]) // 0/0: NaN
	}
	return shares
}

// diagnostics summarizes a segmentation: each region's convergence, the whole (weighted by each
// region's modelled transcripts; Geweke's z the one farthest from 0), and the cells flagged.
// refinedBy is the run deciding each cell (−1 none; nil: each run's own).
func diagnostics(results []*Result, cells []Cell, in *Inputs, idx *Index, refinedBy []int) *Diagnostics {
	d := &Diagnostics{}
	lim := &d.Limitations
	lim.Thresholds.LowTranscripts, lim.Thresholds.UncertainShare = lowTranscripts, uncertainShare
	lim.LowTranscripts, lim.NoNucleus, lim.Uncertain, lim.LeftOut = []string{}, []string{}, []string{}, []string{}
	all := &d.Summary
	all.GewekeJudged = len(results) > 0
	var weight float64
	regions := map[Region]bool{}
	for _, res := range results {
		s := regionSummary(res, cells, in, lim, refinedBy)
		for _, c := range res.LeftOut {
			if res.decidesCell(in, c) {
				lim.LeftOut = append(lim.LeftOut, idx.CellIDs[c])
			}
		}
		rd := RegionDiagnostics{Region: res.Region, Summary: s, Trace: res.Trace}
		if res.patched() {
			own := res.Own
			rd.Patch = &own
			all.Patches++
		}
		d.Regions = append(d.Regions, rd)
		w := float64(s.Transcripts)
		weight += w
		all.Sweeps, all.Samples = s.Sweeps, s.Samples
		all.SplitHalf += mul(w, s.SplitHalf)
		all.MeanProbability += mul(w, s.MeanProbability)
		all.LowProbability += mul(w, s.LowProbability)
		all.LogLikelihoodDrift += mul(w, s.LogLikelihoodDrift)
		all.FalseDecodeShare += mul(w, s.FalseDecodeShare)
		all.AmbientShare += mul(w, s.AmbientShare)
		all.InBackground += mul(w, s.InBackground)
		all.NegativeControls += s.NegativeControls
		if math.Abs(s.GewekeZ) > math.Abs(all.GewekeZ) {
			all.GewekeZ = s.GewekeZ
		}
		all.GewekeJudged = all.GewekeJudged && s.GewekeJudged
		all.MovedLast += s.MovedLast
		all.Refined += s.Refined
		all.Held += s.Held
		all.LowTranscripts += s.LowTranscripts
		all.NoNucleus += s.NoNucleus
		all.Uncertain += s.Uncertain
		all.LeftOut += s.LeftOut
		if !res.patched() || !regions[res.Whole] {
			regions[res.Whole] = true
			all.Regions++
		}
		all.Transcripts += s.Transcripts
	}
	if weight > 0 {
		all.SplitHalf /= weight
		all.MeanProbability /= weight
		all.LowProbability /= weight
		all.LogLikelihoodDrift /= weight
		all.FalseDecodeShare /= weight
		all.AmbientShare /= weight
		all.InBackground /= weight
	}
	return d
}

// agreement counts the modelled transcripts two runs of a region put in the same cell (or both
// in the background), of those either puts in a refined cell.
func agreement(a, b *Result, in *Inputs) (same, all int) {
	cellOf := func(res *Result, i int) (int, bool) {
		k := res.Assigned[i]
		if k == 0 {
			return -1, false
		}
		return res.Kept[k-1], !res.Frozen[k-1]
	}
	for i := range a.Assigned {
		ca, fa := cellOf(a, i)
		cb, fb := cellOf(b, i)
		if !fa && !fb || !a.decidesTranscript(in, i) {
			continue
		}
		all++
		if ca == cb {
			same++
		}
	}
	return same, all
}

// gewekeZ is Geweke's convergence diagnostic of a trace (Geweke J. Evaluating the accuracy of
// sampling-based approaches to the calculation of posterior moments. In: Bayesian Statistics 4,
// Oxford University Press 1992:169–193): the means of its first tenth and its last half compared,
// each mean's standard error from batch means (batches of a fifth of the segment, at least 2), so
// that the samples' autocorrelation is allowed for. It is not judged (0, false) for a trace too
// short, under minGeweke samples; a constant trace has settled (0, true).
func gewekeZ(x []float64) (float64, bool) {
	n := len(x)
	a, b := x[:n/10], x[n/2:]
	meanSE := func(v []float64) (float64, float64, bool) {
		batches := 5
		size := len(v) / batches
		if size < 2 {
			return 0, 0, false
		}
		var mean float64
		means := make([]float64, batches)
		for k := range means {
			for _, y := range v[k*size : (k+1)*size] {
				means[k] += y / float64(size)
			}
			mean += means[k] / float64(batches)
		}
		var ss float64
		for _, m := range means {
			ss += mul(m-mean, m-mean)
		}
		return mean, math.Sqrt(ss / float64(batches-1) / float64(batches)), true
	}
	ma, sa, okA := meanSE(a)
	mb, sb, okB := meanSE(b)
	if !okA || !okB {
		return 0, false
	}
	if sa+sb == 0 {
		return 0, true
	}
	return (mb - ma) / math.Sqrt(mul(sa, sa)+mul(sb, sb)), true
}

// minGeweke is the fewest samples Geweke's diagnostic judges: a first tenth of five batches of
// two.
const minGeweke = 100

// Footprint is a run's 2D segmentation at its lattice, x then y: each column's dominant label
// (Lattice.Consensus) of the summarized labels, with every hole filled. A region of the
// background enclosed by one cell (not reaching the region's edge, every neighbour that cell's)
// is the cell's: a posterior summary leaves holes where a cell's voxels happened to hold no
// transcript, but cells have none (this package's choice, docs/method.md).
func Footprint(res *Result) []uint32 {
	l := res.State.Lat
	cols := l.Consensus(res.Labels)
	seen := make([]bool, len(cols))
	var region, stack []int
	for start, c := range cols {
		if c != 0 || seen[start] {
			continue
		}
		region, stack = region[:0], append(stack[:0], start)
		seen[start] = true
		edge := false
		var around uint32
		mixed := false
		for len(stack) > 0 {
			p := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			region = append(region, p)
			x, y := p%l.W, p/l.W
			if x == 0 || y == 0 || x == l.W-1 || y == l.H-1 {
				edge = true
			}
			for _, q := range [4][2]int{{x - 1, y}, {x + 1, y}, {x, y - 1}, {x, y + 1}} {
				if q[0] < 0 || q[1] < 0 || q[0] >= l.W || q[1] >= l.H {
					continue
				}
				n := q[1]*l.W + q[0]
				switch {
				case cols[n] == 0:
					if !seen[n] {
						seen[n] = true
						stack = append(stack, n)
					}
				case around == 0:
					around = cols[n]
				case cols[n] != around:
					mixed = true
				}
			}
		}
		if !edge && !mixed && around != 0 {
			for _, p := range region {
				cols[p] = around
			}
		}
	}
	return cols
}
