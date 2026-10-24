package situseg

// A segmentation from end to end: its spec (the settings, what it covers, the preset they came
// from), checked before it runs against the inputs and the computer, then run (each region
// sampled, the results stitched into one segmentation of the whole table).

import (
	"context"
	"errors"
	"fmt"
	"math"
	"runtime"
	"runtime/debug"
	"sync"
)

// Spec is a segmentation's settings: what a preset holds and a run records.
type Spec struct {
	// Preset is the preset the options came from ("" for none; Expert once they differ from it in
	// more than the seed).
	Preset  string  `json:"preset,omitempty"`
	Options Options `json:"options"`
	Area    Area    `json:"area"`
	// SeedCheck runs the regions again with the next seed and reports the two runs' agreement.
	SeedCheck bool `json:"seed_check,omitempty"`
	// GCPercent is the garbage collector's target while it runs (0: 50; ConserveMemory): it trades
	// memory against time, and never changes the result.
	GCPercent int `json:"gc_percent,omitempty"`
}

// DefaultSpec is Balanced on the whole tissue.
func DefaultSpec() Spec {
	return Spec{Preset: PresetBalanced, Options: DefaultOptions(), Area: Area{Mode: AreaWhole}}
}

// SpecOf is a preset's spec: its options, area mode and seed check.
func SpecOf(preset string) (Spec, error) {
	p, ok := PresetByID(preset)
	if !ok {
		return Spec{}, fmt.Errorf("situseg: no preset %q", preset)
	}
	s := DefaultSpec()
	s.Preset, s.Options, s.SeedCheck, s.GCPercent = p.ID, p.Options, p.SeedCheck, p.GCPercent
	if p.Area != "" {
		s.Area.Mode = p.Area
	}
	return s, nil
}

// Resolve names the preset the spec's settings are: its own while they are the preset's but for
// the seed and the area, Expert otherwise.
func (s *Spec) Resolve() {
	if s.Preset == "" || s.Preset == PresetExpert {
		return
	}
	p, ok := PresetByID(s.Preset)
	if !ok {
		return
	}
	o := s.Options
	o.Seed, o.Workers = p.Options.Seed, 0
	if !sameSettings(o, p.Options) || s.SeedCheck != p.SeedCheck {
		s.Preset = PresetExpert
	}
}

func sameSettings(a, b Options) bool {
	if len(a.Stages) != len(b.Stages) {
		return false
	}
	for i := range a.Stages {
		if a.Stages[i] != b.Stages[i] {
			return false
		}
	}
	return a.Model == b.Model && a.Sampler == b.Sampler && a.BlockUm == b.BlockUm && a.Samples == b.Samples && a.Seed == b.Seed
}

// Check refuses a spec a segmentation cannot run, saying why in a sentence.
func (s Spec) Check() error {
	if err := s.Options.Check(); err != nil {
		return err
	}
	if s.GCPercent < 0 || s.GCPercent > 400 {
		return errors.New("situseg: the memory target is from 1 to 400 (0: the default)")
	}
	switch s.Area.Mode {
	case "", AreaWhole, AreaSections:
	case AreaRegions:
		if len(s.Area.Regions) == 0 {
			return errors.New("situseg: no region was chosen")
		}
		for _, r := range s.Area.Regions {
			if r.Empty() || math.IsInf(r.Area(), 0) || math.IsNaN(r.Area()) {
				return errors.New("situseg: a region needs its lower corner below and left of its upper one")
			}
		}
	default:
		return errors.New("situseg: the area is whole, sections or regions")
	}
	if s.Area.GapUm < 0 {
		return errors.New("situseg: the gap between sections is positive")
	}
	if s.Area.PatchTranscripts != 0 && (s.Area.PatchTranscripts < 10_000 || s.Area.PatchTranscripts > 100_000_000) {
		return errors.New("situseg: a patch holds 10,000 to 100,000,000 transcripts")
	}
	if s.Preset != "" {
		p, ok := PresetByID(s.Preset)
		if !ok {
			return fmt.Errorf("situseg: no preset %q", s.Preset)
		}
		if p.Area == AreaRegions && s.Area.Mode != AreaRegions {
			return fmt.Errorf("situseg: the %s preset segments chosen regions: choose one", p.Name)
		}
	}
	return nil
}

// Finding is one of a preflight's problems (which stop the run) or risks (which it reports).
type Finding struct {
	Class   string `json:"class"`
	Message string `json:"message"`
	Remedy  string `json:"remedy,omitempty"`
}

// Preflight is what a segmentation will model and cover, and what it will cost, before it runs.
type Preflight struct {
	Spec Spec `json:"spec"`
	// The table's rows; those of negative controls, below the quality and without a place; the
	// gene transcripts modeled (before the regions) and those with a prior cell.
	Rows        int `json:"rows"`
	Controls    int `json:"controls"`
	LowQuality  int `json:"belowQuality"`
	Outside     int `json:"outside"`
	Transcripts int `json:"transcripts"`
	WithPrior   int `json:"withPrior"`
	Genes       int `json:"genes"`
	// The cells it starts from, and those with a nucleus (−1: not known).
	Cells   int      `json:"cells"`
	Nuclei  int      `json:"cellsWithNucleus"`
	Regions []Region `json:"regions"`
	// Patches are the runs, when a region holds too many transcripts for one (Units): each patch's
	// square and the region whose cells it refines.
	Patches []Unit `json:"patches,omitempty"`
	// Estimate is its cost; with a second seed, Seconds doubles.
	Estimate Estimate `json:"estimate"`
	// MemoryBytes is the computer's memory (0: unknown).
	MemoryBytes int64     `json:"memoryBytes"`
	Problems    []Finding `json:"problems"`
	Risks       []Finding `json:"risks"`
}

// Thresholds of the preflight: a run whose estimated memory passes memoryStop of the computer's
// is refused (the estimate is within about 15%, and the computer runs other programs), past
// memoryWarn it is a risk.
const (
	memoryStop = 0.85
	memoryWarn = 0.5
)

// Plan is a spec's preflight on a table's inputs, with the computer's memory (0: unknown) and the
// threads it runs on.
func Plan(in *Inputs, idx *Index, spec Spec, cores int, memory int64) Preflight {
	spec.Resolve()
	p := Preflight{Spec: spec, Rows: idx.Total, Controls: idx.Controls, LowQuality: idx.LowQuality, Outside: idx.Outside,
		Transcripts: len(idx.Rows), WithPrior: idx.WithPrior, Genes: in.Genes, Cells: in.Cells, Nuclei: -1, MemoryBytes: memory,
		Problems: []Finding{}, Risks: []Finding{}}
	if in.Nucleus != nil {
		p.Nuclei = 0
		for _, has := range in.Nucleus {
			if has {
				p.Nuclei++
			}
		}
	}
	problem := func(class, remedy, format string, args ...any) {
		p.Problems = append(p.Problems, Finding{Class: class, Message: fmt.Sprintf(format, args...), Remedy: remedy})
	}
	risk := func(class, format string, args ...any) {
		p.Risks = append(p.Risks, Finding{Class: class, Message: fmt.Sprintf(format, args...)})
	}
	if err := spec.Check(); err != nil {
		problem("settings", "", "%s.", plain(err))
		return p
	}
	regions, err := in.Regions(spec.Area)
	if err != nil {
		f := Classify(err)
		if f == nil {
			f = &Failed{Class: ClassRegion, Remedy: remedyRegion, Err: err}
		}
		problem(f.Class, f.Remedy, "%s.", plain(err))
		return p
	}
	p.Regions = regions
	units, _ := in.Units(spec.Area)
	for _, u := range units {
		if u.Own != u.Free {
			p.Patches = units
			break
		}
	}
	p.Estimate = EstimateCost(in, idx, units, spec.Options, cores, spec.GCPercent)
	if spec.SeedCheck {
		p.Estimate.Seconds *= 2
	}
	for i, r := range regions {
		cells := 0
		for c := 0; c < in.Cells; c++ {
			if r.Contains(in.SeedX[c], in.SeedY[c]) {
				cells++
			}
		}
		if cells == 0 {
			problem(ClassRegion, remedyRegion, "Region %d, (%.0f, %.0f)–(%.0f, %.0f) µm, holds no cell's seed.", i+1, r.X0, r.Y0, r.X1, r.Y1)
		}
	}
	e := p.Estimate
	switch {
	case memory > 0 && float64(e.PeakBytes) > memoryStop*float64(memory):
		problem(ClassMemory, remedyMemory, "The segmentation needs about %s of memory; the computer has %s.", Memory(e.PeakBytes), Memory(memory))
	case memory > 0 && float64(e.PeakBytes) > memoryWarn*float64(memory):
		risk(ClassMemory, "The segmentation needs about %s of memory, more than half the computer's %s: close other programs while it runs.", Memory(e.PeakBytes), Memory(memory))
	}
	if spec.Options.Samples < minGeweke {
		risk("convergence", "Only %d sweeps are summarized: convergence is judged on %d or more, so this run's will not be.", spec.Options.Samples, minGeweke)
	}
	if in.MinQV < 20 && idx.LowQuality+len(idx.Rows) > 0 {
		risk("quality", "Transcripts below Q20 are modeled (%g or more): more of them were misread.", in.MinQV)
	}
	if p.Nuclei >= 0 && p.Cells > 0 && float64(p.Cells-p.Nuclei) > 0.1*float64(p.Cells) {
		risk("prior", "%d of the %d cells have no transcript over a nucleus: no nuclear prior holds them.", p.Cells-p.Nuclei, p.Cells)
	}
	if e.Cells > 0 && float64(e.Transcripts)/float64(e.Cells+e.Held) < 20 {
		risk("transcripts", "About %.0f transcripts per cell: with so few, boundaries come more from the prior than from expression.", float64(e.Transcripts)/float64(e.Cells+e.Held))
	}
	if spec.Preset == PresetExpert {
		if x, ok := PresetByID(PresetExpert); ok {
			risk("settings", "%s", x.Warning)
		}
	}
	return p
}

// plain is an error's message without the package's prefix, as a sentence.
func plain(err error) string {
	m := err.Error()
	if len(m) > 9 && m[:9] == "situseg: " {
		m = m[9:]
	}
	if m != "" && m[0] >= 'a' && m[0] <= 'z' {
		m = string(m[0]-'a'+'A') + m[1:]
	}
	return m
}

// Memory is an amount of memory as computers state theirs, in binary units (a GB of 2³⁰ bytes).
func Memory(n int64) string {
	switch {
	case n >= 1<<30:
		return fmt.Sprintf("%.1f GB", float64(n)/(1<<30))
	case n >= 1<<20:
		return fmt.Sprintf("%.0f MB", float64(n)/(1<<20))
	case n >= 1<<10:
		return fmt.Sprintf("%.0f kB", float64(n)/(1<<10))
	}
	return fmt.Sprintf("%d bytes", n)
}

// Duration is a time in a unit people read.
func Duration(s float64) string {
	switch {
	case s >= 5400:
		return fmt.Sprintf("%.1f hours", s/3600)
	case s >= 90:
		return fmt.Sprintf("%.0f minutes", s/60)
	}
	if n := math.Round(math.Max(1, s)); n != 1 {
		return fmt.Sprintf("%.0f seconds", n)
	}
	return "1 second"
}

// Job is how a segmentation runs: on how many threads (0: every core), whether past the memory the
// computer has (0: unknown), and where it reports.
type Job struct {
	Workers  int
	Memory   int64
	Force    bool
	Progress func(float64, string)
	Log      func(level, format string, args ...any)
}

// Segment segments a table as a spec says: its preflight checked against the computer's memory,
// each region sampled (twice with a seed check), the sampler's state let go and the results
// stitched. Its failures are classified (Classify).
func Segment(ctx context.Context, t *Table, in *Inputs, idx *Index, spec Spec, job Job) (*Segmentation, error) {
	progress, log := job.Progress, job.Log
	if progress == nil {
		progress = func(float64, string) {}
	}
	if log == nil {
		log = func(string, string, ...any) {}
	}
	spec.Resolve()
	if err := spec.Check(); err != nil {
		return nil, err
	}
	defer ConserveMemory(spec.GCPercent)()
	cores := job.Workers
	if cores <= 0 {
		cores = runtime.NumCPU()
	}
	pre := Plan(in, idx, spec, cores, job.Memory)
	for _, f := range pre.Problems {
		if f.Class == ClassMemory && job.Force {
			log("warning", "%s Running anyway, as asked.", f.Message)
			continue
		}
		return nil, &Failed{Class: f.Class, Remedy: f.Remedy, Err: errors.New(f.Message)}
	}
	for _, f := range pre.Risks {
		log("warning", "%s", f.Message)
	}
	e := pre.Estimate
	runsOf := "region(s)"
	if len(pre.Patches) > 0 {
		runsOf = "patches"
	}
	log("info", "%d %s: %d cells refined, %d held about them, %d transcripts modeled; about %s of memory and %s on %d threads.",
		e.Regions, runsOf, e.Cells, e.Held, e.Transcripts, Memory(e.PeakBytes), Duration(e.Seconds), cores)
	opt := spec.Options
	opt.Workers = cores
	runs := 1
	if spec.SeedCheck {
		runs = 2
	}
	sample := func(o Options, k int) ([]*Result, error) {
		return RunArea(ctx, in, spec.Area, o, func(f float64, message string) {
			if runs > 1 {
				message = fmt.Sprintf("Seed %d: %s", o.Seed, message)
			}
			progress(0.9*(float64(k)+f)/float64(runs), message)
		})
	}
	results, err := sample(opt, 0)
	if err != nil {
		return nil, err
	}
	rec := Recorded{Area: spec.Area, Preset: spec.Preset}
	if spec.Preset != "" {
		rec.PresetsVersion = PresetsVersion
	}
	if spec.SeedCheck {
		second := opt
		second.Seed++
		if rec.Second, err = sample(second, 1); err != nil {
			return nil, err
		}
		rec.SecondSeed = second.Seed
	}
	for _, res := range append(append([]*Result(nil), results...), rec.Second...) {
		res.Release()
	}
	runtime.GC() // the sampler's state goes before the results are stitched
	left := 0
	for _, res := range results {
		left += len(res.LeftOut)
	}
	if left > 0 {
		log("warning", "%d cells could not be kept: their seeds share a voxel with another's at the finest stage.", left)
	}
	progress(0.9, "Stitching the segmentation")
	seg, err := Stitch(t, in, idx, results, opt, rec)
	if err != nil {
		return nil, err
	}
	s := seg.Diagnostics.Summary
	log("info", "%d cells; %.1f%% of modeled transcripts in the same cell over both halves of the last %d sweeps; mean assignment probability %.3f; Geweke's z %.1f.",
		len(seg.Cells), 100*s.SplitHalf, s.Samples, s.MeanProbability, s.GewekeZ)
	if s.SeedAgreement != nil {
		log("info", "Seeds %d and %d put %.1f%% of the transcripts in refined cells in the same cell.", opt.Seed, rec.SecondSeed, 100*(*s.SeedAgreement))
	}
	progress(1, "Done")
	return seg, nil
}

// Release lets the sampler's state and the model's parameters go, keeping the lattice and the
// hyperparameters the outputs record: the results are stitched without them.
func (r *Result) Release() {
	if r.State != nil {
		r.State = &State{Lat: r.State.Lat}
	}
	if r.Model != nil {
		r.Model = &Model{H: r.Model.H}
	}
}

// defaultGC is a segmentation's target for the garbage collector.
const defaultGC = 50

var gcHold struct {
	sync.Mutex
	held   []int
	before int
}

// ConserveMemory has the garbage collector keep the heap within percent more than its live size
// (0: 50; the Go default is 100) until the answered function is called: on a standard simulated
// tissue, at 50, a segmentation's peak memory fell by a quarter and its time rose by about 1%. A
// lower setting of the user's own (GOGC) is kept; while several hold it, the lowest holds.
func ConserveMemory(percent int) (restore func()) {
	if percent <= 0 {
		percent = defaultGC
	}
	gcHold.Lock()
	defer gcHold.Unlock()
	if len(gcHold.held) == 0 {
		gcHold.before = debug.SetGCPercent(-1)
		debug.SetGCPercent(gcHold.before)
	}
	gcHold.held = append(gcHold.held, percent)
	apply := func() {
		target := gcHold.before
		for _, p := range gcHold.held {
			if target < 0 || p < target {
				target = p
			}
		}
		if gcHold.before < 0 {
			target = gcHold.before
		}
		debug.SetGCPercent(target)
	}
	apply()
	var once sync.Once
	return func() {
		once.Do(func() {
			gcHold.Lock()
			defer gcHold.Unlock()
			for i, p := range gcHold.held {
				if p == percent {
					gcHold.held = append(gcHold.held[:i], gcHold.held[i+1:]...)
					break
				}
			}
			if len(gcHold.held) == 0 {
				debug.SetGCPercent(gcHold.before)
				return
			}
			apply()
		})
	}
}
