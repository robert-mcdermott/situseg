// Command situseg segments cells from the transcripts of imaging-based spatial transcriptomics
// (10x Genomics Xenium, NanoString CosMx, Vizgen MERSCOPE, or any table of molecules with places,
// genes and the instrument's cells), by sampling the cells' boundaries from their posterior under
// a model of each cell's expression (after Jones et al., Nat Methods 2025).
//
//	situseg [flags] transcripts.parquet|.csv|.csv.gz|.tsv
//
// It writes the cells, each transcript's cell and its probability, the cell-by-gene counts, the
// cells' outlines and the run's diagnostics to the output folder (situseg-output by default).
package main

import (
	"context"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"os"
	"os/signal"
	"path/filepath"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
	"time"

	"github.com/robert-mcdermott/situseg"
	"github.com/robert-mcdermott/situseg/internal/sysinfo"
)

const usageHead = `situseg segments cells from the transcripts of imaging-based spatial transcriptomics, sampling
the cells' boundaries from their posterior under a model of each cell's expression (after Jones
et al., "Cell simulation as cell segmentation", Nat Methods 2025).

Usage:
  situseg [flags] transcripts.parquet|.csv|.csv.gz|.tsv
  situseg presets      the presets and what each is for
  situseg simulate     write a simulated tissue's transcripts, with the truth
  situseg score        score segmentations of a simulated tissue against its truth
  situseg version

The table's columns are recognized for Xenium (transcripts.parquet, transcripts.csv.gz), CosMx
(tx_file.csv) and MERSCOPE (detected_transcripts.csv); name them with the column flags for any
other table. Every flag may come before or after the table.

Flags:
`

// regions collects --region flags.
type regions []situseg.Region

func (r *regions) String() string { return fmt.Sprint(*r) }

func (r *regions) Set(s string) error {
	parts := strings.Split(s, ",")
	if len(parts) != 4 {
		return errors.New("a region is x0,y0,x1,y1 in µm")
	}
	var v [4]float64
	for i, p := range parts {
		f, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil {
			return fmt.Errorf("a region is x0,y0,x1,y1 in µm: %w", err)
		}
		v[i] = f
	}
	*r = append(*r, situseg.Region{X0: v[0], Y0: v[1], X1: v[2], Y1: v[3]})
	return nil
}

func main() {
	if len(os.Args) > 1 {
		switch os.Args[1] {
		case "presets":
			printPresets()
			return
		case "simulate", "score":
			command := simulateCommand
			if os.Args[1] == "score" {
				command = scoreCommand
			}
			if err := command(os.Args[2:]); err != nil {
				fmt.Fprintln(os.Stderr, "situseg:", err)
				os.Exit(1)
			}
			return
		case "version", "-version", "--version":
			fmt.Println("situseg", version())
			return
		}
	}
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "situseg:", err)
		if f := situseg.Classify(err); f != nil && f.Remedy != "" {
			fmt.Fprintln(os.Stderr, f.Remedy)
		}
		os.Exit(1)
	}
}

func version() string {
	if info, ok := debug.ReadBuildInfo(); ok && info.Main.Version != "" && info.Main.Version != "(devel)" {
		return info.Main.Version
	}
	return "devel"
}

func run(args []string) error {
	fs := flag.NewFlagSet("situseg", flag.ExitOnError)
	fs.Usage = func() {
		fmt.Fprint(fs.Output(), usageHead)
		fs.PrintDefaults()
	}
	// The table.
	xenium := fs.Bool("xenium", false, "the table is Xenium's transcripts.parquet or transcripts.csv.gz")
	cosmx := fs.Bool("cosmx", false, "the table is CosMx's tx_file.csv (places in pixels of 0.12028 µm)")
	merscope := fs.Bool("merscope", false, "the table is MERSCOPE's detected_transcripts.csv")
	colX := fs.String("x", "", "the column of x (µm, or units of -scale)")
	colY := fs.String("y", "", "the column of y")
	colZ := fs.String("z", "", "the column of z (optional)")
	colGene := fs.String("gene", "", "the column of the gene (or control) name")
	colCell := fs.String("cell", "", "the column of the instrument's cell")
	colNucleus := fs.String("nucleus", "", "the column saying whether a transcript is over its cell's nucleus (optional)")
	colQV := fs.String("qv", "", "the column of the quality score (optional)")
	colID := fs.String("id", "", "the column of the transcript id (optional)")
	colGroup := fs.String("cell-group", "", "a column that makes cell ids unique with the cell column (optional)")
	noCell := fs.String("no-cell", "", "comma-separated values of the cell column meaning no cell")
	nuclearValues := fs.String("nuclear-values", "", "comma-separated values of the nucleus column meaning over the nucleus")
	controls := fs.String("controls", "", "a regular expression of the negative controls' names")
	scale := fs.Float64("scale", 0, "multiply x and y by this to make µm")
	zscale := fs.Float64("z-scale", 0, "multiply z by this to make µm")
	minQV := fs.Float64("min-qv", 20, "model gene transcripts of this quality or more (when the table has quality scores)")
	// The start.
	pixelUm := fs.Float64("pixel-um", 0.5, "the segmentation's pixel (µm)")
	radiusUm := fs.Float64("start-radius-um", 2.5, "each starting cell is the pixels within this distance of its transcripts (µm)")
	// The sampler.
	preset := fs.String("preset", situseg.PresetBalanced, "balanced, conservative, flexible, low-memory or roi-pilot (situseg presets)")
	voxelUm := fs.Float64("voxel-um", 0, "the voxel size (µm) of a one-stage schedule (default: the preset's)")
	layers := fs.Int("layers", 0, "the voxel layers in depth (default: the preset's)")
	sweeps := fs.Int("sweeps", 0, "the sweeps (default: the preset's)")
	samples := fs.Int("samples", 0, "the last sweeps summarized (default: the preset's)")
	seed := fs.Uint64("seed", 0, "the seed (default: the preset's); the same seed gives the same result")
	nuclearPrior := fs.Float64("nuclear-prior", 0, "the prior probability that a transcript over a nucleus is its cell's (default: the preset's)")
	perimeter := fs.Float64("perimeter-factor", 0, "each cell's perimeter within this many times the least for its area (default: the preset's)")
	noDrift := fs.Bool("no-drift", false, "keep every transcript where it was observed (no repositioning)")
	seedCheck := fs.Bool("seed-check", false, "run again with the next seed and report the agreement")
	// The area.
	var regs regions
	fs.Var(&regs, "region", "refine only this rectangle, x0,y0,x1,y1 in µm (repeatable)")
	sections := fs.Bool("sections", false, "refine the tissue's separate pieces apart")
	patch := fs.Int("patch-transcripts", 0, "refine in patches of about this many transcripts (default 2,000,000)")
	// The run.
	threads := fs.Int("threads", 0, "threads (0: every core)")
	force := fs.Bool("force", false, "run even when the estimate is beyond the computer's memory")
	out := fs.String("o", "situseg-output", "the output folder")
	parquetOut := fs.Bool("parquet", false, "write cells, transcripts and boundaries as Parquet (default CSV and GeoJSON)")
	labels := fs.Bool("labels", false, "also write the label image (cell_labels.tif)")
	quiet := fs.Bool("quiet", false, "print only warnings and errors")

	// Flags may come before or after the table.
	var positional []string
	for len(args) > 0 {
		fs.Parse(args)
		args = fs.Args()
		if len(args) > 0 {
			positional = append(positional, args[0])
			args = args[1:]
		}
	}
	if len(positional) != 1 {
		fs.Usage()
		os.Exit(2)
	}
	path := positional[0]
	logf := func(level, format string, a ...any) {
		if *quiet && level == "info" {
			return
		}
		prefix := ""
		if level != "info" {
			prefix = level + ": "
		}
		fmt.Fprintf(os.Stderr, prefix+format+"\n", a...)
	}
	began := time.Now()

	// The columns.
	var cols situseg.Columns
	instrument := "custom"
	switch {
	case *xenium:
		cols, instrument = situseg.XeniumColumns(), "xenium"
	case *cosmx:
		cols, instrument = situseg.CosMxColumns(), "cosmx"
	case *merscope:
		cols, instrument = situseg.MERSCOPEColumns(), "merscope"
	case *colX == "":
		header, err := situseg.TableHeader(path)
		if err != nil {
			return err
		}
		var ok bool
		if cols, instrument, ok = situseg.DetectColumns(header); !ok {
			return fmt.Errorf("%s: the columns (%s) are not Xenium's, CosMx's or MERSCOPE's: name them with -x, -y, -gene and -cell", path, strings.Join(header, ", "))
		}
	}
	for _, f := range []struct {
		v   string
		dst *string
	}{{*colX, &cols.X}, {*colY, &cols.Y}, {*colZ, &cols.Z}, {*colGene, &cols.Feature}, {*colCell, &cols.Cell},
		{*colNucleus, &cols.Nuclear}, {*colQV, &cols.QV}, {*colID, &cols.ID}, {*colGroup, &cols.CellGroup}, {*controls, &cols.Controls}} {
		if f.v != "" {
			*f.dst = f.v
		}
	}
	if *noCell != "" {
		cols.NoCell = strings.Split(*noCell, ",")
	}
	if *nuclearValues != "" {
		cols.NuclearValues = strings.Split(*nuclearValues, ",")
	}
	if *scale != 0 {
		cols.Scale = *scale
	}
	if *zscale != 0 {
		cols.ZScale = *zscale
	}

	// The table and the inputs.
	logf("info", "Reading %s (%s columns)", path, instrument)
	table, err := situseg.ReadTable(path, cols)
	if err != nil {
		return err
	}
	controlFeatures := 0
	for _, c := range table.Control {
		if c {
			controlFeatures++
		}
	}
	logf("info", "%d rows: %d features (%d negative controls), %d cells.", table.Rows(), len(table.Features), controlFeatures, len(table.CellIDs))
	if table.Nuclear == nil {
		logf("warning", "The table does not say which transcripts lie over a nucleus: every transcript the instrument gave a cell has it as its prior cell.")
	}
	readSeconds := time.Since(began).Seconds()
	workers := *threads
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	start := situseg.StartOptions{MinQV: *minQV, PixelUm: *pixelUm, RadiusUm: *radiusUm, Workers: workers}
	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt)
	defer stop()
	in, idx, err := situseg.InputsFromTable(ctx, table, start)
	if err != nil {
		return err
	}
	logf("info", "%d gene transcripts of quality %g or more to model, %d genes, %d cells to start from.", len(idx.Rows), *minQV, in.Genes, in.Cells)

	// The spec.
	spec, err := situseg.SpecOf(*preset)
	if err != nil {
		return err
	}
	o := &spec.Options
	if *voxelUm > 0 || *layers > 0 || *sweeps > 0 {
		st := o.Stages[len(o.Stages)-1]
		if *voxelUm > 0 {
			st.VoxelUm = *voxelUm
		}
		if *layers > 0 {
			st.Layers = *layers
		}
		if *sweeps > 0 {
			st.Sweeps = *sweeps
			o.Samples = min(o.Samples, *sweeps)
		}
		o.Stages = []situseg.Stage{st}
	}
	if *samples > 0 {
		o.Samples = *samples
	}
	if *seed > 0 {
		o.Seed = *seed
	}
	if *nuclearPrior > 0 {
		o.Model.NuclearPrior = *nuclearPrior
	}
	if *perimeter > 0 {
		o.Sampler.PerimeterFactor = *perimeter
	}
	if *noDrift {
		o.Sampler.Reposition = false
	}
	spec.SeedCheck = spec.SeedCheck || *seedCheck
	switch {
	case len(regs) > 0:
		spec.Area = situseg.Area{Mode: situseg.AreaRegions, Regions: regs}
	case *sections:
		spec.Area.Mode = situseg.AreaSections
	}
	spec.Area.PatchTranscripts = *patch
	spec.Resolve()
	if spec.Preset == situseg.PresetExpert {
		logf("info", "Settings changed from the %s preset: the run is Expert.", *preset)
	}

	// The run.
	memory, _ := sysinfo.Memory()
	plan := situseg.Plan(in, idx, spec, workers, memory)
	lastProgress := time.Time{}
	job := situseg.Job{Workers: workers, Memory: memory, Force: *force, Log: logf,
		Progress: func(f float64, message string) {
			if *quiet || time.Since(lastProgress) < 2*time.Second && f < 1 {
				return
			}
			lastProgress = time.Now()
			fmt.Fprintf(os.Stderr, "%5.1f%%  %s\n", 100*f, message)
		}}
	sampled := time.Now()
	seg, err := situseg.Segment(ctx, table, in, idx, spec, job)
	if err != nil {
		return err
	}
	sampleSeconds := time.Since(sampled).Seconds()
	written := time.Now()
	if err := seg.WriteFiles(*out, table, situseg.WriteOptions{Parquet: *parquetOut, Labels: *labels, CreatedBy: "situseg " + version()}); err != nil {
		return err
	}
	record := map[string]any{
		"software":    "situseg " + version(),
		"method":      situseg.Method + " " + situseg.MethodVersion,
		"citation":    situseg.Citation,
		"command":     strings.Join(os.Args, " "),
		"input":       map[string]any{"path": path, "instrument": instrument, "columns": cols, "rows": table.Rows()},
		"start":       start,
		"spec":        spec,
		"recorded":    seg.Recorded,
		"preflight":   plan,
		"threads":     workers,
		"seconds":     map[string]float64{"read": readSeconds, "segment": sampleSeconds, "write": time.Since(written).Seconds(), "total": time.Since(began).Seconds()},
		"diagnostics": seg.Diagnostics.Summary,
	}
	raw, _ := json.MarshalIndent(record, "", " ")
	if err := os.WriteFile(filepath.Join(*out, "run.json"), append(raw, '\n'), 0o644); err != nil {
		return err
	}
	s := seg.Diagnostics.Summary
	assigned := 0
	for _, c := range seg.Assignment {
		if c >= 0 {
			assigned++
		}
	}
	logf("info", "Wrote %s: %d cells, %.1f%% of rows in a cell, in %s.", *out, len(seg.Cells), 100*float64(assigned)/float64(max(1, len(seg.Assignment))), situseg.Duration(time.Since(began).Seconds()))
	if s.GewekeJudged && (s.GewekeZ < -2 || s.GewekeZ > 2) {
		logf("info", "The log-likelihood was still changing over the summarized sweeps (Geweke's z %.1f; within ±2 once settled): more sweeps (-sweeps) may sharpen the summary.", s.GewekeZ)
	}
	return nil
}

func printPresets() {
	for _, p := range situseg.Presets() {
		fmt.Printf("%s (%s)\n", p.ID, p.Name)
		fmt.Printf("  %s\n", p.About)
		if p.Changes != "" {
			fmt.Printf("  Changes: %s\n", p.Changes)
		}
		if p.Warning != "" {
			fmt.Printf("  Warning: %s\n", p.Warning)
		}
		if p.Measured != "" {
			fmt.Printf("  Measured: %s\n", p.Measured)
		}
		fmt.Println()
	}
}
