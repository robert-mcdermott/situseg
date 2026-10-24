// Package simulate makes tissues whose every cell and molecule is known, as a table of transcripts
// an instrument would write with the truth beside it, so that a segmentation can be scored
// against the truth.
//
// A tissue is a field of cells of several types, in domains of different composition. Each cell
// grows from its nucleus to a size, elongation and irregular outline of its type, competing for
// space with its neighbours; some cells' nuclei lie outside the section and are not seen. Each
// cell expresses the panel's genes with negative-binomial counts around its type's profile, in
// proportion to its area; its molecules lie within it, in the section's depth, and some are
// observed displaced (diffused), a few far. Background molecules, ambient ones of every gene and
// false decodes of every codeword (genes and negative controls alike), fall anywhere. Every
// molecule has a quality score, lower for background. The instrument segments by expanding the
// nuclei it sees, with a distance of its own (instrument.go), and gives each transcript the cell
// holding the pixel it is observed in.
//
// The simulator shares no code with the segmentation it scores. The same parameters and seed
// give the same tissue, molecule for molecule, on any processor (internal/random,
// internal/dmath).
package simulate

import (
	"context"
	"errors"
	"fmt"
	"math"
	"runtime"
	"sync"

	"github.com/robert-mcdermott/situseg/internal/random"
)

// CellType is a kind of cell: its share of each domain comes from Params.Domains.
type CellType struct {
	Name string `json:"name"`
	// AreaUm2 is the median cell area (µm²) and AreaSpread the log-normal σ around it.
	AreaUm2    float64 `json:"area_um2"`
	AreaSpread float64 `json:"area_spread"`
	// Elongation is the largest ratio of a cell's long axis to its short one (drawn from 1 up).
	Elongation float64 `json:"elongation"`
	// NucleusShare is the nucleus's share of the cell's radius.
	NucleusShare float64 `json:"nucleus_share"`
	// Density is the expected gene transcripts of a cell per 100 µm² of its area.
	Density float64 `json:"density_per_100um2"`
	// Markers is how many genes mark the type, each expressed Fold times its base level.
	Markers int     `json:"markers"`
	Fold    float64 `json:"marker_fold"`
	// Dispersion is the negative binomial's r (variance m + m²/r).
	Dispersion float64 `json:"dispersion"`
	// HeightUm is the cell's extent in z.
	HeightUm float64 `json:"height_um"`
}

// Params are a simulation's every setting.
type Params struct {
	Seed     uint64  `json:"seed"`
	WidthUm  float64 `json:"width_um"`
	HeightUm float64 `json:"height_um"`
	DepthUm  float64 `json:"depth_um"`
	PixelUm  float64 `json:"pixel_um"`
	// Fill is the share of the field the cells' nominal areas cover (gaps appear between them).
	Fill float64 `json:"fill"`
	// Domains is the number of tissue domains, each a Voronoi region with its own mix of types,
	// one of which dominates it.
	Domains int        `json:"domains"`
	Types   []CellType `json:"types"`
	// Genes is the panel's genes (the types' markers among them); Controls the negative control
	// probes and Codewords the negative control codewords.
	Genes     int `json:"genes"`
	Controls  int `json:"control_probes"`
	Codewords int `json:"control_codewords"`
	// HiddenNuclei is the share of cells whose nucleus lies outside the section.
	HiddenNuclei float64 `json:"hidden_nuclei"`
	// NuclearEnrichment is how much denser a cell's molecules are over its nucleus.
	NuclearEnrichment float64 `json:"nuclear_enrichment"`
	// Diffusion: a share of molecules (WideShare) are displaced by a wide normal of σ WideUm in
	// x and y, the rest by a narrow one of σ NarrowUm; z by ZUm.
	WideShare float64 `json:"diffusion_wide_share"`
	WideUm    float64 `json:"diffusion_wide_um"`
	NarrowUm  float64 `json:"diffusion_narrow_um"`
	ZUm       float64 `json:"diffusion_z_um"`
	// The background, molecules from no cell, of two kinds:
	//   - Ambient: gene molecules per 100 µm², shared among the genes in proportion to their mean
	//     expression (RNA loose in the tissue);
	//   - FalseDecodes: molecules per 100 µm² on every codeword alike, genes and controls (decoding
	//     errors and probes bound off target), which is what the negative controls measure.
	Ambient      float64 `json:"ambient_per_100um2"`
	FalseDecodes float64 `json:"false_decodes_per_100um2"`
	// LowQuality is the share of cells' molecules with a quality below 20; background's is
	// BackgroundLowQuality.
	LowQuality           float64 `json:"low_quality"`
	BackgroundLowQuality float64 `json:"background_low_quality"`
	// ExpansionUm is how far the simulated instrument grows the nuclei it sees into cells.
	ExpansionUm float64 `json:"instrument_expansion_um"`
}

// seed is a stream's seed.
func (p Params) seed(uint64) uint64 { return p.Seed }

// MaxPixels is the largest image a simulation makes (its memory grows with it).
const MaxPixels = 50e6

// Presets are named starting points: small (fast, for tests and trials) and standard.
var Presets = map[string]Params{"small": small(), "standard": standard()}

func standardTypes() []CellType {
	return []CellType{
		{Name: "epithelial-like", AreaUm2: 110, AreaSpread: 0.25, Elongation: 1.4, NucleusShare: 0.5, Density: 280, Markers: 8, Fold: 12, Dispersion: 4, HeightUm: 8},
		{Name: "fibroblast-like", AreaUm2: 90, AreaSpread: 0.3, Elongation: 3.2, NucleusShare: 0.42, Density: 160, Markers: 8, Fold: 12, Dispersion: 3, HeightUm: 6},
		{Name: "immune-like", AreaUm2: 45, AreaSpread: 0.2, Elongation: 1.15, NucleusShare: 0.62, Density: 260, Markers: 8, Fold: 14, Dispersion: 5, HeightUm: 6},
		{Name: "endothelial-like", AreaUm2: 70, AreaSpread: 0.25, Elongation: 4, NucleusShare: 0.4, Density: 150, Markers: 6, Fold: 12, Dispersion: 3, HeightUm: 5},
	}
}

func standard() Params {
	return Params{
		Seed: 20261014, WidthUm: 500, HeightUm: 500, DepthUm: 10, PixelUm: 0.2125, Fill: 0.85,
		Domains: 6, Types: standardTypes(), Genes: 120, Controls: 20, Codewords: 41,
		HiddenNuclei: 0.04, NuclearEnrichment: 1.5,
		WideShare: 0.2, WideUm: 3, NarrowUm: 0.5, ZUm: 0.5,
		Ambient: 3, FalseDecodes: 0.025, LowQuality: 0.08, BackgroundLowQuality: 0.35,
		ExpansionUm: 5,
	}
}

func small() Params {
	p := standard()
	p.WidthUm, p.HeightUm, p.Domains = 160, 120, 3
	return p
}

// Check refuses parameters that cannot make a tissue, saying which.
func (p *Params) Check() error {
	var bad []string
	need := func(ok bool, what string) {
		if !ok {
			bad = append(bad, what)
		}
	}
	need(p.WidthUm >= 20 && p.HeightUm >= 20 && p.WidthUm <= 20000 && p.HeightUm <= 20000, "the field is 20 µm to 20 mm a side")
	need(p.DepthUm > 0 && p.DepthUm <= 100, "the depth is above 0 and at most 100 µm")
	need(p.PixelUm >= 0.05 && p.PixelUm <= 2, "the pixel is 0.05 to 2 µm")
	need(p.Fill > 0 && p.Fill <= 1.5, "the fill is above 0 and at most 1.5")
	need(p.Domains >= 1 && p.Domains <= 1000, "1 to 1,000 domains")
	need(len(p.Types) >= 1 && len(p.Types) <= 50, "1 to 50 cell types")
	markers := 0
	for _, t := range p.Types {
		markers += t.Markers
		need(t.Name != "", "every cell type has a name")
		need(t.AreaUm2 >= 4 && t.AreaUm2 <= 10000, fmt.Sprintf("%s: an area of 4 to 10,000 µm²", t.Name))
		need(t.AreaSpread >= 0 && t.AreaSpread <= 2, fmt.Sprintf("%s: an area spread of 0 to 2", t.Name))
		need(t.Elongation >= 1 && t.Elongation <= 10, fmt.Sprintf("%s: an elongation of 1 to 10", t.Name))
		need(t.NucleusShare > 0 && t.NucleusShare < 1, fmt.Sprintf("%s: a nucleus share between 0 and 1", t.Name))
		need(t.Density >= 0 && t.Density <= 10000, fmt.Sprintf("%s: a density of 0 to 10,000", t.Name))
		need(t.Markers >= 0 && t.Fold >= 1, fmt.Sprintf("%s: markers and a fold of 1 or more", t.Name))
		need(t.Dispersion > 0, fmt.Sprintf("%s: a positive dispersion", t.Name))
		need(t.HeightUm > 0, fmt.Sprintf("%s: a positive height", t.Name))
	}
	need(p.Genes >= markers && p.Genes >= 1 && p.Genes <= 20000, "at least as many genes as markers, at most 20,000")
	need(p.Controls >= 0 && p.Codewords >= 0 && p.Controls+p.Codewords <= 1000, "0 to 1,000 controls")
	need(p.HiddenNuclei >= 0 && p.HiddenNuclei < 1, "a hidden-nucleus share below 1")
	need(p.NuclearEnrichment >= 0.1 && p.NuclearEnrichment <= 20, "a nuclear enrichment of 0.1 to 20")
	need(p.WideShare >= 0 && p.WideShare <= 1 && p.WideUm >= 0 && p.NarrowUm >= 0 && p.ZUm >= 0, "diffusion shares and widths of 0 or more")
	need(p.Ambient >= 0 && p.FalseDecodes >= 0, "background densities of 0 or more")
	need(p.LowQuality >= 0 && p.LowQuality <= 1 && p.BackgroundLowQuality >= 0 && p.BackgroundLowQuality <= 1, "quality shares of 0 to 1")
	need(p.ExpansionUm >= 0 && p.ExpansionUm <= 50, "an instrument expansion of 0 to 50 µm")
	if pixels := math.Ceil(p.WidthUm/p.PixelUm) * math.Ceil(p.HeightUm/p.PixelUm); pixels > MaxPixels {
		bad = append(bad, fmt.Sprintf("at most %d million pixels (%.0f million asked; about 1.5 mm a side at 0.2125 µm, which needs about 4 GB of memory)", int(MaxPixels/1e6), pixels/1e6))
	}
	if len(bad) > 0 {
		return fmt.Errorf("the simulation's parameters: %v", bad)
	}
	return nil
}

var errStopped = errors.New("the simulation was stopped")

// Tissue is a simulated tissue's transcripts, as an instrument writes them, with the truth.
type Tissue struct {
	Params Params
	// The table: each molecule's observed place (µm), feature (an index into Features; Control
	// marks the negative controls), quality, transcript id, the instrument's cell (an index into
	// CellIDs, −1 none) and whether it lies over a nucleus the instrument sees.
	X, Y, Z  []float32
	Feature  []uint32
	Features []string
	Control  []bool
	QV       []float32
	ID       []uint64
	Cell     []int32
	CellIDs  []string
	Nuclear  []bool
	// The truth: each molecule's true cell (an index into TrueIDs, −1 none), true place and
	// origin (FromCell, FromBackground, FromControlProbe, FromControlCodeword); each true cell's
	// type (an index into Params.Types) and whether its nucleus is in the section.
	TrueCell       []int32
	TrueIDs        []string
	TrueX, TrueY   []float32
	Origin         []uint8
	TrueType       []int32
	TrueHasNucleus []bool
}

// Simulate makes the tissue p describes.
func Simulate(ctx context.Context, p Params) (*Tissue, error) {
	if err := p.Check(); err != nil {
		return nil, err
	}
	stop := func() bool { return ctx.Err() != nil }
	t := newTissue(p)
	layout := random.New(p.seed(streamLayout), streamLayout)
	t.makeDomains(layout)
	t.placeCells(layout)
	if err := t.rasterize(stop); err != nil {
		return nil, err
	}
	t.placeNuclei(layout)
	if len(t.cells) == 0 {
		return nil, errors.New("the simulation made no cells: raise the fill or the cells' areas")
	}
	pn := makePanel(p, random.New(p.seed(streamPanel), streamPanel), t.typeShares())
	mols, err := t.molecules(pn, stop)
	if err != nil {
		return nil, err
	}
	in := t.segmentLikeTheInstrument()
	if len(in.cells) == 0 {
		return nil, errors.New("no nucleus of the simulated tissue is in the section: lower the hidden-nucleus share")
	}
	out := &Tissue{Params: p, Features: pn.names}
	for f := range pn.names {
		out.Control = append(out.Control, f >= pn.genes)
	}
	for k := range in.cells {
		out.CellIDs = append(out.CellIDs, fmt.Sprintf("cell-%d", k+1))
	}
	for k := range t.cells {
		out.TrueIDs = append(out.TrueIDs, fmt.Sprintf("true-%d", k+1))
		out.TrueType = append(out.TrueType, int32(t.cells[k].typ))
		out.TrueHasNucleus = append(out.TrueHasNucleus, t.cells[k].nucleus)
	}
	n := len(mols)
	out.X, out.Y, out.Z = make([]float32, n), make([]float32, n), make([]float32, n)
	out.Feature, out.QV, out.ID = make([]uint32, n), make([]float32, n), make([]uint64, n)
	out.Cell, out.Nuclear = make([]int32, n), make([]bool, n)
	out.TrueCell, out.TrueX, out.TrueY, out.Origin = make([]int32, n), make([]float32, n), make([]float32, n), make([]uint8, n)
	ps := p.PixelUm
	for i, m := range mols {
		px := min(t.w-1, max(0, int(float64(m.x)/ps)))
		py := min(t.h-1, max(0, int(float64(m.y)/ps)))
		pixel := py*t.w + px
		out.X[i], out.Y[i], out.Z[i] = m.x, m.y, m.z
		out.Feature[i], out.QV[i], out.ID[i] = uint32(m.feature), m.qv, uint64(i+1)
		out.Cell[i] = int32(in.labels[pixel]) - 1
		out.Nuclear[i] = in.nuclei[pixel] != 0
		out.TrueCell[i] = int32(m.cell) - 1
		out.TrueX[i], out.TrueY[i], out.Origin[i] = m.tx, m.ty, m.origin
	}
	return out, nil
}

// typeShares are the types' shares of the cells' expected transcripts, which weigh the
// background's genes.
func (t *tissue) typeShares() []float64 {
	p := t.p
	share := make([]float64, len(p.Types))
	var total float64
	for k := range t.cells {
		c := &t.cells[k]
		e := mul(p.Types[c.typ].Density, float64(len(c.pixels)))
		share[c.typ] += e
		total += e
	}
	for k := range share {
		if total > 0 {
			share[k] /= total
		} else {
			share[k] = 1 / float64(len(share))
		}
	}
	return share
}

// molecules draws every cell's molecules and every background tile's, in parallel, each from its
// own stream, merged in order: the same for any number of processors.
func (t *tissue) molecules(pn *panel, stop func() bool) ([]molecule, error) {
	cells, tiles := len(t.cells), t.tiles()
	parts := make([][]molecule, cells+tiles)
	jobs := make(chan int, 256)
	var wg sync.WaitGroup
	for w := 0; w < runtime.GOMAXPROCS(0); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for j := range jobs {
				if j < cells {
					parts[j] = t.cellMolecules(j, pn, nil)
				} else {
					parts[j] = t.background(j-cells, pn, nil)
				}
			}
		}()
	}
	stopped := false
	for j := range parts {
		if j%256 == 0 && stop() {
			stopped = true
			break
		}
		jobs <- j
	}
	close(jobs)
	wg.Wait()
	if stopped {
		return nil, errStopped
	}
	n := 0
	for _, part := range parts {
		n += len(part)
	}
	out := make([]molecule, 0, n)
	for _, part := range parts {
		out = append(out, part...)
	}
	return out, nil
}
