# The method

situseg segments cells after Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR, Newell EW.
Cell simulation as cell segmentation. *Nat Methods* 2025;22(6):1331–1342.
doi:10.1038/s41592-025-02697-0.

A cell is a region of space whose molecules are a Poisson process with its own rate for each
gene, over a background with a rate for each gene everywhere. A segmentation is the labeling of
space that makes the observed molecules most probable, given what cells' expression and sizes are
like. situseg samples labelings from their posterior, as the publication does, and reports the
posterior's summary. The model and its updates are written out below (the notation is this
document's), then the sampler, the start, the outputs and the validation. Where the publication
leaves a choice open, or where situseg departs from it, the text says so.

## The model

**Space.** The region $R$ (a patch of the section, of volume $|R|$ µm³) is a lattice of voxels,
each of volume $\nu$ (its side in x and y times the section's depth divided by the number of z
layers). A labeling $\sigma$ gives each voxel $v$ a cell $\sigma(v) \in \{1, \dots, C\}$ or the
background, $\sigma(v) = 0$. Cell $c$'s volume is $V_c = \nu \, |\{v : \sigma(v) = c\}|$. A run
can sample in stages of finer voxels, each splitting the last's without changing which cell owns
a region.

**Molecules.** Transcript $i$ (gene transcripts of quality 20 or more, by default) has gene $g_i$
and lies in voxel $v(i)$; write $\sigma(i) = \sigma(v(i))$. The molecules of gene $g$ are a Poisson
process of intensity $\lambda_{0g} + \lambda_{\sigma(x) g}$ at $x$ (per µm³; the second term is
absent over the background): the background's rate $\lambda_{0g}$ everywhere and cell $c$'s rate
$\lambda_{cg}$ within it. The log-likelihood of the labeling and the rates is

$$\log p(\text{molecules} \mid \sigma, \lambda) = \sum_i \log\bigl(\lambda_{0 g_i} + \lambda_{\sigma(i) g_i}\bigr) - |R| \sum_g \lambda_{0g} - \sum_c V_c \Lambda_c, \qquad \Lambda_c = \sum_g \lambda_{cg}.$$

Each molecule in a cell came from the cell or from the background: the attribution $a_i$, with
$P(a_i = \text{cell}) = \lambda_{\sigma(i) g_i} / (\lambda_{0 g_i} + \lambda_{\sigma(i) g_i})$ (the
superposition of two Poisson processes). Given the attributions, cell $c$'s count of gene $g$,
$n_{cg}$, is Poisson with mean $V_c \lambda_{cg}$, and the background's, $n_{0g}$, Poisson with mean
$|R| \lambda_{0g}$.

**The background, of two kinds.** Molecules from no cell arise two ways: decoding errors and probes
bound off target, which fall on every codeword alike (false decodes), and RNA loose in the tissue,
which falls on each gene by its own rate (ambient RNA). So $\lambda_{0g} = f + a_g$: a false-decode
rate $f$ shared by every codeword and an ambient rate $a_g$ for each gene. The panel's
negative-control codewords (probes and codewords that match no gene, $C$ of them) are false
decodes only: their $n_c$ transcripts in the region are Poisson with mean $|R|\,C f$, evidence of
the background that no cell can take. *An extension of situseg's:* the publication has one
background rate per gene and does not model the negative controls.

**Expression.** Each cell belongs to one of $K$ components (broad expression states),
$z_c \in \{1, \dots, K\}$, with probabilities $\pi \sim \mathrm{Dirichlet}(\alpha, \dots, \alpha)$.
Given its component, a cell's rates are gamma (shape $r$, rate $\beta$):

$$\lambda_{cg} \mid z_c = k \sim \mathrm{Gamma}(r_{kg}, \beta_{kg}), \quad r_{kg} \sim \mathrm{Gamma}(a_r, b_r), \quad \beta_{kg} \sim \mathrm{Gamma}(a_\beta, \phi_g), \quad \phi_g \sim \mathrm{Gamma}(e_\phi, f_\phi), \quad f, a_g \sim \mathrm{Gamma}(a_0, b_0).$$

A cell's count, its rate integrated out, is negative binomial:
$n_{cg} \sim \mathrm{NB}(r_{kg}, p_{ckg})$, $p_{ckg} = V_c / (V_c + \beta_{kg})$. The gene-level
$\phi_g$ shares each gene's scale across components, so a component with few cells is drawn toward
the gene's typical rate.

*situseg's choice:* the publication parameterizes the negative binomial's probability
logistically and samples it with Pólya–Gamma variables. situseg places a gamma prior on the rate
$\beta_{kg}$ instead, conjugate given the cells' rates, so that every expression update is an exact
draw and the sampler's target is exactly the model's posterior, which the tests check by
enumeration.

**Volume.** A cell's volume is log-normal by component, $\log V_c \mid z_c = k \sim \mathcal{N}(m_k, s_k^2)$,
with the conjugate normal–inverse-gamma prior $s_k^2 \sim \mathrm{InvGamma}(a_s, b_s)$,
$m_k \mid s_k^2 \sim \mathcal{N}(m_0, s_k^2/\kappa_0)$.

**The nuclear prior.** A transcript with a prior cell $q_i$ (over that cell's nucleus; see the
start, below) weighs the labeling by $\rho$ if it lies in $q_i$ and $1 - \rho$ if it does not
($\rho = 0.8$, the publication's).

**The labeling's prior** is
$p(\sigma \mid z, m, s) \propto I(\sigma) \prod_c \mathrm{LogNormal}(V_c; m_{z_c}, s^2_{z_c}) \prod_{i : q_i \ne 0} \rho^{[\sigma(i) = q_i]} (1 - \rho)^{[\sigma(i) \ne q_i]}$,
where $I(\sigma)$ is 1 when every cell is connected, within its perimeter bound and non-empty, and
0 otherwise. *The one approximation:* this prior's normalizing constant depends on the volume
parameters and the components and cannot be computed (the "doubly intractable" problem). The
components and the volume parameters are updated as though the cells' volumes were independent
log-normal draws. With one volume prior for every component, held fixed (`exact`), the sampler is
exact; the tests that check exactness use that mode.

**The updates.** One sweep of the parameters, given the labeling, draws in this order:

1. **Attributions:** each molecule in a cell from the cell with probability
   $\lambda_{\sigma(i) g_i} / (\lambda_{0 g_i} + \lambda_{\sigma(i) g_i})$; the counts follow.
2. **Background:** each gene's $n_{0g}$ splits into false decodes and ambient RNA,
   $k_g \sim \mathrm{Binomial}(n_{0g}, f / (f + a_g))$; then
   $a_g \sim \mathrm{Gamma}(a_0 + n_{0g} - k_g,\; b_0 + |R|)$ and
   $f \sim \mathrm{Gamma}(a_0 + n_c + \sum_g k_g,\; b_0 + |R|(C + G))$.
3. **Components**, each cell's rates integrated out:
   $P(z_c = k) \propto \pi_k \, \mathrm{LogNormal}(V_c; m_k, s_k^2) \prod_g \mathrm{NB}(n_{cg}; r_{kg}, p_{ckg})$.
4. **Shapes**, the rates integrated out, by Chinese-restaurant-table augmentation (Zhou M, Carin
   L. *IEEE TPAMI* 2015;37:307–320): $\ell_{cg} \sim \mathrm{CRT}(n_{cg}, r_{kg})$, then
   $r_{kg} \sim \mathrm{Gamma}\bigl(a_r + \sum_c \ell_{cg},\; b_r + \sum_{c : z_c = k} \log(1 + V_c/\beta_{kg})\bigr)$.
5. **Rates:** $\lambda_{cg} \sim \mathrm{Gamma}(r_{z_c g} + n_{cg},\; \beta_{z_c g} + V_c)$.
6. **Rate parameters:** $\beta_{kg} \sim \mathrm{Gamma}\bigl(a_\beta + N_k r_{kg},\; \phi_g + \sum_{c : z_c = k} \lambda_{cg}\bigr)$.
7. **Gene scales:** $\phi_g \sim \mathrm{Gamma}(e_\phi + K a_\beta,\; f_\phi + \sum_k \beta_{kg})$.
8. **Proportions:** $\pi \sim \mathrm{Dirichlet}(\alpha + N_1, \dots, \alpha + N_K)$.
9. **Volumes:** each component's normal–inverse-gamma posterior given its cells' log volumes
   (skipped in `exact` mode).

Steps 3 and 4 integrate the cells' rates out, and step 5 draws them anew before anything
conditions on them: a partially collapsed Gibbs sampler, valid in this order (van Dyk DA, Park T.
*JASA* 2008;103:790–796). Every cell's draws come from its own random stream, named by the sweep,
the step and the cell, so cells are updated in parallel with the same result on any number of
threads; logarithms and exponentials are computed the same way on every processor, so the same
seed gives the same bytes everywhere.

**Defaults:**

| Parameter | Default | Why |
| --- | --- | --- |
| Components $K$ | 10 | Enough for a panel's broad cell types; empty components fall back to the prior |
| $\alpha$ | 1 | A uniform prior over the proportions |
| $a_r, b_r$ | 2, 0.5 | Shapes of mean 4 (sd 2.8): overdispersed, without ruling out near-Poisson genes |
| $a_\beta$ | 1 | An exponential prior on each component's rate parameter about the gene's scale |
| $e_\phi, f_\phi$ | 1, 1 µm³ | Weak: rate parameters are hundreds to thousands of µm³ |
| $a_0, b_0$ | 1, 1 µm³ | Weak (on $f$ and each $a_g$) |
| $m_0$ | $\log$ of the median initial cell volume | Centers the volume prior on the cells the sampler starts from |
| $\kappa_0, a_s, b_s$ | 1, 2, 0.25 | A prior mean of 0.25 for $s_k^2$: volumes within a factor of about 1.6 |
| $\rho$ | 0.8 | The publication's |

## The sampler

Each sweep draws the attributions and the parameters, then the labeling, then the transcripts'
places.

**The labeling's moves** (Metropolis–Hastings). A voxel $v$ is chosen. With probability
$p_b = 0.05$ the background is proposed for it (a "bubble", which lets gaps open between cells);
otherwise the label of one of its face neighbors, chosen uniformly. The proposal's probability both
ways is local, $q(a \to b) \propto p_b [b = 0] + (1 - p_b)\, n_b(v)/n(v)$, and the move is accepted
with probability $\min\{1, e^{\Delta}\, q(b \to a)/q(a \to b)\}$, $\Delta$ the change in the
log-posterior (attributions integrated out).

**The constraints,** each holding in both directions so that a move is allowed exactly when its
reverse is:
- **No cell is emptied.**
- **Connectivity.** Among the 26 voxels about $v$, those of the cell losing $v$ and those of the
  cell gaining it must each form exactly one face-connected piece touching $v$'s faces (a local
  test of simple points: Kong TY, Rosenfeld A. *Comput Vis Graph Image Process* 1989;48:357–393).
  *situseg's choice:* the test applies to both cells; a test of the losing cell alone allows moves
  whose reverse it forbids.
- **Perimeter.** In each layer a cell's perimeter, in voxel edges, is at most $f$ times the least
  any shape of its area can have, $2\lceil 2\sqrt{n}\rceil$ for $n$ voxels (Harary F, Harborth H.
  *J Combin Inform System Sci* 1976;1:1–8), $f = 1.6$. *situseg's choice:* the publication bounds
  the perimeter against a circle's; in voxel edges a digital disc's perimeter already exceeds the
  circle's by about 4/π. 1.6 is the tightest bound that still admits an ellipse four times as long
  as wide at any angle. A cell that starts beyond the bound keeps its own starting ratio.

**Parallel and deterministic.** The region is cut into square blocks (40 µm, all layers), colored
in a 2 × 2 checkerboard shifted each sweep by an offset drawn from the sweep's stream. In each of
four phases the blocks of one color are updated in parallel, each proposing as many moves as it
has voxels from its own random stream. A cell may move in a phase only if its bounding box, grown
by one voxel, meets exactly one of the phase's blocks, so blocks change disjoint cells and a phase
is the same as its blocks updated one after another.

**Repositioning.** Each transcript's true place in x and y is its observed place displaced by a
normal of σ 4 µm with probability 0.2, else of σ 0.5 µm (the publication's mixture). Each sweep
proposes a place for every transcript from that prior and accepts it with probability
$\min\{1, \Lambda(\text{new})/\Lambda(\text{old})\}$, $\Lambda$ the gene's intensity at the
place's voxel times the transcript's nuclear prior there.

**Defaults:** one stage of 1 µm voxels in one layer, 600 sweeps, the last 500 summarized. The
publication's schedule (4 µm and one layer, then 2 µm and two, then 1 µm and four) is available
as a setting; on simulated tissues the single fine stage lands closer to the truth.

## The start

situseg starts from the cells the instrument gave the transcripts (a table's cell column: Xenium's
`cell_id`, CosMx's `cell_ID` within its field of view, MERSCOPE's `cell_id`), drawn from those
transcripts alone:

- **The start raster** (0.5 µm pixels): every pixel holding transcripts the instrument gave a
  cell is that cell's (the cell most of them were given), and every other pixel the cell of the
  nearest such pixel within 2.5 µm, by an exact Euclidean distance transform (Felzenszwalb PF,
  Huttenlocher DP. *Theory of Computing* 2012;8:415–428): each cell's region about its
  transcripts, split between neighbors along the line midway between them.
- **On the first lattice,** each voxel column takes the label most of its pixels have, and each
  cell keeps its piece holding its seed.
- **Seeds:** the centroid of the cell's transcripts over its nucleus (the table's nucleus column:
  Xenium's `overlaps_nucleus`, CosMx's `CellComp`), else of all of its transcripts.
- **Prior cells:** a transcript over its cell's nucleus has that cell as its prior; in a table
  without a nucleus column, every transcript the instrument gave a cell has it as its prior.
- **Small cells:** a cell left without a column takes its seed's column; if a cell still has none
  (two seeds in one voxel), the run starts at the next finer stage, and cells lost even at the
  finest are left out and counted.

On the standard simulated tissue this start gives the same accuracy as starting from the
instrument's cell and nucleus masks (below).

## Regions, sections and patches

- **A region's cells** are those whose seed lies in it. They are refined on a lattice covering the
  region and a margin of 20 µm about it; the cells about it are held where they started, so the
  region's edge never cuts a cell. The rest of the tissue keeps the cells it started from.
- **Sections** are the tissue's separate pieces: cells grouped where their boxes come within a gap
  (100 µm unless set), each refined apart, with a model of its own.
- **Patches:** a region of more than twice 2 million modeled transcripts (`-patch-transcripts`) is
  cut into equal rectangles; each refines the cells whose seeds lie within 40 µm of its rectangle,
  holding those about them, one patch at a time. Each pixel and transcript is taken from the patch
  whose rectangle holds it. The patches depend on the data and the settings only, so a run is still
  the same bytes for the same seed. (After Blampey et al. 2024, *Nat Commun* 15:4981.)

## What is reported

The posterior's summary over the last sweeps, not one sample: each voxel's most frequent label,
and each transcript's most frequent cell with its posterior probability (the share of samples).

- **The 2D footprint:** each voxel column's dominant label across layers, with every hole filled
  (a region of the background enclosed by one cell is the cell's: a summary leaves holes where no
  transcript happened to lie, but cells have none).
- **The label image and outlines:** the footprint at the start raster's resolution; outlines traced
  from it and simplified to half a voxel (Douglas–Peucker). With several layers, each layer's
  outlines too.
- **Each transcript's cell:** a modeled transcript in a refined cell, its most frequent cell; every
  other (controls, low quality, held cells, beyond the regions), the cell holding the pixel it was
  observed in.
- **Counts:** gene transcripts of quality 20 or more (all, for a table without quality scores), by
  assigned cell.
- **Each cell's flags:** fewer than 10 gene transcripts (its boundary comes more from the prior
  than from the transcripts); uncertain boundaries (its pixels its own in fewer than 75% of the
  summarized sweeps, on average: its boundary certainty); no transcript over a nucleus.
- **Convergence:** Geweke's diagnostic of the log-likelihood over the summarized sweeps (within ±2
  it has settled), the split-half agreement (each transcript's most frequent cell over each half of
  the summarized sweeps), the voxels moved per sweep, and the background's shares. With the
  defaults the log-likelihood is usually still drifting slowly (Geweke's z about 12 to 15) while
  the summary barely changes with longer runs; more sweeps make sure.

Limitations: no cell absent from the start can be created, and boundaries between cells of like
expression, or of cells with few transcripts, are uncertain.

## Validation

**The sampler** (the package's tests, fixed seeds):
- **Each conjugate update against its analytic posterior** (Kolmogorov–Smirnov): the cells' and the
  components' rates, the gene scales, the proportions and the volume prior; the background's data
  augmentation against the exact mixture posterior, enumerated.
- **The whole sweep** by Geweke's joint distribution test (Geweke J. *J Am Stat Assoc*
  2004;99:799–804): chains started from the prior with data drawn from the model stay distributed
  as the prior.
- **Detailed balance, exactly:** on lattices small enough to list every labeling (3 × 3, and
  2 × 2 × 2), every move from every reachable state balances to $10^{-9}$, and no move is allowed
  one way only; 20,000 chains started from exact draws end distributed as the target.
- **Moves:** each move's computed change equals the target recomputed afresh.
- **Determinism:** a run gives the same digest on 1, 4 and 16 threads, pinned so that every
  processor is checked against the same bytes.

**Against the truth** (`situseg simulate`, `situseg score`; the standard simulated tissue: 500 ×
500 µm, 3,030 true cells, 2,925 seen by the simulated instrument, 394,507 gene transcripts of
quality 20 or more). A transcript is in its true cell when the segmentation's cell holding it is
matched to its true cell, cells matched one to one where each is the other's best overlap; the
Jaccard index is of a matched pair's transcripts.

| Segmentation | In their true cell | Median Jaccard | Time | Memory |
| --- | --- | --- | --- | --- |
| The simulated instrument (a 5 µm nucleus expansion; the start) | 82.5% | 0.701 | | |
| situseg, Balanced, from the transcripts table | **86.8%** | **0.770** | 55 s | 127 MB |
| The same model started from the instrument's cell and nucleus masks (SituCraft 0.9.0) | 86.7% | 0.770 | | |
| Proseg 3.2.0, `--xenium --voxel-layers 1` (three runs) | 78.5–78.7% | 0.640–0.644 | 56 s | 127 MB |
| Proseg 3.2.0, `--xenium` (its defaults, four layers) | 80.5% | 0.683 | | |

Times and peak resident memory are on a 10-core laptop under sustained load, situseg's and
Proseg's runs back to back.

The presets on the same tissue (`TestPresetsOnSimulatedTissue`, `SITUSEG_PRESETS=standard`;
"outlines" is the mean ratio of each cell's outline, in pixel edges, to the least its area can
have):

| Preset | What it changes from Balanced | In their true cell | Median Jaccard | Outlines | Peak heap |
| --- | --- | --- | --- | --- | --- |
| Balanced (the defaults) | | 86.8% | 0.768 | 1.31 | 111 MB |
| Conservative boundary | perimeters within 1.3 times the least; gaps proposed at 0.1; nuclear prior 0.9 | 86.4% | 0.767 | 1.21 | 111 MB |
| Flexible morphology | perimeters within 2.5 times the least; nuclear prior 0.7 | 85.3% | 0.746 | 1.61 | 112 MB |
| Low memory | 2 µm voxels; the collector's target 20% | 85.5% | 0.748 | 1.35 | 82 MB |
| ROI pilot (a 200 µm square in the middle) | a region, run with two seeds | 86.7% in the region (Balanced on the whole tissue: 87.3%) | 0.770 (0.772) | 1.36 | 73 MB |

Every preset puts more transcripts in their true cells than the instrument. The pilot's two seeds
put 93.3% of the transcripts in refined cells in the same cell.

**On real data** (10x Genomics Xenium human lung, 2 FOVs, example data under CC BY 4.0: 11,898
cells, 746,213 modeled transcripts; no truth): situseg took 4 minutes and 307 MB on 10 cores, as
the mask-started model does (4 minutes, 279 MB, back to back). The two put 86.3% of transcripts in
the same cell, as many as two seeds of the mask-started model do (86.1%); each agrees with the
instrument on about 86%, and Proseg 3.2.0 with any of the three on about 60%.

The simulated tissues come from a simulator that shares no code with the segmentation; the
simulated instrument segments with a distance of its own (a 5-7-11 chamfer transform). They are
drawn from a model of the same family as the method's, which favors it: the comparison with
Proseg on real data, where no truth is known, shows how differently the two methods draw cells,
not which is right.
