# situseg

Transcript-informed cell segmentation for imaging-based spatial transcriptomics, as a command and
a Go library. situseg refines the cells an instrument drew (10x Genomics Xenium, NanoString CosMx,
Vizgen MERSCOPE, or any table of molecules with places, genes and cells) by sampling their
boundaries from the posterior of a model of each cell's expression, after

> Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR, Newell EW. Cell simulation as cell
> segmentation. *Nat Methods* 2025;22(6):1331–1342. doi:10.1038/s41592-025-02697-0

It needs only the instrument's table of transcripts: no images, no masks, no GPU. It is one
self-contained program in pure Go, with no dependencies beyond the Go standard library; the same
seed and settings give the same segmentation, byte for byte, on any number of threads and any
processor.

What it writes is read directly by other tools: per-cell and per-transcript tables (CSV or
Parquet), the cell-by-gene counts in the 10x Matrix Market layout (`scanpy.read_10x_mtx`,
Seurat's `Read10X`), the cells' outlines as GeoJSON polygons (or Xenium-style
`cell_boundaries.parquet`), a 32-bit label TIFF, and diagnostics of the run's convergence.

## Install

```bash
go install github.com/robert-mcdermott/situseg/cmd/situseg@latest
```

or as a library:

```bash
go get github.com/robert-mcdermott/situseg
```

## Use

```bash
situseg transcripts.parquet                 # a Xenium bundle's transcripts: columns recognized
situseg -o out tx_file.csv                  # CosMx (places in pixels of 0.12028 µm, cells per field of view)
situseg detected_transcripts.csv            # MERSCOPE
situseg -x X -y Y -gene target -cell cell \ # any other table: name its columns
        -nucleus in_nucleus molecules.csv.gz
situseg -region 1000,1000,1500,1400 transcripts.parquet   # refine only a region (µm)
situseg -preset conservative -seed-check transcripts.parquet
situseg presets                             # what each preset is for
```

It writes, to `situseg-output` unless `-o` says otherwise:

| File | What it holds |
| --- | --- |
| `cells.csv.gz` | each cell: id (the instrument's), centroid and area (µm), transcripts, genes, boundary certainty, flags |
| `transcripts.csv.gz` | each row of the input: its cell (empty: none) and that cell's posterior probability (empty: assigned by place), with its place, gene and the instrument's cell |
| `cell_feature_matrix/` | the counts of gene transcripts (Q20 or more) as `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz` |
| `cell_boundaries.geojson` | each cell's outline (µm) |
| `diagnostics.json` | convergence, the background's shares, the cells flagged |
| `run.json` | the command, every setting, the input's columns, the preflight and the time each step took |

`-parquet` writes the cell and transcript tables and the boundaries as Parquet instead, `-labels`
adds `cell_labels.tif`. With several voxel layers (`-layers`), `cell_boundaries_layers.geojson`
holds each layer's outlines.

As a library:

```go
table, err := situseg.ReadTable("transcripts.parquet", situseg.XeniumColumns())
in, idx, err := situseg.InputsFromTable(ctx, table, situseg.DefaultStartOptions())
seg, err := situseg.Segment(ctx, table, in, idx, situseg.DefaultSpec(), situseg.Job{})
// seg.Cells, seg.Assignment (each row's cell), seg.Probability, seg.Counts, seg.Outlines, seg.Labels
err = seg.WriteFiles("out", table, situseg.WriteOptions{})
```

## How well it works

On a simulated tissue with known truth (`situseg simulate -preset standard`, then
`situseg score`; 3,030 cells, 394,507 gene transcripts of quality 20 or more):

| Segmentation | Transcripts in their true cell | Median Jaccard of matched cells | Time | Memory |
| --- | --- | --- | --- | --- |
| The simulated instrument (its start) | 82.5% | 0.701 | | |
| **situseg** (Balanced, the defaults) | **86.8%** | **0.770** | 55 s | 127 MB |
| Proseg 3.2.0, one voxel layer (`--xenium --voxel-layers 1`) | 78.5% | 0.642 | 56 s | 127 MB |
| Proseg 3.2.0, its defaults (`--xenium`) | 80.5% | 0.683 | | |

(Times and peak resident memory on a 10-core laptop under sustained load, the two runs back to
back.) The simulator shares no code with the segmentation, but it draws tissues from a model of
the same family as the method's, which favors it. On the Xenium human lung 2-FOV data set (10x
Genomics example data, CC BY 4.0; 11,898 cells, 746,213 modeled transcripts) situseg took 4
minutes and 307 MB, as SituCraft's refinement does (4 minutes, 279 MB, back to back), and the two
put 86.3% of transcripts in the same cell, as much as two seeds of either agree.
[docs/method.md](docs/method.md) writes out the model, the sampler, the start, the outputs and
the validation.

Reproduce the benchmark:

```bash
situseg simulate -preset standard standard.parquet
situseg -o standard-seg standard.parquet
situseg score standard.parquet standard-seg
```

## Notes

- situseg refines the cells the instrument drew; it cannot create a cell absent from the start.
- Boundaries between cells of like expression, and of cells with few transcripts, are uncertain:
  `cells.csv.gz` flags them, and `-seed-check` measures how much of the result is chance.
- The whole table is held in memory (about 35 bytes a row) together with what the sampler holds
  (about 130 bytes a modeled transcript); a region of more than 4 million modeled transcripts is
  refined in patches. The preflight estimates the memory before a run starts and refuses a run
  the computer cannot hold (`-force` runs it anyway).
- The CosMx and MERSCOPE column presets follow those instruments' documented file layouts; they
  have been tested on files of that layout, not yet on a real run's.

## Provenance

situseg is the transcript-informed segmentation of
[SituCraft](https://github.com/robert-mcdermott/situcraft), a workbench for imaging-based spatial
transcriptomics, as a program and package of its own.

## License

Apache License 2.0; see [LICENSE](LICENSE) and [NOTICE](NOTICE).
