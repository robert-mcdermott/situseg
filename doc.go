// Package situseg segments cells from the transcripts of imaging-based spatial transcriptomics
// (10x Genomics Xenium, NanoString CosMx, Vizgen MERSCOPE, or any table of molecules with places,
// genes and the instrument's cells), after Jones DC, Elz AE, Hadadianpour A, Ryu H, Glass DR,
// Newell EW. Cell simulation as cell segmentation. Nat Methods 2025;22(6):1331–1342.
// doi:10.1038/s41592-025-02697-0: cells are regions of space whose molecules are Poisson processes
// with rates of their own, over a background, and a segmentation is sampled from its posterior.
// What is reported is the posterior's summary over the last sweeps: each voxel's most frequent
// cell, and each transcript's with its probability.
//
// A run reads a table (ReadTable), makes the inputs from it (InputsFromTable: the cells to start
// from are the instrument's, drawn from the transcripts it gave each), and segments (Segment),
// which samples each region (Run, RunArea) and stitches the results into one Segmentation
// (Stitch); WriteFiles writes it in formats other tools read. The same seed and settings give the
// same segmentation, byte for byte, on any number of threads and any processor.
//
// docs/method.md writes the model, the sampler and their validation out.
package situseg
