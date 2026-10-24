# Changelog

## 0.1.0 (unreleased)

The first release as a program and package of its own; the method was SituCraft's
transcript-informed segmentation.

- Segmentation after Jones et al. (Nat Methods 2025), from a table of transcripts alone: the
  cells to start from are drawn from the transcripts the instrument gave each cell.
- Tables of Xenium (Parquet or CSV), CosMx and MERSCOPE recognized by their columns, and any
  other table by named columns; Parquet, CSV and TSV, gzipped or not.
- Outputs: cells and transcripts as CSV or Parquet, the counts as a 10x Matrix Market folder, the
  outlines as GeoJSON (or Xenium-style Parquet), a 32-bit label TIFF, diagnostics and the run's
  record.
- Presets (balanced, conservative, flexible, low-memory, roi-pilot), regions, sections, patches
  for large tissues, a second seed's agreement, a preflight estimate of memory and time.
- `situseg simulate` and `situseg score`: simulated tissues with their truth, and segmentations
  scored against it.
