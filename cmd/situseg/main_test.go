package main

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"

	"github.com/robert-mcdermott/situseg"
)

// The command simulates a tissue, segments it (from Parquet and from CSV, recognizing Xenium's
// columns), records the run and scores it.
func TestCommand(t *testing.T) {
	dir := t.TempDir()
	pq, csv := filepath.Join(dir, "sim.parquet"), filepath.Join(dir, "sim.csv.gz")
	for _, out := range []string{pq, csv} {
		if err := simulateCommand([]string{"-preset", "small", "-width-um", "60", "-height-um", "50", out}); err != nil {
			t.Fatal(err)
		}
	}
	segPq, segCSV := filepath.Join(dir, "seg-pq"), filepath.Join(dir, "seg-csv")
	if err := run([]string{"-quiet", "-sweeps", "40", "-samples", "30", "-o", segPq, pq}); err != nil {
		t.Fatal(err)
	}
	if err := run([]string{csv, "-quiet", "-sweeps", "40", "-samples", "30", "-parquet", "-labels", "-o", segCSV}); err != nil {
		t.Fatal(err)
	}
	for _, f := range []string{situseg.FileCells, situseg.FileTranscripts, situseg.FileBoundaries, situseg.FileDiagnostics, "run.json",
		filepath.Join(situseg.DirMatrix, "matrix.mtx.gz")} {
		if _, err := os.Stat(filepath.Join(segPq, f)); err != nil {
			t.Error(err)
		}
	}
	for _, f := range []string{situseg.FileCellsPq, situseg.FileTranscriptsPq, situseg.FileBoundariesPq, situseg.FileLabels} {
		if _, err := os.Stat(filepath.Join(segCSV, f)); err != nil {
			t.Error(err)
		}
	}
	var record struct {
		Input struct {
			Instrument string `json:"instrument"`
			Rows       int    `json:"rows"`
		} `json:"input"`
		Spec situseg.Spec `json:"spec"`
	}
	raw, err := os.ReadFile(filepath.Join(segPq, "run.json"))
	if err != nil {
		t.Fatal(err)
	}
	if err := json.Unmarshal(raw, &record); err != nil {
		t.Fatal(err)
	}
	if record.Input.Instrument != "xenium" || record.Input.Rows == 0 || record.Spec.Preset != situseg.PresetExpert {
		t.Errorf("recorded %+v", record)
	}
	if err := scoreCommand([]string{pq, segPq, segCSV}); err != nil {
		t.Fatal(err)
	}
	printPresets()
}
