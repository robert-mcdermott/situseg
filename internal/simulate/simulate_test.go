package simulate

import (
	"context"
	"reflect"
	"testing"
)

// The same parameters give the same tissue; the instrument sees most cells, and most molecules
// are cells'.
func TestSimulateIsDeterministic(t *testing.T) {
	p := Presets["small"]
	p.WidthUm, p.HeightUm = 80, 60
	a, err := Simulate(context.Background(), p)
	if err != nil {
		t.Fatal(err)
	}
	b, err := Simulate(context.Background(), p)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(a, b) {
		t.Fatal("two simulations of the same parameters differ")
	}
	if len(a.TrueIDs) < 30 || len(a.CellIDs) < len(a.TrueIDs)*8/10 || len(a.CellIDs) > len(a.TrueIDs) {
		t.Errorf("%d true cells, %d seen by the instrument", len(a.TrueIDs), len(a.CellIDs))
	}
	fromCells := 0
	for _, o := range a.Origin {
		if o == FromCell {
			fromCells++
		}
	}
	if fromCells < len(a.Origin)*9/10 {
		t.Errorf("%d of %d molecules from cells", fromCells, len(a.Origin))
	}
	p.Seed++
	c, err := Simulate(context.Background(), p)
	if err != nil {
		t.Fatal(err)
	}
	if reflect.DeepEqual(a.X, c.X) {
		t.Error("another seed gave the same tissue")
	}
}

// Scores: a segmentation that is the truth scores 1; one that puts every transcript in one cell
// scores far less.
func TestCompare(t *testing.T) {
	p := Presets["small"]
	p.WidthUm, p.HeightUm = 60, 40
	tis, err := Simulate(context.Background(), p)
	if err != nil {
		t.Fatal(err)
	}
	counted := func(i int) bool { return !tis.Control[tis.Feature[i]] && tis.QV[i] >= 20 }
	s := Compare(tis.TrueCell, tis.TrueCell, len(tis.TrueIDs), len(tis.TrueIDs), counted)
	if s.InTrueCell != 1 || s.MedianJaccard != 1 || s.MatchedCells == 0 {
		t.Errorf("the truth against itself: %+v", s)
	}
	one := make([]int32, len(tis.TrueCell))
	if s := Compare(tis.TrueCell, one, len(tis.TrueIDs), 1, counted); s.InTrueCell > 0.2 || s.MatchedCells != 1 {
		t.Errorf("one cell for everything: %+v", s)
	}
}
