package simulate

// Scoring a segmentation against the truth, transcript by transcript: the cells of the two are
// matched one to one where each is the other's best overlap (the cell sharing most of its counted
// transcripts), and a counted transcript the truth puts in a cell is in its true cell when the
// segmentation puts it in the cell matched to it.

import "sort"

// Score is how a segmentation does against the truth.
type Score struct {
	// InTrueCell is the share of the counted transcripts the truth puts in a cell that the
	// segmentation puts in the cell matched to it.
	InTrueCell float64 `json:"inTrueCell"`
	// MatchedCells of the TrueCells are matched; MedianJaccard is the median, over them, of the
	// two cells' transcripts' Jaccard index (shared over either's).
	MatchedCells  int     `json:"matchedCells"`
	TrueCells     int     `json:"trueCells"`
	MedianJaccard float64 `json:"medianJaccard"`
	// Assigned and TrueAssigned are the shares of the counted transcripts in a cell.
	Assigned     float64 `json:"assigned"`
	TrueAssigned float64 `json:"trueAssigned"`
}

// Compare scores pred against truth (each transcript's cell, −1 none; trueCells and predCells
// the numbers of cells), over the transcripts counted says to count.
func Compare(truth, pred []int32, trueCells, predCells int, counted func(i int) bool) Score {
	type pair struct{ t, p int32 }
	overlap := map[pair]int{}
	tSize, pSize := make([]int, trueCells), make([]int, predCells)
	total, truthIn, predIn := 0, 0, 0
	for i := range truth {
		if !counted(i) {
			continue
		}
		total++
		tc, pc := truth[i], pred[i]
		if tc >= 0 {
			tSize[tc]++
			truthIn++
		}
		if pc >= 0 {
			pSize[pc]++
			predIn++
		}
		if tc >= 0 && pc >= 0 {
			overlap[pair{tc, pc}]++
		}
	}
	bestP, bestT := make([]int32, trueCells), make([]int32, predCells)
	bestPN, bestTN := make([]int, trueCells), make([]int, predCells)
	for i := range bestP {
		bestP[i] = -1
	}
	for i := range bestT {
		bestT[i] = -1
	}
	for k, n := range overlap {
		if n > bestPN[k.t] || n == bestPN[k.t] && k.p < bestP[k.t] {
			bestP[k.t], bestPN[k.t] = k.p, n
		}
		if n > bestTN[k.p] || n == bestTN[k.p] && k.t < bestT[k.p] {
			bestT[k.p], bestTN[k.p] = k.t, n
		}
	}
	s := Score{TrueCells: trueCells}
	var jac []float64
	correct := 0
	for t := range tSize {
		p := bestP[t]
		if p < 0 || bestT[p] != int32(t) {
			continue
		}
		s.MatchedCells++
		n := overlap[pair{int32(t), p}]
		correct += n
		jac = append(jac, float64(n)/float64(tSize[t]+pSize[p]-n))
	}
	sort.Float64s(jac)
	if len(jac) > 0 {
		s.MedianJaccard = jac[len(jac)/2]
	}
	if truthIn > 0 {
		s.InTrueCell = float64(correct) / float64(truthIn)
	}
	if total > 0 {
		s.Assigned, s.TrueAssigned = float64(predIn)/float64(total), float64(truthIn)/float64(total)
	}
	return s
}
