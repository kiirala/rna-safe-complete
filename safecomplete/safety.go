package safecomplete /* import "keltainen.duckdns.org/rnafolding/safecomplete" */

import "log"
import "math/big"

import "keltainen.duckdns.org/rnafolding/folding"
import "keltainen.duckdns.org/rnafolding/types"

func TrivialSafety(foldings folding.FoldingSet) []bool {
	out := make([]bool, len(foldings[0]))
	for i := 0; i < len(out); i++ {
		out[i] = true
	}
	for _, f := range foldings {
		for i := 0; i < len(f); i++ {
			if f[i] != foldings[0][i] && (f[i] >= 0 || foldings[0][i] >= 0) {
				out[i] = false
			}
		}
	}
	return out
}

func (p *Predictor) IteratedTrivialSafety(folds *types.FoldTree) []bool {
	out := make([]bool, len(p.Seq.Bases))
	for i := 0; i < len(out); i++ {
		out[i] = true
	}
	var ref folding.FoldingPairs
	for f := range folds.PairArrayIterator(p.Seq) {
		if ref == nil {
			ref = f
		} else {
			for i := 0; i < len(out); i++ {
				if f[i] != ref[i] && (f[i] >= 0 || ref[i] >= 0) {
					out[i] = false
				}
			}
		}
	}
	return out
}

func (p *Predictor) SafetyFromBacktrack() []bool {
	out := make([]bool, len(p.Seq.Bases))
	max := p.Sol[0][len(p.Seq.Bases)-1]
	for i := 0; i < len(out); i++ {
		sum := new(big.Int).Set(p.SingleSafety[i])
		numPairings := 0
		for j := 0; j < i; j++ {
			if p.PairSafety[j][i].Cmp(new(big.Int)) > 0 {
				numPairings++
				sum.Add(sum, p.PairSafety[j][i])
			}
		}
		for j := i + 1; j < len(out); j++ {
			if p.PairSafety[i][j].Cmp(new(big.Int)) > 0 {
				numPairings++
				sum.Add(sum, p.PairSafety[i][j])
			}
		}
		if sum.Cmp(max) != 0 {
			log.Printf("SafetyFromBacktrack: base %d is in %d solutions, expected %d", i, sum, max)
		}
		if (p.SingleSafety[i].Cmp(new(big.Int)) == 0 && numPairings == 1) ||
			(p.SingleSafety[i].Cmp(new(big.Int)) > 0 && numPairings == 0) {
			out[i] = true
		} else {
			out[i] = false
		}
	}
	return out
}
