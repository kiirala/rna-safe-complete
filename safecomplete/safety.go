package safecomplete

import "log"

import "keltainen.duckdns.org/rnafolding/types"

func TrivialSafety(foldings types.FoldingSet) []bool {
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

func (p *Predictor) SafetyFromBacktrack() []bool {
	out := make([]bool, len(p.Seq.Bases))
	max := p.Sol[0][len(p.Seq.Bases)-1]
	for i := 0; i < len(out); i++ {
		sum := p.SingleSafety[i]
		numPairings := 0
		for j := 0; j < i; j++ {
			if p.PairSafety[j][i] > 0 {
				numPairings++
				sum += p.PairSafety[j][i]
			}
		}
		for j := i; j < len(out); j++ {
			if p.PairSafety[i][j] > 0 {
				numPairings++
				sum += p.PairSafety[i][j]
			}
		}
		if sum != max {
			log.Printf("SafetyFromBacktrack: base %d is in %d solutions, expected %d", i, sum, max)
		}
		if (p.SingleSafety[i] == 0 && numPairings == 1) || (p.SingleSafety[i] > 0 && numPairings == 0) {
			out[i] = true
		} else {
			out[i] = false
		}
	}
	return out
}
