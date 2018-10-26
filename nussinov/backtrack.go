package nussinov /* import "keltainen.duckdns.org/rnafolding/nussinov" */

import "keltainen.duckdns.org/rnafolding/folding"

func (p *Predictor) Backtrack() (int, folding.FoldingPairs) {
	paired := make(folding.FoldingPairs, len(p.V))
	for i := 0; i < len(p.V); i++ {
		paired[i] = -1
	}
	p.recursiveBacktrack(0, len(p.V)-1, paired)
	return p.V[0][len(p.V)-1], paired
}

func (p *Predictor) JoinedBacktrack(i, j int) (int, folding.FoldingPairs) {
	paired := make(folding.FoldingPairs, len(p.V))
	pairs := 0
	for i := 0; i < len(p.V); i++ {
		paired[i] = -1
	}
	if p.Seq.CanPair(i, j, p.MinHairpin) {
		paired[i] = j
		paired[j] = i
		pairs++
	}
	p.recursiveBacktrack(i+1, j-1, paired)
	pairs += p.V[i+1][j-1]
	p.complementaryBacktrack(i, j, paired)
	pairs += p.W[i][j]
	return pairs, paired
}

func (p *Predictor) recursiveBacktrack(i, j int, paired folding.FoldingPairs) {
	if i >= j {
		return
	} else if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
		paired[i] = j
		paired[j] = i
		p.recursiveBacktrack(i+1, j-1, paired)
	} else {
		for k := i; k < j; k++ {
			if p.V[i][j] == p.V[i][k]+p.V[k+1][j] {
				p.recursiveBacktrack(i, k, paired)
				p.recursiveBacktrack(k+1, j, paired)
				break
			}
		}
	}
}

func (p *Predictor) complementaryBacktrack(i, j int, paired folding.FoldingPairs) {
	if j-i >= len(p.V)-2 {
		return
	}
	if i == 0 {
		p.recursiveBacktrack(j+1, len(p.V)-1, paired)
	} else if j == len(p.V)-1 {
		p.recursiveBacktrack(0, i-1, paired)
	} else if p.Seq.CanPair(i-1, j+1, p.MinHairpin) && p.W[i][j] == p.W[i-1][j+1]+1 {
		paired[i-1] = j + 1
		paired[j+1] = i - 1
		p.complementaryBacktrack(i-1, j+1, paired)
	} else {
		for k := j + 1; k < i+len(p.V)-1; k++ {
			if k < len(p.V) && p.W[i][j] == p.V[j+1][k]+p.W[i][k] {
				p.recursiveBacktrack(j+1, k, paired)
				p.complementaryBacktrack(i, k, paired)
				break
			} else if k >= len(p.V) && p.W[i][j] == p.W[k+1-len(p.V)][j]+p.V[k+1-len(p.V)][i-1] {
				p.complementaryBacktrack(k+1-len(p.V), j, paired)
				p.recursiveBacktrack(k+1-len(p.V), i-1, paired)
				break
			}
		}
	}
}
