package safecomplete

import "keltainen.duckdns.org/rnafolding/types"

func (p *Predictor) BacktrackAll() *types.Folding {
	folding := &types.Folding{}
	p.recursiveBacktrack(0, len(p.V)-1, folding)
	return folding
}

func (p *Predictor) recursiveBacktrack(i, j int, folding *types.Folding) {
	if i >= j {
		if i == j {
			folding.Free = append(folding.Free, i)
		}
		return
	}
	numFound := 0
	if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
		numFound++
	}
	for k := i; k < j; k++ {
		if p.V[i][j] == p.V[i][k]+p.W[k+1][j] && p.W[k+1][j] > 0 {
			numFound++
		}
	}
	if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
		branch := folding
		if numFound > 1 {
			branch = &types.Folding{}
			folding.Branches = append(folding.Branches, branch)
		}
		branch.Pairs = append(branch.Pairs, types.Pair{I: i, J: j})
		p.recursiveBacktrack(i+1, j-1, branch)
	}
	for k := i; k < j; k++ {
		if p.V[i][j] == p.V[i][k]+p.W[k+1][j] && p.W[k+1][j] > 0 {
			branch := folding
			if numFound > 1 {
				branch = &types.Folding{}
				folding.Branches = append(folding.Branches, branch)
			}
			branch.JoinPrefix = &types.Folding{}
			p.recursiveBacktrack(i, k, branch.JoinPrefix)

			branch.JoinSuffix = &types.Folding{}
			branch.JoinSuffix.Pairs = append(branch.JoinSuffix.Pairs, types.Pair{I: k + 1, J: j})
			p.recursiveBacktrack(k+2, j-1, branch.JoinSuffix)
		}
	}
	if numFound == 0 {
		for k := i; k < j; k++ {
			if p.V[i][j] == p.V[i][k]+p.V[k+1][j] {
				folding.JoinPrefix = &types.Folding{}
				p.recursiveBacktrack(i, k, folding.JoinPrefix)
				folding.JoinSuffix = &types.Folding{}
				p.recursiveBacktrack(k+1, j, folding.JoinSuffix)
				break
			}
		}
	}
}
