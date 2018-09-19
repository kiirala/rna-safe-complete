package safecomplete

import "log"
import "math/big"

import "keltainen.duckdns.org/rnafolding/types"

func (p *Predictor) BacktrackAll() *types.FoldTree {
	if p.Sol == nil {
		log.Fatal("Must run safecomplete.CountSolutions before safecomplete.BacktrackAll")
	}
	numBases := len(p.Seq.Bases)
	p.PairSafety = make([][]*big.Int, numBases)
	for i := 0; i < len(p.Seq.Bases); i++ {
		p.PairSafety[i] = make([]*big.Int, numBases)
	}
	p.SingleSafety = make([]*big.Int, numBases)
	folding := &types.FoldTree{}
	p.recursiveBacktrack(0, numBases-1, folding, p.Sol[0][numBases-1])
	return folding
}

func (p *Predictor) recursiveBacktrack(i, j int, folding *types.FoldTree, numSols *big.Int) {
	if i >= j {
		if i == j {
			p.SingleSafety[i].Add(p.SingleSafety[i], numSols)
			folding.Free = append(folding.Free, i)
		}
		return
	}

	solScale := new(big.Int)
	solMod := new(big.Int)
	solScale.DivMod(numSols, p.Sol[i][j], solMod)
	if solMod.Cmp(big.NewInt(0)) != 0 {
		log.Printf("recursiveBacktrack(%d, %d, ...): numSols=%v is not divisible by p.Sol=%v",
			i, j, numSols, p.Sol[i][j])
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
	if p.V[i][j] == p.V[i][j-1] {
		numFound++
	}

	if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
		branch := folding
		if numFound > 1 {
			branch = &types.FoldTree{}
			folding.Branches = append(folding.Branches, branch)
		}
		partNumSols := new(big.Int).Mul(p.Sol[i+1][j-1], solScale)
		p.PairSafety[i][j].Add(p.PairSafety[i][j], partNumSols)
		branch.Pairs = append(branch.Pairs, types.Pair{I: i, J: j})
		p.recursiveBacktrack(i+1, j-1, branch, partNumSols)
	}

	for k := i; k < j; k++ {
		if p.V[i][j] == p.V[i][k]+p.W[k+1][j] && p.W[k+1][j] > 0 {
			partNumSols := new(big.Int).Mul(p.Sol[i][k], p.Sol[k+2][j-1])
			partNumSols.Mul(partNumSols, solScale)
			branch := folding
			if numFound > 1 {
				branch = &types.FoldTree{}
				folding.Branches = append(folding.Branches, branch)
			}
			branch.JoinPrefix = &types.FoldTree{}
			p.recursiveBacktrack(i, k, branch.JoinPrefix, partNumSols)

			branch.JoinSuffix = &types.FoldTree{}
			p.PairSafety[k+1][j].Add(p.PairSafety[k+1][j], partNumSols)
			branch.JoinSuffix.Pairs = append(branch.JoinSuffix.Pairs, types.Pair{I: k + 1, J: j})
			p.recursiveBacktrack(k+2, j-1, branch.JoinSuffix, partNumSols)
		}
	}

	if p.V[i][j] == p.V[i][j-1] {
		partNumSols := new(big.Int).Mul(p.Sol[i][j-1], solScale)
		branch := folding
		if numFound > 1 {
			branch = &types.FoldTree{}
			folding.Branches = append(folding.Branches, branch)
		}
		branch.JoinPrefix = &types.FoldTree{}
		p.recursiveBacktrack(i, j-1, branch.JoinPrefix, partNumSols)
		branch.JoinSuffix = &types.FoldTree{}
		p.recursiveBacktrack(j, j, branch.JoinSuffix, partNumSols)
	}

	folding.SimplifyImmediate()
}
