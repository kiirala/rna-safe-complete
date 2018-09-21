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
	for i := 0; i < numBases; i++ {
		p.PairSafety[i] = make([]*big.Int, numBases)
		for j := i + 1; j < numBases; j++ {
			p.PairSafety[i][j] = new(big.Int)
		}
	}
	p.SingleSafety = make([]*big.Int, numBases)
	for i := 0; i < numBases; i++ {
		p.SingleSafety[i] = new(big.Int)
	}
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

	splits := p.findSplits(i, j)

	for _, s := range splits {
		branch := folding
		if len(splits) > 1 {
			branch = &types.FoldTree{}
			folding.Branches = append(folding.Branches, branch)
		}

		partNumSols := new(big.Int).Mul(p.Sol[s.Pre.I][s.Pre.J], solScale)
		if s.Suf != nil {
			partNumSols.Mul(partNumSols, p.Sol[s.Suf.I][s.Suf.J])
		}

		if s.Pair != nil {
			p.PairSafety[s.Pair.I][s.Pair.J].Add(p.PairSafety[s.Pair.I][s.Pair.J], partNumSols)
			branch.Pairs = append(branch.Pairs, types.Pair{I: s.Pair.I, J: s.Pair.J})
		}

		if s.Suf != nil {
			branch.JoinPrefix = &types.FoldTree{}
			p.recursiveBacktrack(s.Pre.I, s.Pre.J, branch.JoinPrefix, partNumSols)
			branch.JoinSuffix = &types.FoldTree{}
			p.recursiveBacktrack(s.Suf.I, s.Suf.J, branch.JoinSuffix, partNumSols)
		} else {
			p.recursiveBacktrack(s.Pre.I, s.Pre.J, branch, partNumSols)
		}

	}

	folding.SimplifyImmediate()
}
