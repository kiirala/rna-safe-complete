package safecomplete

import "log"
import "math/big"

import "keltainen.duckdns.org/rnafolding/types"

func (p *Predictor) BacktrackAll() *types.FoldTree {
	p.BacktrackSafety()
	return p.BacktrackFolding()
}

func (p *Predictor) BacktrackFolding() *types.FoldTree {
	numBases := len(p.Seq.Bases)
	folding := &types.FoldTree{}
	p.recursiveFolding(0, numBases-1, folding)
	return folding
}

func (p *Predictor) BacktrackSafety() {
	if p.Sol == nil {
		log.Fatal("Must run safecomplete.CountSolutions before safecomplete.SafetyBacktrack")
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

	p.recursiveSafety(0, numBases-1, p.Sol[0][numBases-1])
}

func (p *Predictor) recursiveFolding(i, j int, folding *types.FoldTree) {
	if i >= j {
		if i == j {
			folding.Free = append(folding.Free, i)
		}
		return
	}

	splits := p.findSplits(i, j)

	for _, s := range splits {
		branch := folding
		if len(splits) > 1 {
			branch = &types.FoldTree{}
			folding.Branches = append(folding.Branches, branch)
		}

		if s.Pair != nil {
			branch.Pairs = append(branch.Pairs, types.Pair{I: s.Pair.I, J: s.Pair.J})
		}

		if s.Suf != nil {
			branch.JoinPrefix = &types.FoldTree{}
			p.recursiveFolding(s.Pre.I, s.Pre.J, branch.JoinPrefix)
			branch.JoinSuffix = &types.FoldTree{}
			p.recursiveFolding(s.Suf.I, s.Suf.J, branch.JoinSuffix)
		} else {
			p.recursiveFolding(s.Pre.I, s.Pre.J, branch)
		}

	}

	folding.SimplifyImmediate()
}

func (p *Predictor) recursiveSafety(i, j int, numSols *big.Int) {
	if i >= j {
		if i == j {
			p.SingleSafety[i].Add(p.SingleSafety[i], numSols)
		}
		return
	}

	solScale := new(big.Int)
	solMod := new(big.Int)
	solScale.DivMod(numSols, p.Sol[i][j], solMod)
	if solMod.Cmp(big.NewInt(0)) != 0 {
		log.Printf("recursiveFolding(%d, %d, ...): numSols=%v is not divisible by p.Sol=%v",
			i, j, numSols, p.Sol[i][j])
	}

	splits := p.findSplits(i, j)

	for _, s := range splits {
		partNumSols := new(big.Int).Mul(p.Sol[s.Pre.I][s.Pre.J], solScale)
		if s.Suf != nil {
			partNumSols.Mul(partNumSols, p.Sol[s.Suf.I][s.Suf.J])
		}

		if s.Pair != nil {
			p.PairSafety[s.Pair.I][s.Pair.J].Add(p.PairSafety[s.Pair.I][s.Pair.J], partNumSols)
		}

		if s.Suf != nil {
			p.recursiveSafety(s.Pre.I, s.Pre.J, partNumSols)
			p.recursiveSafety(s.Suf.I, s.Suf.J, partNumSols)
		} else {
			p.recursiveSafety(s.Pre.I, s.Pre.J, partNumSols)
		}
	}
}
