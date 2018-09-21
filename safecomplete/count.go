package safecomplete

import "log"
import "math/big"

func (p *Predictor) CountSolutions() {
	numBases := len(p.V)
	p.Sol = make([][]*big.Int, numBases)
	for i := 0; i < len(p.V); i++ {
		p.Sol[i] = make([]*big.Int, numBases)
		p.Sol[i][i] = big.NewInt(1)
		if i > 0 {
			p.Sol[i][i-1] = big.NewInt(1)
		}
	}

	for l := 2; l <= numBases; l++ {
		for i := 0; i <= numBases-l; i++ {
			j := i + l - 1
			numFound := big.NewInt(0)
			if p.V[i][j] == p.V[i][j-1] {
				numFound.Add(numFound, p.Sol[i][j-1])
			}
			// Have to check CanPair here: p.W[i][j] might be less than p.V[i][j]
			if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
				numFound.Add(numFound, p.Sol[i+1][j-1])
			}
			for k := i; k < j; k++ {
				if p.V[i][j] == p.V[i][k]+p.W[k+1][j] && p.W[k+1][j] > 0 {
					// k+1 and i are paired, so check p.Sol for bases inside that pair
					numFound.Add(numFound, new(big.Int).Mul(p.Sol[i][k], p.Sol[k+2][j-1]))
				}
			}
			p.Sol[i][j] = numFound
		}
	}
}

func (p *Predictor) CountPairings() {
	if p.Sol == nil {
		log.Fatal("Must run safecomplete.CountSolutions before safecomplete.CountPairings")
	}
	numBases := len(p.Seq.Bases)
	p.PairSafety = make([][]*big.Int, numBases)
	for i := 0; i < numBases; i++ {
		p.PairSafety[i] = make([]*big.Int, numBases)
		for j := i + 1; j < numBases; j++ {
			p.PairSafety[i][j] = new(big.Int)
		}
	}

	p.UsedBy = make([][]*big.Int, numBases)
	for i := 0; i < numBases; i++ {
		p.UsedBy[i] = make([]*big.Int, numBases)
		for j := i; j < numBases; j++ {
			p.UsedBy[i][j] = new(big.Int)
		}
	}
	p.UsedBy[0][numBases-1] = p.Sol[0][numBases-1]

	for l := numBases; l >= 2; l-- {
		for i := 0; i <= numBases-l; i++ {
			j := i + l - 1

			numSols := p.UsedBy[i][j]
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

				p.UsedBy[s.Pre.I][s.Pre.J].Add(p.UsedBy[s.Pre.I][s.Pre.J], partNumSols)
				if s.Suf != nil {
					p.UsedBy[s.Suf.I][s.Suf.J].Add(p.UsedBy[s.Suf.I][s.Suf.J], partNumSols)
				}
			}
		}
	}

	p.SingleSafety = make([]*big.Int, numBases)
	for i := 0; i < numBases; i++ {
		p.SingleSafety[i] = p.UsedBy[i][i]
	}
}
