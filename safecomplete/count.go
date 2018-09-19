package safecomplete

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
