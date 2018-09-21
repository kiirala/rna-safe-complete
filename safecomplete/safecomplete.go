package safecomplete

import "math/big"

import "keltainen.duckdns.org/rnafolding/base"

type Predictor struct {
	Seq          *base.Sequence
	V            [][]int
	W            [][]int
	Sol          [][]*big.Int
	PairSafety   [][]*big.Int
	SingleSafety []*big.Int

	MinHairpin int
}

type pair struct {
	I int
	J int
}

type split struct {
	Pre  *pair
	Suf  *pair
	Pair *pair
}

func (p *Predictor) findSplits(i, j int) []split {
	var ret []split

	if p.Seq.CanPair(i, j, p.MinHairpin) && p.V[i][j] == p.V[i+1][j-1]+1 {
		ret = append(ret, split{&pair{i + 1, j - 1}, nil, &pair{i, j}})
	}
	for k := i; k < j; k++ {
		if p.V[i][j] == p.V[i][k]+p.W[k+1][j] && p.W[k+1][j] > 0 {
			ret = append(ret, split{&pair{i, k}, &pair{k + 2, j - 1}, &pair{k + 1, j}})
		}
	}
	if p.V[i][j] == p.V[i][j-1] {
		ret = append(ret, split{&pair{i, j - 1}, &pair{j, j}, nil})
	}

	return ret
}
