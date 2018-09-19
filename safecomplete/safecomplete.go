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
