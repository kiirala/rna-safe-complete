package safecomplete

import "keltainen.duckdns.org/rnafolding/base"

type Predictor struct {
	Seq          *base.Sequence
	V            [][]int
	W            [][]int
	Sol          [][]int
	PairSafety   [][]int
	SingleSafety []int

	MinHairpin int
}
