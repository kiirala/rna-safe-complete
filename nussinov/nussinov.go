package nussinov /* import "keltainen.duckdns.org/rnafolding/nussinov" */

import "keltainen.duckdns.org/rnafolding/base"

type Predictor struct {
	Seq *base.Sequence
	V   [][]int
	W   [][]int

	MinHairpin int
}
