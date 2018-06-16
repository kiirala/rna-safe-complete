package safecomplete

import "keltainen.duckdns.org/rnafolding/base"

func FillArray(seq *base.Sequence, v [][]int) [][]int {
	numBases := len(seq.Bases)
	w := make([][]int, numBases)
	for i := 0; i < numBases; i++ {
		w[i] = make([]int, numBases)
	}

	for l := 2; l <= numBases; l++ {
		for i := 0; i <= numBases-l; i++ {
			j := i + l - 1
			if seq.CanPair(i, j) {
				w[i][j] = v[i+1][j-1] + 1
			}
		}
	}

	return w
}
