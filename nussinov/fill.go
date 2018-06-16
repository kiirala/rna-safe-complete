package nussinov

import "keltainen.duckdns.org/rnafolding/base"

func FillArray(sequence *base.Sequence) [][]int {
	numBases := len(sequence.Bases)
	v := make([][]int, numBases)
	for i := 0; i < numBases; i++ {
		v[i] = make([]int, numBases)
	}

	for i := 0; i < numBases; i++ {
		v[i][i] = 0
		if i > 0 {
			v[i-1][i] = 0
		}
	}

	for l := 1; l < numBases; l++ {
		for i := 0; i < numBases-l; i++ {
			j := i + l
			best := 0
			for k := i; k < j; k++ {
				val := v[i][k] + v[k+1][j]
				if val > best {
					best = val
				}
			}
			if sequence.CanPair(i, j) {
				if val := v[i+1][j-1] + 1; val > best {
					best = val
				}
			}
			v[i][j] = best
		}
	}

	return v
}

func FillComplementary(sequence *base.Sequence, v [][]int) [][]int {
	numBases := len(sequence.Bases)
	w := make([][]int, numBases)
	for i := 0; i < numBases; i++ {
		w[i] = make([]int, numBases)
	}
	w[0][numBases-1] = 0
	w[0][numBases-2] = 0
	w[1][numBases-1] = 0
	for l := 2; l < numBases; l++ {
		for j := numBases - l - 1; j < numBases; j++ {
			i := j + l + 1
			best := 0
			for k := j + 1; k < i-1; k++ {
				var val int
				if k < numBases {
					val = v[j+1][k] + w[i-numBases][k]
				} else {
					val = w[k+1-numBases][j] + v[k+1-numBases][i-1-numBases]
				}
				if val > best {
					best = val
				}
			}
			if sequence.CanPair((i-1)%numBases, (j+1)%numBases) {
				var val int
				if i-1 == numBases-1 {
					val = v[j+2][numBases-2] + 1
				} else if j+1 == numBases {
					val = v[1][i-2-numBases] + 1
				} else {
					val = w[i-1-numBases][j+1] + 1
				}
				if val > best {
					best = val
				}
			}
			w[i-numBases][j] = best
		}
	}
	return w
}
