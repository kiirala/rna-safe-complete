package safecomplete /* import "keltainen.duckdns.org/rnafolding/safecomplete" */

func (p *Predictor) FillArray() {
	numBases := len(p.Seq.Bases)

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
			if p.Seq.CanPair(i, j, p.MinHairpin) {
				if val := v[i+1][j-1] + 1; val > best {
					best = val
				}
			}
			v[i][j] = best
		}
	}

	p.V = v

	w := make([][]int, numBases)
	for i := 0; i < numBases; i++ {
		w[i] = make([]int, numBases)
	}

	for l := 2; l <= numBases; l++ {
		for i := 0; i <= numBases-l; i++ {
			j := i + l - 1
			if p.Seq.CanPair(i, j, p.MinHairpin) {
				w[i][j] = p.V[i+1][j-1] + 1
			}
		}
	}

	p.W = w
}
