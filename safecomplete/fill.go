package safecomplete

func (p *Predictor) FillArray() {
	numBases := len(p.Seq.Bases)
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
