// Package wuchty is all-solutions dynamic programming predictor.
// Source:
// Wuchty, Fontana, Hofacker, Schuster - Complete suboptimal folding of RNA and the stability of secondary structures
// Biopolymers, 1999, vol. 49, no. 2, pages 145–165.
package wuchty

import "log"

import "keltainen.duckdns.org/rnafolding/base"

type Predictor struct {
	Seq *base.Sequence
	V   [][]int

	MinHairpin int
	MaxStack   int
}

func (p *Predictor) FillArray() {
	numBases := len(p.Seq.Bases)
	v := make([][]int, numBases)
	for i := 0; i < numBases; i++ {
		v[i] = make([]int, numBases)
	}

	for i := numBases - 1; i >= 0; i-- {
		for j := i + 1; j < numBases; j++ {
			best := v[i][j-1]
			for l := i; l < j; l++ {
				if p.Seq.CanPair(l, j, p.MinHairpin) {
					val := 1 + v[l+1][j-1]
					if i < l-1 {
						val += v[i][l-1]
					}
					if val > best {
						best = val
					}
				}
			}
			v[i][j] = best
		}
	}

	p.V = v
}

type pair struct {
	i int
	j int
}

type state struct {
	intervals []pair
	pairs     []pair
}

func pairsToArray(pp []pair, l int) []int {
	a := make([]int, l)
	for i := 0; i < len(a); i++ {
		a[i] = -1
	}
	for _, p := range pp {
		a[p.i] = p.j
		a[p.j] = p.i
	}
	return a
}

func (p *Predictor) BacktrackAll() [][]int {
	var out [][]int
	stack := []state{
		state{
			[]pair{pair{0, len(p.Seq.Bases) - 1}},
			nil,
		},
	}

	for len(stack) > 0 {
		//log.Printf("Pop from stack, length %d", len(stack))
		var s state
		stack, s = stack[:len(stack)-1], stack[len(stack)-1]
		//log.Printf("Stack length %d", len(stack))
		if len(s.intervals) == 0 {
			out = append(out, pairsToArray(s.pairs, len(p.Seq.Bases)))
			//log.Printf("Found pairing with %d pairs", len(s.pairs))
			continue
		}
		//log.Printf("Stack length after %d", len(stack))
		numfound := 0
		var interval pair
		s.intervals, interval = s.intervals[:len(s.intervals)-1], s.intervals[len(s.intervals)-1]
		//log.Printf("Examining interval [%d, %d], %d pairs found", interval.i, interval.j, len(s.pairs))
		if interval.i < interval.j && p.V[interval.i][interval.j] == p.V[interval.i][interval.j-1] {
			stack = append(stack, state{append(s.intervals, pair{interval.i, interval.j - 1}), s.pairs})
			numfound++
		}

		//if numfound > 0 {
		//	continue
		//}

		for l := interval.i; l < interval.j; l++ {
			if !p.Seq.CanPair(l, interval.j, p.MinHairpin) {
				continue
			}
			val := 1 + p.V[l+1][interval.j-1]
			if interval.i < l-1 {
				val += p.V[interval.i][l-1]
			}
			if p.V[interval.i][interval.j] == val {
				stack = append(stack,
					state{
						append(s.intervals, pair{interval.i, l - 1}, pair{l + 1, interval.j - 1}),
						append(s.pairs, pair{l, interval.j})})
				numfound++
			}
		}

		if numfound == 0 {
			//log.Print("Default push to stack")
			stack = append(stack, s)
		}
		if len(stack) >= p.MaxStack {
			log.Printf("Stack overflow in wuchty.BacktrackAll: %d items, max %d", len(stack), p.MaxStack)
			break
		}
	}
	return out
}
