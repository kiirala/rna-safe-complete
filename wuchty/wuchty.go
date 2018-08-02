// Package wuchty is all-solutions dynamic programming predictor.
// Source:
// Wuchty, Fontana, Hofacker, Schuster - Complete suboptimal folding of RNA and the stability of secondary structures
// Biopolymers, 1999, vol. 49, no. 2, pages 145–165.
package wuchty

import "log"
import "fmt"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

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
	structure string
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

func refinedStructure(s state, i, p []pair, str string) state {
	ints := make([]pair, len(s.intervals), len(s.intervals)+len(i))
	copy(ints, s.intervals)
	ints = append(ints, i...)
	pairs := make([]pair, len(s.pairs), len(s.pairs)+len(p))
	copy(pairs, s.pairs)
	pairs = append(pairs, p...)
	return state{ints, pairs, s.structure + str}
}

func (p *Predictor) BacktrackAll() folding.FoldingSet {
	var out folding.FoldingSet
	stack := []state{
		state{
			[]pair{pair{0, len(p.Seq.Bases) - 1}},
			nil,
			"",
		},
	}

	for len(stack) > 0 {
		var s state
		stack, s = stack[:len(stack)-1], stack[len(stack)-1]
		if len(s.intervals) == 0 {
			out = append(out, pairsToArray(s.pairs, len(p.Seq.Bases)))
			//log.Println(s.structure)
			continue
		}
		numfound := 0
		var interval pair
		s.intervals, interval = s.intervals[:len(s.intervals)-1], s.intervals[len(s.intervals)-1]
		if interval.i < interval.j && p.V[interval.i][interval.j] == p.V[interval.i][interval.j-1] {
			stack = append(stack, refinedStructure(s, []pair{pair{interval.i, interval.j - 1}}, nil, fmt.Sprintf("d(%d,%d) ", interval.i, interval.j-1)))
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
				stack = append(stack, refinedStructure(
					s,
					[]pair{pair{interval.i, l - 1}, pair{l + 1, interval.j - 1}},
					[]pair{pair{l, interval.j}},
					fmt.Sprintf("p(%d,%d) j(%d,%d)(%d,%d) ", l, interval.j, interval.i, l-1, l+1, interval.j-1)))
				numfound++
			}
		}

		if numfound == 0 {
			stack = append(stack, state{s.intervals, s.pairs, s.structure + fmt.Sprintf("pop(%d,%d) ", interval.i, interval.j)})
		}
		if len(stack) >= p.MaxStack {
			log.Printf("Stack overflow in wuchty.BacktrackAll: %d items, max %d", len(stack), p.MaxStack)
			break
		}
	}
	return out
}
