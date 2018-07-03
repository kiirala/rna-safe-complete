package types

import "keltainen.duckdns.org/rnafolding/base"

type Pair struct {
	I int
	J int
}

type Folding struct {
	Pairs      []Pair
	Free       []int
	JoinPrefix *Folding
	JoinSuffix *Folding
	Branches   []*Folding
}

func branchless(f *Folding) bool {
	return f.JoinPrefix == nil && len(f.Branches) == 0
}

func liftBranches(p *Folding, c *Folding) {
	p.Branches = c.Branches
	p.JoinPrefix = c.JoinPrefix
	p.JoinSuffix = c.JoinSuffix
}

func (f *Folding) CollapseTree() {
	for _, b := range f.Branches {
		b.CollapseTree()
	}
	if f.JoinPrefix != nil {
		f.JoinPrefix.CollapseTree()
		f.JoinSuffix.CollapseTree()

		if branchless(f.JoinPrefix) || branchless(f.JoinSuffix) {
			pref := f.JoinPrefix
			suff := f.JoinSuffix
			if !branchless(pref) {
				liftBranches(f, pref)
			} else if !branchless(suff) {
				liftBranches(f, suff)
			}
			f.Pairs = append(f.Pairs, pref.Pairs...)
			f.Free = append(f.Free, pref.Free...)
			f.JoinPrefix = nil
			f.Pairs = append(f.Pairs, suff.Pairs...)
			f.Free = append(f.Free, suff.Free...)
			f.JoinSuffix = nil
		}
	}
}

func (f *Folding) CountSolutions() int {
	if len(f.Branches) > 0 {
		sum := 0
		for _, b := range f.Branches {
			sum += b.CountSolutions()
		}
		return sum
	}
	if f.JoinPrefix != nil {
		return f.JoinPrefix.CountSolutions() * f.JoinSuffix.CountSolutions()
	}
	return 1
}

type PairSet struct {
	s map[Pair]bool
}

type FreeSet struct {
	s map[int]bool
}

func newPairSet(pp []Pair) *PairSet {
	s := map[Pair]bool{}
	for _, p := range pp {
		s[p] = true
	}
	return &PairSet{s: s}
}

func newFreeSet(ff []int) *FreeSet {
	s := map[int]bool{}
	for _, f := range ff {
		s[f] = true
	}
	return &FreeSet{s: s}
}

func (set *PairSet) intersect(pp []Pair) *PairSet {
	s := map[Pair]bool{}
	for _, p := range pp {
		if set.s[p] {
			s[p] = true
		}
	}
	return &PairSet{s: s}
}

func (set *FreeSet) intersect(pp []int) *FreeSet {
	s := map[int]bool{}
	for _, p := range pp {
		if set.s[p] {
			s[p] = true
		}
	}
	return &FreeSet{s: s}
}

func (set *PairSet) addTo(pp []Pair) []Pair {
	ret := pp
	for p := range set.s {
		ret = append(ret, p)
	}
	return ret
}

func (set *FreeSet) addTo(pp []int) []int {
	ret := pp
	for p := range set.s {
		ret = append(ret, p)
	}
	return ret
}

func (set *PairSet) removedFrom(pp []Pair) []Pair {
	var ret []Pair
	for _, p := range pp {
		if !set.s[p] {
			ret = append(ret, p)
		}
	}
	return ret
}

func (set *FreeSet) removedFrom(pp []int) []int {
	var ret []int
	for _, p := range pp {
		if !set.s[p] {
			ret = append(ret, p)
		}
	}
	return ret
}

func (f *Folding) LiftCommon() {
	for _, b := range f.Branches {
		b.LiftCommon()
	}
	if f.JoinPrefix != nil {
		f.JoinPrefix.LiftCommon()
		f.JoinSuffix.LiftCommon()
	}

	if len(f.Branches) > 0 {
		pairs := newPairSet(f.Branches[0].Pairs)
		free := newFreeSet(f.Branches[0].Free)
		for _, b := range f.Branches {
			pairs = pairs.intersect(b.Pairs)
			free = free.intersect(b.Free)
		}
		f.Pairs = pairs.addTo(f.Pairs)
		f.Free = free.addTo(f.Free)
		for _, b := range f.Branches {
			b.Pairs = pairs.removedFrom(b.Pairs)
			b.Free = free.removedFrom(b.Free)
		}
	}
}

func joinArrays(a, b []int) []int {
	ret := make([]int, len(a))
	copy(ret, a)
	for i, x := range b {
		if x >= 0 {
			ret[i] = x
		}
	}
	return ret
}

func (f *Folding) GeneratePairArrays(seq *base.Sequence) [][]int {
	var pp [][]int
	for _, b := range f.Branches {
		pp = append(pp, b.GeneratePairArrays(seq)...)
	}
	if f.JoinPrefix != nil {
		pref := f.JoinPrefix.GeneratePairArrays(seq)
		suff := f.JoinSuffix.GeneratePairArrays(seq)
		for _, pa := range pref {
			for _, sa := range suff {
				pp = append(pp, joinArrays(pa, sa))
			}
		}
	}
	if len(pp) == 0 {
		pp = make([][]int, 1)
		pp[0] = make([]int, len(seq.Bases))
		for i := 0; i < len(pp[0]); i++ {
			pp[0][i] = -9
		}
	}
	for _, p := range pp {
		for _, pair := range f.Pairs {
			p[pair.I] = pair.J
			p[pair.J] = pair.I
		}
		for _, free := range f.Free {
			p[free] = -1
		}
	}
	return pp
}
