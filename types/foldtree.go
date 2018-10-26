package types /* import "keltainen.duckdns.org/rnafolding/types" */

import "fmt"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

type Pair struct {
	I int
	J int
}

type FoldTree struct {
	Pairs      []Pair
	Free       []int
	JoinPrefix *FoldTree
	JoinSuffix *FoldTree
	Branches   []*FoldTree
}

func branchless(f *FoldTree) bool {
	return f.JoinPrefix == nil && len(f.Branches) == 0
}

func liftBranches(p *FoldTree, c *FoldTree) {
	p.Branches = c.Branches
	p.JoinPrefix = c.JoinPrefix
	p.JoinSuffix = c.JoinSuffix
}

func (f *FoldTree) CollapseTree() {
	for _, b := range f.Branches {
		b.CollapseTree()
	}
	if f.JoinPrefix != nil {
		f.JoinPrefix.CollapseTree()
		f.JoinSuffix.CollapseTree()
		f.Pairs = append(f.Pairs, f.JoinPrefix.Pairs...)
		f.JoinPrefix.Pairs = nil
		f.Free = append(f.Free, f.JoinPrefix.Free...)
		f.JoinPrefix.Free = nil
		f.Pairs = append(f.Pairs, f.JoinSuffix.Pairs...)
		f.JoinSuffix.Pairs = nil
		f.Free = append(f.Free, f.JoinSuffix.Free...)
		f.JoinSuffix.Free = nil

		if branchless(f.JoinPrefix) || branchless(f.JoinSuffix) {
			pref := f.JoinPrefix
			suff := f.JoinSuffix
			f.JoinPrefix = nil
			f.JoinSuffix = nil
			if !branchless(pref) {
				liftBranches(f, pref)
			} else if !branchless(suff) {
				liftBranches(f, suff)
			}
		}
	}
}

func (f *FoldTree) CountSolutions() int {
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

func (f *FoldTree) LiftCommon() {
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

func (f *FoldTree) SimplifyImmediate() {
	if f.JoinPrefix != nil {
		f.Pairs = append(f.Pairs, f.JoinPrefix.Pairs...)
		f.JoinPrefix.Pairs = nil
		f.Free = append(f.Free, f.JoinPrefix.Free...)
		f.JoinPrefix.Free = nil
		f.Pairs = append(f.Pairs, f.JoinSuffix.Pairs...)
		f.JoinSuffix.Pairs = nil
		f.Free = append(f.Free, f.JoinSuffix.Free...)
		f.JoinSuffix.Free = nil

		if branchless(f.JoinPrefix) || branchless(f.JoinSuffix) {
			pref := f.JoinPrefix
			suff := f.JoinSuffix
			f.JoinPrefix = nil
			f.JoinSuffix = nil
			if !branchless(pref) {
				liftBranches(f, pref)
			} else if !branchless(suff) {
				liftBranches(f, suff)
			}
		}
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

func joinArrays(a, b folding.FoldingPairs) []int {
	ret := make([]int, len(a))
	copy(ret, a)
	for i, x := range b {
		if x >= 0 {
			ret[i] = x
		}
	}
	return ret
}

func (f *FoldTree) GeneratePairArrays(seq *base.Sequence) folding.FoldingSet {
	var pp folding.FoldingSet
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
		pp = make(folding.FoldingSet, 1)
		pp[0] = folding.NewFoldingPairs(len(seq.Bases))
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

func (f *FoldTree) PairArrayIterator(seq *base.Sequence) <-chan folding.FoldingPairs {
	fold := folding.NewFoldingPairs(len(seq.Bases))
	return f.pairArrayIteratorInternal(seq, fold)
}

func (f *FoldTree) pairArrayIteratorInternal(seq *base.Sequence, fold folding.FoldingPairs) <-chan folding.FoldingPairs {
	ch := make(chan folding.FoldingPairs, 10)
	go f.recursivePAIEntry(seq, fold, ch)
	return ch
}

func (f *FoldTree) recursivePAIEntry(seq *base.Sequence, fold folding.FoldingPairs, ch chan folding.FoldingPairs) {
	f.recursivePairArrayIterator(seq, fold, ch)
	close(ch)
}

func (f *FoldTree) recursivePairArrayIterator(seq *base.Sequence, fold folding.FoldingPairs, ch chan folding.FoldingPairs) {
	for _, pair := range f.Pairs {
		fold[pair.I] = pair.J
		fold[pair.J] = pair.I
	}
	for _, free := range f.Free {
		fold[free] = -1
	}

	if len(f.Branches) == 0 && f.JoinPrefix == nil {
		ch <- fold
	}
	for _, b := range f.Branches {
		fcopy := folding.NewFoldingPairs(len(seq.Bases))
		copy(fcopy, fold)
		b.recursivePairArrayIterator(seq, fcopy, ch)
	}
	if f.JoinPrefix != nil {
		fcopy := folding.NewFoldingPairs(len(seq.Bases))
		copy(fcopy, fold)
		for pa := range f.JoinPrefix.pairArrayIteratorInternal(seq, fcopy) {
			fnew := folding.NewFoldingPairs(len(seq.Bases))
			for sa := range f.JoinSuffix.pairArrayIteratorInternal(seq, fnew) {
				ch <- joinArrays(pa, sa)
			}
		}
	}
}

func (f *FoldTree) String() string {
	return recursiveSCFolding(f, 0)
}

func recursiveSCFolding(f *FoldTree, depth int) string {
	out := ""
	indent := ""
	for i := 0; i < depth; i++ {
		indent += "    "
	}
	if len(f.Pairs) > 0 {
		out += indent + "Pairs: "
		var pairs []string
		for _, p := range f.Pairs {
			pairs = append(pairs, fmt.Sprintf("(%d,%d)", p.I, p.J))
		}
		out += strings.Join(pairs, ", ") + "\n"
	}
	if len(f.Free) > 0 {
		out += indent + "Free: "
		var free []string
		for _, i := range f.Free {
			free = append(free, fmt.Sprintf("%d", i))
		}
		out += strings.Join(free, ", ") + "\n"
	}
	if f.JoinPrefix != nil {
		out += indent + "Join prefix:\n"
		out += recursiveSCFolding(f.JoinPrefix, depth+1)
		out += indent + "Join suffix:\n"
		out += recursiveSCFolding(f.JoinSuffix, depth+1)
	}
	for i, b := range f.Branches {
		out += fmt.Sprintf("%sAlternative %d of %d:\n", indent, i+1, len(f.Branches))
		out += recursiveSCFolding(b, depth+1)
	}
	return out
}
