package folding /* import "keltainen.duckdns.org/rnafolding/folding" */

import "log"
import "sort"

type FoldingPairs []int

func NewFoldingPairs(l int) FoldingPairs {
	f := make(FoldingPairs, l)
	for i := 0; i < l; i++ {
		f[i] = -9
	}
	return f
}

type FoldingSet []FoldingPairs

type FoldingOrdering struct{ FoldingSet }

func (f FoldingOrdering) Len() int {
	return len(f.FoldingSet)
}

func (f FoldingOrdering) Swap(i, j int) {
	f.FoldingSet[i], f.FoldingSet[j] = f.FoldingSet[j], f.FoldingSet[i]
}

func (f FoldingOrdering) Less(i, j int) bool {
	a := f.FoldingSet[i]
	b := f.FoldingSet[j]
	if len(a) != len(b) {
		log.Printf("arrLess: array lenths different! len(f[%d])=%d, len(f[%d])=%d", i, len(a), j, len(b))
		return len(a) < len(b)
	}
	for n := 0; n < len(a); n++ {
		if a[n] < b[n] {
			return true
		}
		if a[n] > b[n] {
			return false
		}
	}
	return false
}

func IsSubsetOf(sub FoldingSet, all FoldingSet) bool {
	ssub := make(FoldingSet, len(sub))
	copy(ssub, sub)
	sort.Sort(FoldingOrdering{ssub})
	sall := make(FoldingSet, len(all))
	copy(sall, all)
	sort.Sort(FoldingOrdering{sall})
	iall := 0
	isub := 0
	for iall < len(sall) && isub < len(ssub) {
		if FoldingArraysEqual(sall[iall], ssub[isub]) {
			isub++
		}
		iall++
	}
	if isub < len(ssub) {
		return false
	}
	return true
}

func FoldingSetsEqual(a, b FoldingSet) bool {
	if len(a) != len(b) {
		return false
	}
	sa := make(FoldingSet, len(a))
	copy(sa, a)
	sort.Sort(FoldingOrdering{sa})
	sb := make(FoldingSet, len(b))
	copy(sb, b)
	sort.Sort(FoldingOrdering{sb})
	for i := 0; i < len(sa); i++ {
		if !FoldingArraysEqual(sa[i], sb[i]) {
			log.Printf("foldingSetsEqual: difference found!\na=%v\nb=%v", sa[i], sb[i])
			return false
		}
	}
	return true
}

func FoldingArraysEqual(a, b FoldingPairs) bool {
	if len(a) != len(b) {
		log.Printf("arraysEqual: array lengths different! len(a)=%d, len(b)=%d", len(a), len(b))
		return false
	}
	for i := 0; i < len(a); i++ {
		if a[i] != b[i] && (a[i] >= 0 || b[i] >= 0) {
			return false
		}
	}
	return true
}

func FoldingInArray(f FoldingPairs, ff FoldingSet) bool {
	for _, o := range ff {
		if FoldingArraysEqual(f, o) {
			return true
		}
	}
	return false
}
