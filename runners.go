package main /* import "keltainen.duckdns.org/rnafolding" */

import "fmt"
import "math/big"
import "reflect"
import "sort"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"
import "keltainen.duckdns.org/rnafolding/types"
import "keltainen.duckdns.org/rnafolding/nussinov"
import "keltainen.duckdns.org/rnafolding/wuchty"
import "keltainen.duckdns.org/rnafolding/safecomplete"

func runNussinovZuker(seq *base.Sequence) (*nussinov.Predictor, int, folding.FoldingSet) {
	nu := &nussinov.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
	}
	nu.FillArray()
	nu.FillComplementary()
	flen, fold := nu.Backtrack()
	zukerOptimals := folding.FoldingSet{fold}
	for i := 0; i < len(seq.Bases); i++ {
		for j := i + 1; j < len(seq.Bases); j++ {
			solen, sopairs := nu.JoinedBacktrack(i, j)
			if solen == flen && !folding.FoldingInArray(sopairs, zukerOptimals) {
				zukerOptimals = append(zukerOptimals, sopairs)
			}
		}
	}
	return nu, flen, zukerOptimals
}

func sanityNussinovZuker(zukerOptimals folding.FoldingSet, numOptimalPairs int) string {
	var out []string
	for _, f := range zukerOptimals {
		if fPairs := countPairs(f); fPairs != numOptimalPairs {
			out = append(out,
				fmt.Sprintf("Sanity check failed! Optimal Zuker folding has %d pairs, expected %d", fPairs, numOptimalPairs))
		}
	}
	return strings.Join(out, "\n")
}

func runWuchty(seq *base.Sequence) (*wuchty.Predictor, folding.FoldingSet) {
	wu := &wuchty.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
		MaxStack:   10000,
	}
	wu.FillArray()
	wuFoldings := wu.BacktrackAll()
	return wu, wuFoldings
}

func sanityWuchty(nu *nussinov.Predictor, wu *wuchty.Predictor, numOptimalPairs int, zukerOptimals folding.FoldingSet, wuFoldings folding.FoldingSet) string {
	var out []string
	if !reflect.DeepEqual(nu.V, wu.V) {
		out = append(out, fmt.Sprintf("Sanity check failed! Nussinov folding and Wuchty folding produced different DP arrays!"))
	}
	if !folding.IsSubsetOf(zukerOptimals, wuFoldings) {
		out = append(out, fmt.Sprintf("Sanity check failed! Zuker method found solutions that Wuchty method didn't"))
	}
	if sanity := allFoldingsSanity(nu.Seq, wuFoldings); sanity != "" {
		out = append(out, fmt.Sprint("Sanity check failed!\n", sanity))
	}
	for _, f := range wuFoldings {
		if sanity := singleFoldingSanity(nu.Seq, f, numOptimalPairs); len(sanity) > 0 {
			out = append(out, fmt.Sprint("Sanity check failed!\n", sanity))
		}
	}
	return strings.Join(out, "\n")
}

func runSafeComplete(seq *base.Sequence) (*safecomplete.Predictor, *types.FoldTree) {
	sc := &safecomplete.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
	}
	sc.FillArray()

	scFoldings := sc.BacktrackFolding()
	return sc, scFoldings
}

func sanitySafeComplete(sc *safecomplete.Predictor, scFoldings *types.FoldTree, numOptimalPairs int, scPairArrays folding.FoldingSet, wuFoldings folding.FoldingSet) string {
	var out []string
	if numSol := scFoldings.CountSolutions(); sc.Sol[0][len(sc.Sol)-1].Cmp(big.NewInt(int64(numSol))) != 0 {
		out = append(out, fmt.Sprintf("Sanity check failed! Solution count matrix shows %d solutions, folding tree %d solutions", sc.Sol[0][len(sc.Sol)-1], numSol))
	}

	if sanity := allFoldingsSanity(sc.Seq, scPairArrays); sanity != "" {
		out = append(out, fmt.Sprint("Sanity check failed!\n", sanity, "\n"))
	}
	for _, f := range scPairArrays {
		if sanity := singleFoldingSanity(sc.Seq, f, numOptimalPairs); len(sanity) > 0 {
			out = append(out, fmt.Sprint("Sanity check failed!\n", sanity, "\n"))
		}
	}
	if !folding.FoldingSetsEqual(wuFoldings, scPairArrays) {
		out = append(out, fmt.Sprint("Sanity check failed! Wuchty method and safe & complete method produced different foldings"))
	}
	return strings.Join(out, "\n")
}

func singleFoldingSanity(seq *base.Sequence, f folding.FoldingPairs, numPairs int) string {
	var errs []string
	for i, j := range f {
		if j < 0 {
			continue
		}
		if i == j {
			errs = append(errs, fmt.Sprintf("Base %d (%s) is paired with itself", i, seq.Bases[i].ToCode()))
		}
		if f[j] != i {
			k := f[j]
			if k >= 0 && k < len(seq.Bases) {
				errs = append(errs, fmt.Sprintf("Non-symmetric pair: %d (%s) -> %d (%s) but %[3]d (%s) -> %d (%s)", i, seq.Bases[i].ToCode(), j, seq.Bases[j].ToCode(), k, seq.Bases[k].ToCode()))
			} else {
				errs = append(errs, fmt.Sprintf("Non-symmetric pair: %d (%s) -> %d (%s) but %[3]d (%s) -> %d (outside sequence)", i, seq.Bases[i].ToCode(), j, seq.Bases[j].ToCode(), k))
			}
		}
		if i < j && !seq.CanPair(i, j, *minhairpin) {
			errs = append(errs, fmt.Sprintf("%d (%s) and %d (%s) are paired, but not a valid base pair", i, seq.Bases[i].ToCode(), j, seq.Bases[j].ToCode()))
		}
	}
	if fPairs := countPairs(f); fPairs != numPairs {
		errs = append(errs, fmt.Sprintf("Folding has %d pairs, expected %d", fPairs, numPairs))
	}
	return strings.Join(errs, "\n")
}

func allFoldingsSanity(seq *base.Sequence, foldings folding.FoldingSet) string {
	var errs []string
	ff := make(folding.FoldingSet, len(foldings))
	copy(ff, foldings)
	sort.Sort(folding.FoldingOrdering{ff})
	for i := 0; i < len(foldings)-1; i++ {
		if folding.FoldingArraysEqual(foldings[i], foldings[i+1]) {
			errs = append(errs, fmt.Sprintf("At least two foldings are exactly the same"))
			break
		}
	}
	return strings.Join(errs, "\n")
}
