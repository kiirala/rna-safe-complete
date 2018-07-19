package main

import "fmt"
import "log"
import "reflect"
import "sort"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/types"
import "keltainen.duckdns.org/rnafolding/nussinov"
import "keltainen.duckdns.org/rnafolding/wuchty"
import "keltainen.duckdns.org/rnafolding/safecomplete"

func runNussinovZuker(seq *base.Sequence) (*nussinov.Predictor, int, types.FoldingSet) {
	nu := &nussinov.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
	}
	nu.FillArray()
	nu.FillComplementary()
	flen, folding := nu.Backtrack()
	zukerOptimals := types.FoldingSet{folding}
	for i := 0; i < len(seq.Bases); i++ {
		for j := i + 1; j < len(seq.Bases); j++ {
			solen, sopairs := nu.JoinedBacktrack(i, j)
			if solen == flen && !types.FoldingInArray(sopairs, zukerOptimals) {
				zukerOptimals = append(zukerOptimals, sopairs)
			}
		}
	}
	return nu, flen, zukerOptimals
}

func sanityNussinovZuker(zukerOptimals types.FoldingSet, numOptimalPairs int) {
	for _, f := range zukerOptimals {
		if fPairs := countPairs(f); fPairs != numOptimalPairs {
			log.Printf("Sanity check failed! Optimal Zuker folding has %d pairs, expected %d\n", fPairs, numOptimalPairs)
		}
	}
}

func runWuchty(seq *base.Sequence) (*wuchty.Predictor, types.FoldingSet) {
	wu := &wuchty.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
		MaxStack:   10000,
	}
	wu.FillArray()
	wuFoldings := wu.BacktrackAll()
	return wu, wuFoldings
}

func sanityWuchty(nu *nussinov.Predictor, wu *wuchty.Predictor, numOptimalPairs int, zukerOptimals types.FoldingSet, wuFoldings types.FoldingSet) {
	if !reflect.DeepEqual(nu.V, wu.V) {
		log.Printf("Sanity check failed! Nussinov folding and Wuchty folding produced different DP arrays!\n")
	}
	if !types.IsSubsetOf(zukerOptimals, wuFoldings) {
		log.Printf("Sanity check failed! Zuker method found solutions that Wuchty method didn't\n")
	}
	if sanity := allFoldingsSanity(nu.Seq, wuFoldings); sanity != "" {
		log.Print("Sanity check failed!\n", sanity, "\n")
	}
	for _, f := range wuFoldings {
		if sanity := singleFoldingSanity(nu.Seq, f, numOptimalPairs); len(sanity) > 0 {
			log.Print("Sanity check failed!\n", sanity, "\n")
		}
	}
}

func runSafeComplete(seq *base.Sequence, nu *nussinov.Predictor) (*safecomplete.Predictor, *types.FoldTree) {
	sc := &safecomplete.Predictor{
		Seq:        seq,
		V:          nu.V,
		MinHairpin: *minhairpin,
	}
	sc.FillArray()

	sc.CountSolutions()

	scFoldings := sc.BacktrackAll()
	return sc, scFoldings
}

func sanitySafeComplete(sc *safecomplete.Predictor, scFoldings *types.FoldTree, numOptimalPairs int, scPairArrays types.FoldingSet, wuFoldings types.FoldingSet) {
	if numSol := scFoldings.CountSolutions(); sc.Sol[0][len(sc.Sol)-1] != numSol {
		log.Printf("Sanity check failed! Solution count matrix shows %d solutions, folding tree %d solutions", sc.Sol[0][len(sc.Sol)-1], numSol)
	}

	if sanity := allFoldingsSanity(sc.Seq, scPairArrays); sanity != "" {
		log.Print("Sanity check failed!\n", sanity, "\n")
	}
	for _, f := range scPairArrays {
		if sanity := singleFoldingSanity(sc.Seq, f, numOptimalPairs); len(sanity) > 0 {
			log.Print("Sanity check failed!\n", sanity, "\n")
		}
	}
	if !types.FoldingSetsEqual(wuFoldings, scPairArrays) {
		log.Print("Sanity check failed! Wuchty method and safe & complete method produced different foldings")
	}
}

func singleFoldingSanity(seq *base.Sequence, f types.FoldingPairs, numPairs int) string {
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
			errs = append(errs, fmt.Sprintf("Non-symmetric pair: %d (%s) -> %d (%s) but %[3]d (%s) -> %d (%s)", i, seq.Bases[i].ToCode(), j, seq.Bases[j].ToCode(), k, seq.Bases[k].ToCode()))
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

func allFoldingsSanity(seq *base.Sequence, foldings types.FoldingSet) string {
	var errs []string
	ff := make(types.FoldingSet, len(foldings))
	copy(ff, foldings)
	sort.Sort(types.FoldingOrdering{ff})
	for i := 0; i < len(foldings)-1; i++ {
		if types.FoldingArraysEqual(foldings[i], foldings[i+1]) {
			errs = append(errs, fmt.Sprintf("At least two foldings are exactly the same"))
			break
		}
	}
	return strings.Join(errs, "\n")
}
