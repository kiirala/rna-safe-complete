package main

import "flag"
import "fmt"
import "log"
import "os"
import "sort"
import "strings"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/trnadb"
import "keltainen.duckdns.org/rnafolding/nussinov"
import "keltainen.duckdns.org/rnafolding/wuchty"
import "keltainen.duckdns.org/rnafolding/safecomplete"
import "keltainen.duckdns.org/rnafolding/format"

var (
	infile     = flag.String("in", "", "Name of input file in FASTA format")
	dbfile     = flag.String("db", "", "Location of tRNA database file")
	trna       = flag.String("trna", "", "Name of tRNA sequence within database file")
	minhairpin = flag.Int("minhairpin", 3, "Minimum number of free bases in hairpin loop")
)

func readFasta() *base.Sequence {
	f, err := os.Open(*infile)
	if err != nil {
		log.Fatalf("Unable to open file \"%s\": %v", *infile, err)
	}
	seq, err := fasta.ReadSequence(f)
	if err != nil {
		log.Fatalf("Error reading FASTA format file \"%s\": %v", *infile, err)
	}
	return seq
}

func readTRNA() *base.Sequence {
	f, err := os.Open(*dbfile)
	if err != nil {
		log.Fatalf("Unable to open tRNA database \"%s\": %v", *dbfile, err)
	}
	seq, err := trnadb.ReadSequence(f, *trna)
	if err != nil {
		log.Fatalf("Error reading tRNA database file \"%s\": %v", *dbfile, err)
	}
	return seq
}

func main() {
	flag.Parse()

	var seq *base.Sequence
	if *infile == "" && *dbfile == "" {
		log.Fatal("Please specify either -in parameter or -db and -trna parameters")
	} else if *infile != "" && *dbfile != "" {
		log.Fatal("Specify either -in or -db parameter, not both")
	} else if *infile != "" {
		seq = readFasta()
	} else {
		seq = readTRNA()
	}
	fmt.Printf("Sequence \"%s\"\n", seq.Comment)
	fmt.Printf("Contains %d bases\n\n", len(seq.Bases))
	fmt.Printf("Folding rules:\n  * hairpin loop must contain at least %d free bases\n\n", *minhairpin)

	nu := &nussinov.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
	}
	nu.FillArray()
	nu.FillComplementary()
	flen, folding := nu.Backtrack()
	fmt.Printf("Optimal folding, %d pairs: %v\n", flen, folding)
	if fPairs := countPairs(folding); fPairs != flen {
		log.Printf("Sanity check failed!\nOptimal folding has %d pairs, expected %d\n", fPairs, flen)
	}
	zukerOptimals := [][]int{folding}
	for i := 0; i < len(seq.Bases); i++ {
		for j := i + 1; j < len(seq.Bases); j++ {
			solen, sopairs := nu.JoinedBacktrack(i, j)
			if solen == flen && !foldingInArray(sopairs, zukerOptimals) {
				zukerOptimals = append(zukerOptimals, sopairs)
			}
		}
	}
	fmt.Printf("Zuker method found %d optimal solutions\n", len(zukerOptimals))
	//for _, f := range zukerOptimals {
	//	fmt.Printf(format.Folding(seq, f))
	//	fmt.Println()
	//}

	wu := &wuchty.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
		MaxStack:   10000,
	}
	wu.FillArray()
	if !reflect.DeepEqual(nu.V, wu.V) {
		log.Printf("Sanity check failed! Nussinov folding and Wuchty folding produced different DP arrays!\n")
	}
	wuFoldings := wu.BacktrackAll()
	fmt.Printf("Wuchty predictor produced %d foldings\n", len(wuFoldings))
	if !isSubsetOf(zukerOptimals, wuFoldings) {
		log.Printf("Sanity check failed! Zuker method found solutions that Wuchty method didn't\n")
	}
	if sanity := allFoldingsSanity(seq, wuFoldings); sanity != "" {
		log.Print("Sanity check failed!\n", sanity, "\n")
	}
	for _, f := range wuFoldings {
		if sanity := singleFoldingSanity(seq, f, flen); len(sanity) > 0 {
			log.Print("Sanity check failed!\n", sanity, "\n")
			log.Println(f)
		}
		//fmt.Println(format.Folding(seq, f))
	}

	sc := &safecomplete.Predictor{
		Seq:        seq,
		V:          nu.V,
		MinHairpin: *minhairpin,
	}
	sc.FillArray()

	sc.CountSolutions()
	fmt.Println(format.Matrix(sc.Sol))

	scFoldings := sc.BacktrackAll()
	if numSol := scFoldings.CountSolutions(); sc.Sol[0][len(sc.Sol)-1] != numSol {
		log.Printf("Sanity check failed! Solution count matrix shows %d solutions, folding tree %d solutions", sc.Sol[0][len(sc.Sol)-1], numSol)
	}

	//fmt.Print("Matrix v:\n", format.Matrix(v), "\n")
	//fmt.Print("Matrix w:\n", format.Matrix(w), "\n")
	fmt.Printf("Found %d solutions in total\n", scFoldings.CountSolutions())
	//fmt.Println(scFoldings)
	scFoldings.CollapseTree()
	fmt.Printf("Found %d solutions after CollapseTree\n", scFoldings.CountSolutions())
	//fmt.Println(scFoldings)
	scFoldings.LiftCommon()
	fmt.Printf("Found %d solutions after LiftCommon\n", scFoldings.CountSolutions())

	scPairArrays := scFoldings.GeneratePairArrays(seq)
	fmt.Printf("Folding tree -> folding arrays conversion produced %d foldings\n", len(scPairArrays))
	if sanity := allFoldingsSanity(seq, scPairArrays); sanity != "" {
		log.Print("Sanity check failed!\n", sanity, "\n")
	}
	for _, f := range scPairArrays {
		if sanity := singleFoldingSanity(seq, f, flen); len(sanity) > 0 {
			log.Print("Sanity check failed!\n", sanity, "\n")
			log.Println(f)
			log.Print(format.Folding(seq, f))
		}
	}
	if !foldingSetsEqual(wuFoldings, scPairArrays) {
		log.Print("Sanity check failed! Wuchty method and safe & complete method produced different foldings")
	}
	fmt.Println(scFoldings)
	safety := safecomplete.TrivialSafety(scPairArrays)
	fmt.Printf(format.FoldingWithSafety(seq, scPairArrays[0], safety))
	numSafe := 0
	for _, s := range safety {
		if s {
			numSafe++
		}
	}
	fmt.Printf("Safe bases %d/%d (%f %%)\n", numSafe, len(safety), float64(numSafe*100)/float64(len(safety)))

	/*
		for _, f := range scPairArrays {
			fmt.Print("\n")
			fmt.Print(f, "\n")
			if sanity := singleFoldingSanity(seq, f, flen); len(sanity) > 0 {
				fmt.Print("Sanity check failed!\n", sanity, "\n")
			}
			fmt.Print(format.Folding(seq, f))
		}
	*/
}

func singleFoldingSanity(seq *base.Sequence, f []int, numPairs int) string {
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

func countPairs(f []int) int {
	c := 0
	for _, n := range f {
		if n >= 0 {
			c++
		}
	}
	return c / 2
}

type FoldingSet [][]int

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

func isSubsetOf(sub [][]int, all [][]int) bool {
	ssub := make([][]int, len(sub))
	copy(ssub, sub)
	sort.Sort(FoldingOrdering{ssub})
	sall := make([][]int, len(all))
	copy(sall, all)
	sort.Sort(FoldingOrdering{sall})
	iall := 0
	isub := 0
	for iall < len(sall) && isub < len(ssub) {
		if foldingArraysEqual(sall[iall], ssub[isub]) {
			isub++
		}
		iall++
	}
	if isub < len(ssub) {
		return false
	}
	return true
}

func foldingSetsEqual(a [][]int, b [][]int) bool {
	if len(a) != len(b) {
		return false
	}
	sa := make([][]int, len(a))
	copy(sa, a)
	sort.Sort(FoldingOrdering{sa})
	sb := make([][]int, len(b))
	copy(sb, b)
	sort.Sort(FoldingOrdering{sb})
	for i := 0; i < len(sa); i++ {
		if !foldingArraysEqual(sa[i], sb[i]) {
			log.Printf("foldingSetsEqual: difference found!\na=%v\nb=%v", sa[i], sb[i])
			return false
		}
	}
	return true
}

func allFoldingsSanity(seq *base.Sequence, foldings [][]int) string {
	var errs []string
	for i := 0; i < len(foldings); i++ {
		for j := i + 1; j < len(foldings); j++ {
			if reflect.DeepEqual(foldings[i], foldings[j]) {
				errs = append(errs, fmt.Sprintf("Foldings %d and %d are exactly the same", i, j))
			}
		}
	}
	return strings.Join(errs, "\n")
}

func foldingArraysEqual(a, b []int) bool {
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

func foldingInArray(f []int, ff [][]int) bool {
	for _, o := range ff {
		if foldingArraysEqual(f, o) {
			return true
		}
	}
	return false
}
