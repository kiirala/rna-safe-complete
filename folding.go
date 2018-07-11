package main

import "flag"
import "fmt"
import "log"
import "os"
import "sort"
import "time"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/trnadb"
import "keltainen.duckdns.org/rnafolding/safecomplete"
import "keltainen.duckdns.org/rnafolding/format"

var (
	infile     = flag.String("in", "", "Name of input file in FASTA format")
	dbfile     = flag.String("db", "", "Location of tRNA database file")
	trna       = flag.String("trna", "", "Name of tRNA sequence within database file")
	all        = flag.Bool("all", false, "Analyze all sequences in tRNA database")
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

func readTRNA() map[string]*base.Sequence {
	f, err := os.Open(*dbfile)
	if err != nil {
		log.Fatalf("Unable to open tRNA database \"%s\": %v", *dbfile, err)
	}
	seqs, err := trnadb.ReadSequences(f)
	if err != nil {
		log.Fatalf("Error reading tRNA database file \"%s\": %v", *dbfile, err)
	}
	return seqs
}

func main() {
	flag.Parse()

	var seq *base.Sequence
	var seqs map[string]*base.Sequence
	if *infile == "" && *dbfile == "" {
		log.Fatal("Please specify either -in parameter or -db parameters")
	} else if *infile != "" && *dbfile != "" {
		log.Fatal("Specify either -in or -db parameter, not both")
	} else if *infile != "" {
		seq = readFasta()
	} else {
		seqs = readTRNA()
		if *trna != "" {
			s, ok := seqs[*trna]
			if !ok {
				log.Fatalf("tRNA sequence \"%s\" not found in %s", *trna, *dbfile)
			}
			seq = s
		} else if !*all {
			log.Fatal("Specify either a tRNA sequence with -trna flag or all sequences with -all flag")
		}
	}

	if *dbfile != "" && *all {
		foldingStats(seqs)
	} else {
		singleFolding(seq)
	}
}

func singleFolding(seq *base.Sequence) {
	fmt.Printf("Sequence \"%s\"\n", seq.Comment)
	fmt.Printf("Contains %d bases\n\n", len(seq.Bases))
	fmt.Printf("Folding rules:\n  * hairpin loop must contain at least %d free bases\n\n", *minhairpin)

	nu, flen, zukerOptimals := runNussinovZuker(seq)
	fmt.Printf("Optimal folding, %d pairs: %v\n", flen, zukerOptimals[0])
	fmt.Printf("Zuker method found %d optimal solutions\n", len(zukerOptimals))
	sanityNussinovZuker(zukerOptimals, flen)
	//for _, f := range zukerOptimals {
	//	fmt.Printf(format.Folding(seq, f))
	//	fmt.Println()
	//}

	wu, wuFoldings := runWuchty(seq)
	fmt.Printf("Wuchty predictor produced %d foldings\n", len(wuFoldings))
	sanityWuchty(nu, wu, flen, zukerOptimals, wuFoldings)

	sc, scFoldings := runSafeComplete(seq, nu)
	//fmt.Println(format.Matrix(sc.Sol))

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
	sanitySafeComplete(sc, scFoldings, flen, scPairArrays, wuFoldings)

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
}

func foldingStats(seqs map[string]*base.Sequence) {
	fmt.Println("# Name NumBases FoldingPairs NumZuker NumAll NumSafeBases TimeNussinov TimeWuchty TimeSafeComplete")

	for name := range seqs {
		seq := seqs[name]

		nuStart := time.Now()
		nu, flen, zukerOptimals := runNussinovZuker(seq)
		nuTime := time.Since(nuStart)
		sanityNussinovZuker(zukerOptimals, flen)

		wuStart := time.Now()
		wu, wuFoldings := runWuchty(seq)
		wuTime := time.Since(wuStart)
		sanityWuchty(nu, wu, flen, zukerOptimals, wuFoldings)

		scStart := time.Now()
		sc, scFoldings := runSafeComplete(seq, nu)
		scFoldings.CollapseTree()
		scFoldings.LiftCommon()
		scPairArrays := scFoldings.GeneratePairArrays(seq)
		scTime := time.Since(scStart)
		sanitySafeComplete(sc, scFoldings, flen, scPairArrays, wuFoldings)

		safety := safecomplete.TrivialSafety(scPairArrays)
		numSafe := 0
		for _, s := range safety {
			if s {
				numSafe++
			}
		}
		fmt.Printf("%s %8d %12d %8d %6d %12d %12.6f %10.6f %16.6f\n",
			name, len(seq.Bases), flen, len(zukerOptimals), len(wuFoldings), numSafe,
			nuTime.Seconds(), wuTime.Seconds(), scTime.Seconds())
	}
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
