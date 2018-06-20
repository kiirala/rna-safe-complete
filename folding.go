package main

import "flag"
import "fmt"
import "log"
import "os"
import "strings"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/nussinov"
import "keltainen.duckdns.org/rnafolding/safecomplete"
import "keltainen.duckdns.org/rnafolding/format"

var (
	infile = flag.String("in", "", "Name of input file in FASTA format")
)

func main() {
	flag.Parse()

	f, err := os.Open(*infile)
	if err != nil {
		log.Fatalf("Unable to open file \"%s\": %v", *infile, err)
	}
	seq, err := fasta.ReadSequence(f)
	if err != nil {
		log.Fatalf("Error reading FASTA format file \"%s\": %v", infile, err)
	}
	fmt.Printf("Sequence \"%s\"\n", seq.Comment)
	fmt.Printf("Contains %d bases\n", len(seq.Bases))

	v := nussinov.FillArray(seq)
	//vcmpl := nussinov.FillComplementary(seq, v)
	flen, folding := nussinov.Backtrack(seq, v)
	fmt.Printf("Optimal folding, %d pairs: %v\n", flen, folding)
	if fPairs := countPairs(folding); fPairs != flen {
		fmt.Printf("Sanity check failed!\nOptimal folding has %d pairs, expected %d\n", fPairs, flen)
	}
	fmt.Printf(format.Folding(seq, folding))
	/*
		for i := 0; i < len(v); i++ {
			for j := i + 1; j < len(v); j++ {
				solen, sopairs := nussinov.JoinedBacktrack(seq, v, vcmpl, i, j)
				fmt.Printf("Joining [%d, %d], %d pairs\n", i, j, solen)
				fmt.Print(format.Folding(seq, sopairs))
			}
		}
	*/
	w := safecomplete.FillArray(seq, v)
	scFoldings := safecomplete.BacktrackAll(seq, v, w)
	//fmt.Printf("Safe and complete version found %d foldings\n", len(allFoldings))
	//if sanity := allFoldingsSanity(seq, allFoldings); sanity != "" {
	//	fmt.Print("Sanity check failed!\n", sanity, "\n")
	//}
	fmt.Print("Matrix v:\n", format.Matrix(v), "\n")
	fmt.Print("Matrix w:\n", format.Matrix(w), "\n")
	fmt.Printf("Found %d solutions in total\n", scFoldings.CountSolutions())
	fmt.Println(format.SCFolding(scFoldings))
	/*for i, f := range allFoldings {
		fmt.Print("\n")
		fmt.Print(f, "\n")
		if sanity := singleFoldingSanity(seq, f, flen); len(sanity) > 0 {
			fmt.Print("Sanity check failed!\n", sanity, "\n")
		}
		fmt.Printf("%s\n", allMakings[i])
		fmt.Print(format.Folding(seq, f))
	}*/
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
		if i < j && !seq.CanPair(i, j) {
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

func allFoldingsSanity(seq *base.Sequence, foldings [][]int) string {
	var errs []string
	// Check all-pairs Nussinov is subset of all solutions
	// Check that safe and complete produces same solutions as existing all-solutions
	// Check that there are no duplicate solutions
	for i := 0; i < len(foldings); i++ {
		for j := i + 1; j < len(foldings); j++ {
			if reflect.DeepEqual(foldings[i], foldings[j]) {
				errs = append(errs, fmt.Sprintf("Foldings %d and %d are exactly the same", i, j))
			}
		}
	}
	return strings.Join(errs, "\n")
}
