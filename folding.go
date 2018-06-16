package main

import "flag"
import "fmt"
import "log"
import "os"

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
	vcmpl := nussinov.FillComplementary(seq, v)
	flen, folding := nussinov.Backtrack(seq, v)
	fmt.Printf("Optimal folding, %d pairs: %v\n", flen, folding)
	fmt.Printf(format.Folding(seq, folding))
	for i := 0; i < len(v); i++ {
		for j := i + 1; j < len(v); j++ {
			solen, sopairs := nussinov.JoinedBacktrack(seq, v, vcmpl, i, j)
			fmt.Printf("Joining [%d, %d], %d pairs\n", i, j, solen)
			fmt.Print(format.Folding(seq, sopairs))
		}
	}
	w := safecomplete.FillArray(seq, v)
	allFoldings, allMakings := safecomplete.BacktrackAll(seq, v, w)
	fmt.Printf("Safe and complete version found %d foldings\n", len(allFoldings))
	fmt.Print(format.Matrix(allFoldings))
	for i, f := range allFoldings {
		fmt.Print("\n")
		fmt.Printf("%s\n", allMakings[i])
		fmt.Print(format.Folding(seq, f))
	}
}
