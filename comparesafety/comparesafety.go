package main

import "flag"
import "fmt"
import "log"
import "os"
import "reflect"
import "time"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/safecomplete"

var (
	infile     = flag.String("in", "", "Name of input file from STRAND database")
	minhairpin = flag.Int("minhairpin", 3, "Minimum number of free bases in hairpin loop")
)

func main() {
	flag.Parse()

	if *infile == "" {
		fmt.Print("Please set the input file with --in flag")
		return
	}

	seq := readFasta(*infile)

	log.Printf("Starting...")
	sFill := time.Now()
	sc := &safecomplete.Predictor{
		Seq:        seq,
		MinHairpin: *minhairpin,
	}
	sc.FillArray()

	sc.CountSolutions()
	tFill := time.Since(sFill)
	log.Printf("Filling DP arrays done in %.0f seconds, %d solutions", tFill.Seconds(), sc.Sol[0][len(sc.Sol[0])-1])

	sBacktrack := time.Now()
	folds := sc.BacktrackAll()
	tBacktrack := time.Since(sBacktrack)
	log.Printf("Prediction done in %.0f seconds", tBacktrack.Seconds())

	sTrivial := time.Now()
	trivialSafety := sc.IteratedTrivialSafety(folds)
	tTrivial := time.Since(sTrivial)
	log.Printf("Trivial safety done")

	sComplex := time.Now()
	complexSafety := sc.SafetyFromBacktrack()
	tComplex := time.Since(sComplex)
	log.Printf("Complex safety done")

	if !reflect.DeepEqual(trivialSafety, complexSafety) {
		log.Printf("Sanity check failed! IteratedTrivialSafety and SafetyFromBacktrack return different values!")
	}

	fmt.Print("# Name NumFolds TimeFill TimeBTrk TimeTriv TimeCplx\n")
	fmt.Printf("%s %8d %8.2f %8.2f %8.2f %8.2f\n", seq.Name, folds.CountSolutions(),
		tFill.Seconds(), tBacktrack.Seconds(), tTrivial.Seconds(), tComplex.Seconds())
}

func readFasta(fname string) *base.Sequence {
	f, err := os.Open(fname)
	if err != nil {
		log.Fatalf("Unable to open file \"%s\": %v", *infile, err)
	}
	seq, err := fasta.ReadSequence(f)
	if err != nil {
		log.Fatalf("Error reading FASTA format file \"%s\": %v", *infile, err)
	}
	return seq
}
