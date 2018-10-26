package main

import "encoding/json"
import "flag"
import "fmt"
import "log"
import "os"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"

var (
	infile  = flag.String("in", "", "Name of input FASTA file")
	outfile = flag.String("out", "", "Name of output JSON file")
)

func main() {
	flag.Parse()

	if *infile == "" {
		fmt.Print("Please set the input file with --in flag")
		return
	}
	if *outfile == "" {
		fmt.Print("Please set the output file with --out flag")
		return
	}

	seq := readFasta(*infile)

	of, err := os.Create(*outfile)
	defer of.Close()
	if err != nil {
		log.Fatalf("Could not open %s for writing: %v", outfile, err)
		return
	}
	d, err := json.Marshal(seq)
	if err != nil {
		log.Fatalf("Could not convert sequence to JSON format: %v", err)
		return
	}
	_, err = of.Write(d)
	if err != nil {
		log.Fatalf("Failed to write to %s: %v", outfile, err)
	}
}

func readFasta(fname string) *base.Sequence {
	f, err := os.Open(fname)
	if err != nil {
		log.Fatalf("Unable to open file \"%s\": %v", fname, err)
	}
	seq, err := fasta.ReadSequence(f)
	if err != nil {
		log.Fatalf("Error reading FASTA format file \"%s\": %v", fname, err)
	}
	return seq
}
