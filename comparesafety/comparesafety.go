package main /* import "keltainen.duckdns.org/rnafolding/comparesafety" */

import "encoding/json"
import "flag"
import "fmt"
import "io/ioutil"
import "log"
import "math/big"
import "os"
import "path"
import "time"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/folding"
import "keltainen.duckdns.org/rnafolding/nussinov"
import "keltainen.duckdns.org/rnafolding/safecomplete"

var (
	infile     = flag.String("in", "", "Name of input file from STRAND database")
	name       = flag.String("name", "", "Name of sequence within STRAND file")
	minhairpin = flag.Int("minhairpin", 3, "Minimum number of free bases in hairpin loop")
	single     = flag.Bool("single", false, "Solve only a single optimal folding, without safety")
	jsonout    = flag.Bool("json", false, "JSON output")
)

func main() {
	flag.Parse()

	if *infile == "" {
		fmt.Print("Please set the input file with --in flag")
		return
	}

	seqs := readFasta(*infile)
	//log.Printf("Read %d sequences from %s", len(seqs), *infile)
	if *name != "" {
		seq, ok := seqs[*name]
		if !ok {
			log.Fatalf("Sequence %s not found in %s", *name, *infile)
		}
		seqs = make(map[string]*base.Sequence)
		seqs[seq.Name] = seq
	}
	//log.Printf("Analyzing %d sequences", len(seqs))

	vals := make(chan retitem, 1)
	if *single {
		go singleAnalyze(seqs, vals)
	} else {
		go analyze(seqs, vals)
	}

	if *jsonout {
		for r := range vals {
			d, err := json.Marshal(r)
			if err != nil {
				log.Printf("Failed to write JSON: %v", err)
			}
			os.Stdout.Write(d)
		}
	} else {
		fmt.Print("# Name          NumBases NumFolds NumPairs TimeFill TimeCplx\n")
		for r := range vals {
			fmt.Printf("%-15s %8d %8s %8d %8.2f %8.2f\n",
				r.Name, r.Bases, approxint(r.NumFolds), r.NumPairs, r.Tfill, r.Tcomplex)
		}
	}
}

type retitem struct {
	Name       string
	Bases      int
	NumFolds   *big.Int
	NumPairs   int
	Tfill      float64
	Tcomplex   float64
	Tbacktrack float64
	Pairs      [][]*big.Int
	Free       []*big.Int
	Sample     folding.FoldingPairs
}

func analyze(seqs map[string]*base.Sequence, ret chan retitem) {
	for _, seq := range seqs {
		//log.Printf("Analyzing %s, with %d bases", seq.Name, len(seq.Bases))
		sFill := time.Now()
		sc := &safecomplete.Predictor{
			Seq:        seq,
			MinHairpin: *minhairpin,
		}
		sc.FillArray()

		tFill := time.Since(sFill)
		//log.Printf("Filling DP arrays done in %.0f seconds", tFill.Seconds())

		sComplex := time.Now()
		sc.CountSolutions()
		sc.CountPairings()
		/*complexSafety :=*/ sc.SafetyFromBacktrack()
		tComplex := time.Since(sComplex)
		//log.Printf("Safety done in %.0f seconds, %d solutions", tComplex.Seconds(), sc.Sol[0][len(sc.Sol[0])-1])

		ret <- retitem{
			seq.Name, len(seq.Bases), sc.Sol[0][len(seq.Bases)-1], sc.V[0][len(seq.Bases)-1],
			tFill.Seconds(), tComplex.Seconds(), 0,
			sc.PairSafety, sc.SingleSafety, nil,
		}
	}
	close(ret)
}

func singleAnalyze(seqs map[string]*base.Sequence, ret chan retitem) {
	for _, seq := range seqs {
		sFill := time.Now()
		p := &nussinov.Predictor{
			Seq:        seq,
			MinHairpin: *minhairpin,
		}
		p.FillArray()
		tFill := time.Since(sFill)

		sBt := time.Now()
		num, fold := p.Backtrack()
		tBt := time.Since(sBt)

		ret <- retitem{seq.Name, len(seq.Bases), big.NewInt(0), num,
			tFill.Seconds(), 0, tBt.Seconds(),
			nil, nil, fold}
	}
	close(ret)
}

func readFasta(fname string) map[string]*base.Sequence {
	stat, err := os.Stat(fname)
	if err != nil {
		log.Fatalf("Unable to stat %v: %v", fname, err)
	}
	if stat.IsDir() {
		ff, err := ioutil.ReadDir(fname)
		if err != nil {
			log.Fatalf("Unable to list directory %v: %v", fname, err)
		}
		ret := make(map[string]*base.Sequence)
		for _, fn := range ff {
			subfname := path.Join(fname, fn.Name())
			f, err := os.Open(subfname)
			if err != nil {
				log.Fatalf("Unable to open file \"%s\": %v", subfname, err)
			}
			seq, err := fasta.ReadSequence(f)
			if err != nil {
				log.Fatalf("Error reading FASTA format file \"%s\": %v", fn, err)
			}
			ret[seq.Name] = seq
		}
		return ret
	} else {
		f, err := os.Open(fname)
		if err != nil {
			log.Fatalf("Unable to open file \"%s\": %v", fname, err)
		}
		seq, err := fasta.ReadMultiSequence(f)
		if err != nil {
			log.Fatalf("Error reading FASTA format file \"%s\": %v", fname, err)
		}
		return seq
	}
}

func approxint(i *big.Int) string {
	if i.Cmp(big.NewInt(100000000)) < 0 {
		return i.String()
	}

	q := new(big.Int).Set(i)
	ten := big.NewInt(10)
	exp := 0
	scale := big.NewInt(1000)
	for q.Cmp(scale) >= 0 {
		exp++
		q.Div(q, ten)
	}

	return fmt.Sprintf("%.2fe%d", float64(q.Int64())/100.0, exp)
}
