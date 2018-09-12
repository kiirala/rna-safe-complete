package main

import "encoding/json"
import "flag"
import "fmt"
import "io/ioutil"
import "log"
import "os"
import "path"
import "reflect"
import "time"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/folding"
import "keltainen.duckdns.org/rnafolding/trnadb"
import "keltainen.duckdns.org/rnafolding/safecomplete"
import "keltainen.duckdns.org/rnafolding/format"

var (
	infile     = flag.String("in", "", "Name of input file in FASTA format")
	stranddir  = flag.String("stranddir", "", "Directory containing a number of dot-bracket files from STRAND DB")
	dbfile     = flag.String("db", "", "Location of tRNA database file")
	trna       = flag.String("trna", "", "Name of tRNA sequence within database file")
	all        = flag.Bool("all", false, "Analyze all sequences in tRNA database")
	minhairpin = flag.Int("minhairpin", 3, "Minimum number of free bases in hairpin loop")
	outdir     = flag.String("outdir", "", "Write result files to given directory")
)

type Rules struct {
	MinHairpin int
}

type Timing struct {
	ZukerSeconds              float64
	WuchtySeconds             float64
	SafeCompleteSeconds       float64
	PairArraysSeconds         float64
	SafeCompleteTotalSeconds  float64
	TrivialSafetySeconds      float64
	SafeCompleteSafetySeconds float64
}

type Counts struct {
	SequenceBases          int
	OptimalPairs           int
	ZukerFoldings          int
	WuchtyFoldings         int
	SafeCompleteFoldings   int
	AfterCollapseTree      int
	AfterLiftCommon        int
	SafeCompletePairArrays int
	SafeBases              int
}

type Sanity struct {
	Zuker        string
	Wuchty       string
	SafeComplete string
	Safety       string
}

type Folding struct {
	Pairing    folding.FoldingPairs
	PairCount  int
	DotBracket string
	AsciiArt   string
}

type Safety struct {
	SafeBase  []bool
	PairCount [][]int
	FreeCount []int
}

type OutputEntry struct {
	Name                    string
	Comment                 string
	Sequence                string
	SequenceWithSafety      string
	Rules                   Rules
	Timing                  Timing
	Counts                  Counts
	Sanity                  Sanity
	ReferenceFolding        Folding
	ZukerFoldings           []Folding
	AllFoldings             []Folding
	SafeCompleteFoldingTree string
	Safety                  Safety
	ReferencePosition       []int
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

func readStrandDir() map[string]*base.Sequence {
	ff, err := ioutil.ReadDir(*stranddir)
	if err != nil {
		log.Fatalf("Unable to list directory %v: %v", *stranddir, err)
	}
	ret := make(map[string]*base.Sequence)
	for _, f := range ff {
		seq := readFasta(path.Join(*stranddir, f.Name()))
		ret[seq.Name] = seq
	}
	return ret
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
	if *infile == "" && *dbfile == "" && *stranddir == "" {
		log.Fatal("Please specify either -in parameter or -db parameters")
	} else if *infile != "" && *dbfile != "" {
		log.Fatal("Specify either -in or -db parameter, not both")
	} else if *infile != "" {
		seq = readFasta(*infile)
	} else if *stranddir != "" {
		seqs = readStrandDir()
	} else {
		seqs = readTRNA()
	}

	if seq != nil {
		// Do nothing
	} else if *trna != "" {
		s, ok := seqs[*trna]
		if !ok {
			log.Fatalf("tRNA sequence \"%s\" not found in %s", *trna, *dbfile)
		}
		seq = s
	} else if !*all {
		log.Fatal("Specify either a tRNA sequence with -trna flag or all sequences with -all flag")
	}

	out := make(chan OutputEntry, 1)
	if seqs != nil && *all {
		go foldingStats(seqs, out)
	} else {
		o := singleFolding(seq)
		out <- o
		close(out)
	}

	if *outdir != "" {
		for o := range out {
			fname := path.Join(*outdir, o.Name+".json")
			f, err := os.Create(fname)
			if err != nil {
				log.Printf("Could not open %s for writing: %v", fname, err)
				continue
			}
			defer f.Close()
			b, err := json.Marshal(o)
			if err != nil {
				log.Print("Count not encode JSON for %s: %v", o.Name, err)
				continue
			}
			_, err = f.Write(b)
			if err != nil {
				log.Print("Failed to write to %s: %v", fname, err)
			}
		}
	}
}

func foldingsToOutputFormat(seq *base.Sequence, ff folding.FoldingSet, safety []bool) []Folding {
	if len(ff) == 1 && ff[0] == nil {
		// Escape hatch for reference foldings, when one does not exist
		return []Folding{{}}
	}
	out := make([]Folding, len(ff))
	for i, f := range ff {
		out[i] = Folding{
			Pairing:    f,
			PairCount:  countPairs(f),
			DotBracket: format.DotBracket(f),
			AsciiArt:   format.FoldingWithSafety(seq, f, safety),
		}
	}
	return out
}

func foldSequence(seq *base.Sequence) OutputEntry {
	nuStart := time.Now()
	nu, flen, zukerOptimals := runNussinovZuker(seq)
	nuTime := time.Since(nuStart)

	wuStart := time.Now()
	wu, wuFoldings := runWuchty(seq)
	wuTime := time.Since(wuStart)

	scStart := time.Now()
	sc, scFoldings := runSafeComplete(seq)
	scCountOriginal := scFoldings.CountSolutions()
	scFoldings.CollapseTree()
	scCountCollapsed := scFoldings.CountSolutions()
	scFoldings.LiftCommon()
	scCountLifted := scFoldings.CountSolutions()
	scTime := time.Since(scStart)

	paStart := time.Now()
	scPairArrays := scFoldings.GeneratePairArrays(seq)
	paTime := time.Since(paStart)
	scTotalTime := time.Since(scStart)

	trsaStart := time.Now()
	safety := safecomplete.TrivialSafety(scPairArrays)
	trsaTime := time.Since(trsaStart)
	newsaStart := time.Now()
	newSafety := sc.SafetyFromBacktrack()
	newsaTime := time.Since(newsaStart)

	var safetySanity string
	if !reflect.DeepEqual(safety, newSafety) {
		safetySanity = "Sanity check failed! TrivialSafety and SafetyFromBacktrack return different values!"
	}
	numSafe := 0
	for _, s := range safety {
		if s {
			numSafe++
		}
	}

	return OutputEntry{
		Name:               seq.Name,
		Comment:            seq.Comment,
		Sequence:           seq.BasesString(),
		SequenceWithSafety: seq.BasesSafetyString(safety),
		Rules: Rules{
			MinHairpin: *minhairpin,
		},
		Timing: Timing{
			ZukerSeconds:              nuTime.Seconds(),
			WuchtySeconds:             wuTime.Seconds(),
			SafeCompleteSeconds:       scTime.Seconds(),
			PairArraysSeconds:         paTime.Seconds(),
			SafeCompleteTotalSeconds:  scTotalTime.Seconds(),
			TrivialSafetySeconds:      trsaTime.Seconds(),
			SafeCompleteSafetySeconds: newsaTime.Seconds(),
		},
		Counts: Counts{
			SequenceBases:          len(seq.Bases),
			OptimalPairs:           flen,
			ZukerFoldings:          len(zukerOptimals),
			WuchtyFoldings:         len(wuFoldings),
			SafeCompleteFoldings:   scCountOriginal,
			AfterCollapseTree:      scCountCollapsed,
			AfterLiftCommon:        scCountLifted,
			SafeCompletePairArrays: len(scPairArrays),
			SafeBases:              numSafe,
		},
		Sanity: Sanity{
			Zuker:        sanityNussinovZuker(zukerOptimals, flen),
			Wuchty:       sanityWuchty(nu, wu, flen, zukerOptimals, wuFoldings),
			SafeComplete: sanitySafeComplete(sc, scFoldings, flen, scPairArrays, wuFoldings),
			Safety:       safetySanity,
		},
		ReferenceFolding:        foldingsToOutputFormat(seq, folding.FoldingSet{seq.ReferenceFolding}, safety)[0],
		ZukerFoldings:           foldingsToOutputFormat(seq, zukerOptimals, safety),
		AllFoldings:             foldingsToOutputFormat(seq, scPairArrays, safety),
		SafeCompleteFoldingTree: scFoldings.String(),
		Safety: Safety{
			SafeBase:  safety,
			PairCount: sc.PairSafety,
			FreeCount: sc.SingleSafety,
		},
		ReferencePosition: seq.ReferencePosition,
	}
}

func singleFolding(seq *base.Sequence) OutputEntry {
	o := foldSequence(seq)

	fmt.Printf("Sequence \"%s\"\n", o.Comment)
	fmt.Printf("Contains %d bases\n\n", o.Counts.SequenceBases)
	fmt.Printf("Folding rules:\n  * hairpin loop must contain at least %d free bases\n\n", o.Rules.MinHairpin)

	fmt.Printf("Optimal folding, %d pairs: %v\n", o.Counts.OptimalPairs, o.ZukerFoldings[0].Pairing)
	fmt.Printf("Zuker method found %d optimal solutions\n", o.Counts.ZukerFoldings)
	if o.Sanity.Zuker != "" {
		log.Print(o.Sanity.Zuker)
	}

	fmt.Printf("Wuchty predictor produced %d foldings\n", o.Counts.WuchtyFoldings)
	if o.Sanity.Wuchty != "" {
		log.Print(o.Sanity.Wuchty)
	}

	fmt.Printf("Found %d solutions in total\n", o.Counts.SafeCompleteFoldings)
	fmt.Printf("Found %d solutions after CollapseTree\n", o.Counts.AfterCollapseTree)
	fmt.Printf("Found %d solutions after LiftCommon\n", o.Counts.AfterLiftCommon)

	fmt.Printf("Folding tree -> folding arrays conversion produced %d foldings\n", o.Counts.SafeCompletePairArrays)
	if o.Sanity.SafeComplete != "" {
		log.Print(o.Sanity.SafeComplete)
	}

	fmt.Println(o.SafeCompleteFoldingTree)
	if o.Sanity.Safety != "" {
		log.Print(o.Sanity.Safety)
	}
	if o.ReferenceFolding.Pairing != nil {
		fmt.Printf("Reference folding: (%d pairs)\n", o.ReferenceFolding.PairCount)
		fmt.Print(format.FoldingWithSafety(seq, o.ReferenceFolding.Pairing, o.Safety.SafeBase))
	}
	fmt.Println("Example folding with maximal pairing:")
	fmt.Print(format.FoldingWithSafety(seq, o.AllFoldings[0].Pairing, o.Safety.SafeBase))
	fmt.Printf("Safe bases %d/%d (%f %%)\n",
		o.Counts.SafeBases, o.Counts.SequenceBases, float64(o.Counts.SafeBases*100)/float64(o.Counts.SequenceBases))
	fmt.Print(o.ReferencePosition)

	return o
}

func foldingStats(seqs map[string]*base.Sequence, out chan OutputEntry) {
	fmt.Println("# Name NumBases FoldingPairs NumZuker  NumAll NumSafeBases TimeNussinov TimeWuchty TimeSafeComplete TimePairArrays TimeTrivSafety TimeSCSafety")

	for name := range seqs {
		seq := seqs[name]
		o := foldSequence(seq)

		if o.Sanity.Zuker != "" {
			log.Print(o.Sanity.Zuker)
		}
		if o.Sanity.Wuchty != "" {
			log.Print(o.Sanity.Wuchty)
		}
		if o.Sanity.SafeComplete != "" {
			log.Print(o.Sanity.SafeComplete)
		}
		if o.Sanity.Safety != "" {
			log.Print(o.Sanity.Safety)
		}

		out <- o
		fmt.Printf("%s %8d %12d %8d %7d %12d %12.6f %10.6f %16.6f %14.6f %14.6f %12.6f\n",
			o.Name, o.Counts.SequenceBases, o.Counts.OptimalPairs,
			o.Counts.ZukerFoldings, o.Counts.WuchtyFoldings, o.Counts.SafeBases,
			o.Timing.ZukerSeconds, o.Timing.WuchtySeconds, o.Timing.SafeCompleteSeconds, o.Timing.PairArraysSeconds,
			o.Timing.TrivialSafetySeconds, o.Timing.SafeCompleteSafetySeconds)
	}
	close(out)
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
