package main /* import "keltainen.duckdns.org/rnafolding/dbtofasta" */

import "flag"
import "log"
import "os"
import "path"

import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/trnadb"

var (
	db     = flag.String("db", "", "Location of the tRNA database file")
	outdir = flag.String("outdir", "", "Write FASTA files to given directory")
)

func main() {
	flag.Parse()

	if *db == "" || *outdir == "" {
		log.Fatalf("Both --db and --outdir parameters have to be set.")
	}

	f, err := os.Open(*db)
	if err != nil {
		log.Fatalf("Unable to open tRNA database \"%s\": %v", *db, err)
	}
	seqs, err := trnadb.ReadSequences(f)
	if err != nil {
		log.Fatalf("Error reading tRNA database file \"%s\": %v", *db, err)
	}

	for _, seq := range seqs {
		ofname := path.Join(*outdir, seq.Name+".fasta")
		of, err := os.Create(ofname)
		if err != nil {
			log.Fatalf("Error opening \"%s\" for writing: %v", ofname, err)
		}
		defer of.Close()
		err = fasta.WriteSequenceFolding(of, seq)
		if err != nil {
			log.Fatalf("Error writing FASTA format to \"%s\": %v", ofname, err)
		}
	}
}
