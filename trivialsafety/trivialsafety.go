package main

import "bufio"
import "io"
import "encoding/json"
import "log"
import "os"
import "strings"

import "keltainen.duckdns.org/rnafolding/fasta"
import "keltainen.duckdns.org/rnafolding/folding"

func main() {
	r := bufio.NewReader(os.Stdin)

	name, seq := readHeader(r)

	folds := make(chan folding.FoldingPairs, 1)
	go getFoldings(r, folds)

	numFolds := 0
	pairs := make([][]int, len(seq))
	free := make([]int, len(seq))
	for i := 0; i < len(seq); i++ {
		pairs[i] = make([]int, len(seq))
	}

	for fold := range folds {
		if len(fold) != len(seq) {
			log.Printf("Sequence length is %d, but folding length %d", len(seq), len(fold))
		}
		numFolds++
		for i, j := range fold {
			if j < 0 {
				free[i]++
			} else if i < j {
				pairs[i][j]++
			}
		}
	}

	j, err := json.Marshal(
		struct {
			Name     string
			Bases    int
			NumFolds int
			Pairs    [][]int
			Free     []int
		}{
			name, len(seq), numFolds, pairs, free,
		})
	if err != nil {
		log.Printf("Failed to write results as JSON: %v", err)
	}
	os.Stdout.Write(j)
	os.Stdout.Write([]byte{'\n'})
}

func readHeader(r *bufio.Reader) (string, string) {
	cmt, err := r.ReadString('\n')
	if err != nil {
		log.Printf("Failed to read comment: %v", err)
		return "", ""
	}
	if len(cmt) < 4 || cmt[0] != '>' {
		log.Printf("Input doesn't begin with comment, got %q instead", cmt)
	}
	name := strings.SplitN(cmt, " ", 3)[0][1:]

	bb, err := r.ReadString('\n')
	if err != nil {
		log.Printf("Failed to read sequence: %v", err)
		return cmt, ""
	}
	if bb[0] != 'A' && bb[0] != 'C' && bb[0] != 'G' && bb[0] != 'U' && bb[0] != 'N' {
		log.Printf("Line 2 of input doesn't look like sequence, got %q instead", bb)
	}
	seq := strings.SplitN(bb, " ", 2)[0]

	return name, seq
}

func getFoldings(r *bufio.Reader, ret chan folding.FoldingPairs) {
	defer close(ret)

	for {
		fs, err := r.ReadString('\n')
		if err != nil && err != io.EOF {
			log.Printf("Failed to read folding: %v", err)
		}
		if len(fs) > 0 {
			db := strings.SplitN(fs, " ", 2)[0]
			fold := fasta.DbToFold(db)
			ret <- fold
		}

		if err == io.EOF {
			return
		}
	}
}
