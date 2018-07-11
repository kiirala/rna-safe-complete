// Package trnadb reads single tRNA sequences from a tRNA database.
// The database can be found at ftp://ftp.ebi.ac.uk/pub/databases/trna/
// For details of the database, see:
// Sprinzl, Steegborn, Hübel, Steinberg - Compilation of tRNA sequences and sequences of tRNA genes
// Nucleic Acids Research, 1996, vol. 24, no. 1, pp. 68 -- 72
package trnadb

import "bufio"
import "io"
import "strings"
import "unicode"

import "keltainen.duckdns.org/rnafolding/base"

// Map from modified bases to non-modified ones.
// This mapping copied from:
// Paul G. Higgs - Thermodynamic Properties of Transfer RNA: A Computational Study
// J. Chem. Soc. Faraday Trans., 1995, vol. 91, no. 16, pp. 2531--2540
var baseMapping = map[rune]base.Base{
	'H': base.A, // Unknown modified adenosine
	'^': base.A, // 2'-O-ribosyladenosine (phosphat)

	'<': base.C, // Unknown modified cytidine
	'B': base.C, // 2'-O-methylcytidine
	'M': base.C, // N4-acetylcytidine
	'?': base.C, // 5-Methylcytidine

	';': base.G, // Unknown modified guanosine
	'L': base.G, // N2-methylguanosine
	'#': base.G, // 2'-O-methylguanosine
	'R': base.G, // N2,N2-dimethylguanosine

	'N': base.U, // Unknown modified uridine
	'J': base.U, // 2'-O-methyluridine
	'P': base.U, // Pseudouridine
	']': base.U, // 1-Methylpseudouridine
	'Z': base.U, // 2'-O-methylpseudouridine
}

func parseLine(s string) (string, []base.Base) {
	comment := s[:36]
	var bases []base.Base
	for _, c := range s[36:] {
		if unicode.IsSpace(c) || c == '-' {
			continue
		}
		if b, ok := baseMapping[c]; ok {
			bases = append(bases, b)
		} else {
			bases = append(bases, base.FromCode(string(c)))
		}
	}
	return comment, bases
}

func ReadSequences(r io.Reader) (map[string]*base.Sequence, error) {
	br := bufio.NewReader(r)
	seqs := map[string]*base.Sequence{}
	for {
		c, err := br.ReadString('\n')
		if err != nil && err != io.EOF {
			return nil, err
		}

		// All RNA sequences in the DB file start with an 'R'
		if len(c) > 0 && c[0] == 'R' {
			comment, bases := parseLine(c)
			name := strings.SplitN(comment, " ", 2)[0]
			seqs[name] = &base.Sequence{
				Comment: comment,
				Bases:   bases,
			}
		}
		if err == io.EOF {
			break
		}
	}
	return seqs, nil
}
