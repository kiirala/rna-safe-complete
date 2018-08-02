// Package trnadb reads single tRNA sequences from a tRNA database.
// The database can be found at ftp://ftp.ebi.ac.uk/pub/databases/trna/
// For details of the database, see:
// Sprinzl, Steegborn, Hübel, Steinberg - Compilation of tRNA sequences and sequences of tRNA genes
// Nucleic Acids Research, 1996, vol. 24, no. 1, pp. 68 -- 72
package trnadb

import "bufio"
import "io"
import "log"
import "strings"
import "unicode"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

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

type baseName int

const (
	b0 = iota // Acceptor stem
	b1
	b2
	b3
	b4
	b5
	b6
	b7

	b8
	b9

	b10 // D-domain
	b11
	b12
	b13
	b14
	b15
	b16
	b17
	b17a
	b18
	b19
	b20
	b20a
	b20b
	b21
	b22
	b23
	b24
	b25

	b26 // anticodon domain
	b27
	b28
	b29
	b30
	b31
	b32
	b33
	b34
	b35
	b36
	b37
	b38
	b39
	b40
	b41
	b42
	b43
	b44

	b45 // variable region
	e11 // extra loop
	e12
	e13
	e14
	e15
	e16
	e17
	e1
	e2
	e3
	e4
	e5
	e27
	e26
	e25
	e24
	e23
	e22
	e21
	b46
	b47
	b48

	b49 // T-domain
	b50
	b51
	b52
	b53
	b54
	b55
	b56
	b57
	b58
	b59
	b60
	b61
	b62
	b63
	b64
	b65

	b66 // acceptor stem
	b67
	b68
	b69
	b70
	b71
	b72

	b73
	b74
	b75
	b76
)

const (
	basesStartColumn = 40
)

func parseLine(s string) (string, []base.Base, []int) {
	sepIdx := strings.LastIndex(s, "\t")
	//log.Printf("%q %q", s[:sepIdx], s[sepIdx+1:])
	//log.Printf("%q", s[sepIdx+1:])
	comment := s[:sepIdx]
	var bases []base.Base
	var pos []int
	for i, c := range s[sepIdx+1:] {
		if unicode.IsSpace(c) || c == '-' {
			continue
		}
		if b, ok := baseMapping[c]; ok {
			bases = append(bases, b)
		} else {
			bases = append(bases, base.FromCode(string(c)))
		}
		pos = append(pos, i)
	}
	return comment, bases, pos
}

func lookup(pos []int, n baseName) int {
	for i, v := range pos {
		if v == int(n) {
			return i
		}
	}
	log.Printf("lookup: position %d not found in %v", n, pos)
	return 0
}

func joinranges(seq []base.Base, f folding.FoldingPairs, s string, pos []int, a5, a3, b5, b3 baseName) {
	a := a5
	b := b3
	for a <= a3 && b >= b5 {
		if s[a] == '=' || s[a] == '*' {
			if s[b] == s[a] {
				apos := lookup(pos, a)
				bpos := lookup(pos, b)
				if s[a] == '=' {
					if !((seq[apos] == base.A && seq[bpos] == base.U) ||
						(seq[apos] == base.U && seq[bpos] == base.A) ||
						(seq[apos] == base.G && seq[bpos] == base.C) ||
						(seq[apos] == base.C && seq[bpos] == base.G)) {
						// Some tRNA sequences in the DB actually do have non-traditional pairs.
						//log.Printf("joinrange: expected Watson-Crick pair, received %d (%v) to %d (%v) (%d-%d, %d-%d)",
						//	a, seq[apos], b, seq[bpos], a5, a3, b5, b3)
					}
				} else {
					if !((seq[apos] == base.U && seq[bpos] == base.G) ||
						(seq[apos] == base.G && seq[bpos] == base.U)) {
						log.Printf("joinrange: expected wobble pair, received %d (%v) to %d (%v) (%d-%d, %d-%d)",
							a, seq[apos], b, seq[bpos], a5, a3, b5, b3)
					}
				}
				f[apos] = bpos
				f[bpos] = apos
				a++
				b--
			} else if !(s[b] == '=' || s[b] == '*') {
				b--
			} else {
				log.Printf("joinrange: pairing mismatch: s[%d]='%c', s[%d]='%c'", a, s[a], b, s[b])
				a++
				b--
			}
		} else {
			a++
		}
	}
}

func parseFold(seq []base.Base, raws string, positions []int) folding.FoldingPairs {
	s := raws[basesStartColumn:]
	f := folding.NewFoldingPairs(len(positions))
	//log.Print(seq)
	//log.Printf("%q", s)
	joinranges(seq, f, s, positions, b1, b7, b66, b72)
	joinranges(seq, f, s, positions, b10, b13, b22, b25)
	joinranges(seq, f, s, positions, b26, b31, b39, b44)
	joinranges(seq, f, s, positions, e11, e17, e27, e21)
	joinranges(seq, f, s, positions, b45, b45, b46, b48)
	joinranges(seq, f, s, positions, b49, b53, b61, b65)
	return f
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
			//log.Printf("%q", c)
			comment, bases, positions := parseLine(c)
			name := strings.SplitN(comment, " ", 2)[0]
			foldLine, err := br.ReadString('\n')
			if err != nil {
				return nil, err
			}
			fold := parseFold(bases, foldLine, positions)
			seqs[name] = &base.Sequence{
				Name:             name,
				Comment:          comment,
				Bases:            bases,
				ReferenceFolding: fold,
			}
		}
		if err == io.EOF {
			break
		}
	}
	return seqs, nil
}
