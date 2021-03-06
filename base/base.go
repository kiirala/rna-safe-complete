package base /* import "keltainen.duckdns.org/rnafolding/base" */

import "strings"

import "keltainen.duckdns.org/rnafolding/folding"

type Base int

const (
	Unknown Base = iota
	A
	C
	G
	U
	Other
)

type Sequence struct {
	Name              string
	Comment           string
	Bases             []Base
	ReferenceFolding  folding.FoldingPairs
	ReferencePosition []int
}

func (b Base) String() string {
	switch b {
	case Unknown:
		return "Unknown"
	case A:
		return "Adenine"
	case C:
		return "Cytocine"
	case G:
		return "Guanine"
	case U:
		return "Uracil"
	case Other:
		return "Other"
	}
	return "Bad"
}

func FromCode(char string) Base {
	switch char {
	case "A", "a":
		return A
	case "C", "c":
		return C
	case "G", "g":
		return G
	case "U", "u", "T", "t":
		return U
	}
	return Other
}

func (b Base) ToCode() string {
	switch b {
	case Unknown:
		return "-"
	case A:
		return "A"
	case C:
		return "C"
	case G:
		return "G"
	case U:
		return "U"
	case Other:
		return "N"
	}
	return "."
}

func (b Base) CodeAndSafety(safe bool) string {
	c := b.ToCode()
	if !safe {
		return strings.ToLower(c)
	}
	return c
}

func (b Base) CanPair(o Base) bool {
	switch b {
	case A:
		if o == U {
			return true
		}
	case U:
		if o == A || o == G {
			return true
		}
	case G:
		if o == C || o == U {
			return true
		}
	case C:
		if o == G {
			return true
		}
	}
	return false
}

func (s Sequence) CanPair(i, j, minloop int) bool {
	return s.Bases[i].CanPair(s.Bases[j]) && (j-i > minloop || j-i < -minloop)
}

func SequenceFromString(s string) *Sequence {
	bases := make([]Base, len(s))
	for i, c := range s {
		bases[i] = FromCode(string(c))
	}
	return &Sequence{Bases: bases}
}

func (s Sequence) BasesString() string {
	out := make([]string, len(s.Bases))
	for i, b := range s.Bases {
		out[i] = b.ToCode()
	}
	return strings.Join(out, "")
}

func (s Sequence) BasesSafetyString(safety []bool) string {
	out := make([]string, len(s.Bases))
	for i, b := range s.Bases {
		out[i] = b.CodeAndSafety(safety[i])
	}
	return strings.Join(out, "")
}
