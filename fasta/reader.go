package fasta

import "bufio"
import "fmt"
import "io"
import "log"
import "strings"
import "unicode"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

func ReadSequence(r io.Reader) (*base.Sequence, error) {
	br := bufio.NewReader(r)
	seq, err := readSingleSequence(br)
	if err == io.EOF {
		err = nil
	}
	return seq, err
}

func ReadMultiSequence(r io.Reader) (map[string]*base.Sequence, error) {
	br := bufio.NewReader(r)
	ret := make(map[string]*base.Sequence)
	for {
		seq, err := readSingleSequence(br)
		if seq != nil {
			ret[seq.Name] = seq
		}
		if err == io.EOF {
			break
		} else if err != nil {
			return ret, err
		}
	}
	return ret, nil
}

func readSingleSequence(br *bufio.Reader) (*base.Sequence, error) {
	read, err := br.ReadString('\n')
	read = strings.TrimSpace(read)
	if err == io.EOF {
		if len(read) == 0 {
			return nil, err
		}
	} else if err != nil {
		return nil, err
	}

	seq := &base.Sequence{}

	var comm []string
	for {
		if len(read) == 0 {
			break
		}
		if read[0] == ';' || read[0] == '>' || read[0] == '#' {
			comm = append(comm, strings.TrimSpace(read[1:]))
		} else {
			break
		}
		read, err = nextLine(br)
		if err != nil {
			return seq, err
		}
		// STRAND files have an empty line between comment and sequence
		if len(read) == 0 {
			read, err = nextLine(br)
			if err != nil {
				return seq, err
			}
			break
		}
	}
	seq.Comment = strings.Join(comm, " / ")
	seq.Comment = strings.TrimPrefix(seq.Comment, "File ")
	seq.Name = strings.SplitN(seq.Comment, " ", 2)[0]

	for {
		if len(read) == 0 {
			break
		}
		if read[0] == '.' || read[0] == '(' || read[0] == '[' {
			break
		}
		if read[0] != ';' && read[0] != '>' {
			for _, code := range read {
				seq.Bases = append(seq.Bases, base.FromCode(string(code)))
			}
			if read[len(read)-1] == '*' {
				seq.Bases = seq.Bases[:len(seq.Bases)-1]
				break
			}
		}
		read, err = nextLine(br)
		if err != nil {
			return seq, err
		}
	}

	if len(read) > 0 && (read[0] == '.' || read[0] == '(' || read[0] == '[') {
		var dotbracket string
		for {
			if len(read) == 0 {
				break
			}
			// Dot-bracket may contain a '>' -- sometimes as first character of a line
			if read[0] != ';' && read[0] != '#' {
				dotbracket += read
			}
			read, err = nextLine(br)
			if err != nil {
				return seq, err
			}
		}
		seq.ReferenceFolding = DbToFold(dotbracket)
	}

	if len(seq.ReferenceFolding) > 0 && len(seq.Bases) != len(seq.ReferenceFolding) {
		log.Printf("Sequence and its folding are different lengths. Bases: %v, folding: %v", seq.Bases, seq.ReferenceFolding)
	}
	return seq, nil
}

func nextLine(r *bufio.Reader) (string, error) {
	read, err := r.ReadString('\n')
	read = strings.TrimSpace(read)
	if err == io.EOF {
		return read, nil
	}
	return read, err
}

type Ref struct {
	pos int
	typ rune
}

func toEndChar(c rune) rune {
	if c == '(' {
		return ')'
	} else if c == '{' {
		return '}'
	} else if c == '[' {
		return ']'
	} else if c == '<' {
		return '>'
	} else if unicode.IsUpper(c) {
		return unicode.ToLower(c)
	}
	return ' '
}

func getWithType(stack []Ref, c rune) (int, []Ref, error) {
	for i := len(stack) - 1; i >= 0; i-- {
		if stack[i].typ == c {
			if i == len(stack)-1 {
				return stack[i].pos, stack[:len(stack)-1], nil
			} else if i == 0 {
				return stack[i].pos, stack[1:], nil
			} else {
				// Grab the .pos before the stack is modified
				r := stack[i].pos
				return r, append(stack[:i], stack[i+1:]...), nil
			}
		}
	}
	return 0, stack, fmt.Errorf("Didn't find matching pair for '%c' (%v) in %v", c, c, stack)
}

func DbToFold(s string) folding.FoldingPairs {
	if s == "" {
		return nil
	}
	var stack []Ref
	fold := make(folding.FoldingPairs, len(s))
	for i, c := range s {
		if c == '.' {
			fold[i] = -1
		} else if ec := toEndChar(c); ec != ' ' {
			stack = append(stack, Ref{i, ec})
		} else {
			var j int
			var err error
			j, stack, err = getWithType(stack, c)
			if err != nil {
				log.Printf("Error at char %d: %v", i, err)
			}
			fold[i] = j
			fold[j] = i
		}
	}
	if len(stack) > 0 {
		log.Printf("Stack not empty at end of parsing folding. Stack: %v", stack)
	}
	return fold
}
