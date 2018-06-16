package fasta

import "io"
import "bufio"
import "errors"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"

func ReadSequence(r io.Reader) (*base.Sequence, error) {
	br := bufio.NewReader(r)
	comment, err := readComment(br)
	if err != nil {
		return nil, err
	}
	bases, err := readBases(br)
	return &base.Sequence{
		Comment: comment,
		Bases:   bases,
	}, err
}

func readComment(r *bufio.Reader) (string, error) {
	c, err := r.ReadString('\n')
	if err != nil && err != io.EOF {
		return "", err
	}
	if c[0] != '>' && c[0] != ';' {
		return "", errors.New("fasta.ReadSequence: no comment at beginning of file")
	}
	return strings.TrimSpace(c[1:]), nil
}

func readBases(r *bufio.Reader) ([]base.Base, error) {
	var bases []base.Base
	for {
		read, err := r.ReadString('\n')
		read = strings.TrimSpace(read)
		if len(read) == 0 {
			if err != nil && err != io.EOF {
				return bases, err
			}
			break
		}
		if read[0] == ';' || read[0] == '>' {
			continue
		}
		for _, code := range read {
			bases = append(bases, base.FromCode(string(code)))
		}
		if read[len(read)-1] == '*' {
			bases = bases[:len(bases)-1]
			break
		}
		if err == io.EOF {
			break
		} else if err != nil {
			return bases, err
		}
	}
	return bases, nil
}
