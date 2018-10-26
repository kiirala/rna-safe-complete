package fasta /* import "keltainen.duckdns.org/rnafolding/fasta" */

import "io"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/format"

func WriteSequenceFolding(of io.Writer, seq *base.Sequence) error {
	_, err := io.WriteString(of, ">"+seq.Comment+"\n")
	if err != nil {
		return err
	}

	bases := seq.BasesString()
	dotbr := format.DotBracket(seq.ReferenceFolding)

	bss := splitLines(bases, 80)
	_, err = io.WriteString(of, bss)
	if err != nil {
		return err
	}

	dbss := splitLines(dotbr, 80)
	_, err = io.WriteString(of, dbss)
	if err != nil {
		return err
	}

	return nil
}

func splitLines(i string, c int) string {
	out := ""
	for len(i) > c {
		out += i[:c] + "\n"
		i = i[c:]
	}
	if len(i) > 0 {
		out += i + "\n"
	}
	return out
}
