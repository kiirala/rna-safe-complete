package fasta

import "testing"
import "strings"
import "bufio"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

func TestReadSequence(t *testing.T) {
	var tests = []struct {
		input   string
		comment string
		bases   []base.Base
		folding folding.FoldingPairs
	}{
		// Only comment
		{">", "", nil, nil},
		{";", "", nil, nil},
		{">Foo", "Foo", nil, nil},
		{";Foo", "Foo", nil, nil},
		{">sys|1234| foo [bar]\n", "sys|1234| foo [bar]", nil, nil},

		// Trivial comment + sequence
		{">Cmt\n", "Cmt", nil, nil},
		{">Cmt\n\n", "Cmt", nil, nil},
		{">Cmt\n\t \n", "Cmt", nil, nil},
		{">Cmt\nACGUT", "Cmt", []base.Base{base.A, base.C, base.G, base.U, base.U}, nil},
		{">Cmt\nacgut", "Cmt", []base.Base{base.A, base.C, base.G, base.U, base.U}, nil},
		{">Cmt\nA\nC\nG\nU\n", "Cmt", []base.Base{base.A, base.C, base.G, base.U}, nil},
		{">Cmt\nABCD", "Cmt", []base.Base{base.A, base.Other, base.C, base.Other}, nil},
		// Comments inside a sequence are ignored
		{">Cmt\nAA\n;CC\n>UU\nGG", "Cmt", []base.Base{base.A, base.A, base.G, base.G}, nil},
		// An asterisk at end of line  ends a sequence
		{">Cmt\nUUAA*\nGG\n", "Cmt", []base.Base{base.U, base.U, base.A, base.A}, nil},
		{">Cmt\nUU\nAA\n*\nGG\n", "Cmt", []base.Base{base.U, base.U, base.A, base.A}, nil},
		{">Cmt\nUU\nAA*GG", "Cmt", []base.Base{base.U, base.U, base.A, base.A, base.Other, base.G, base.G}, nil},
		// An empty line ends a sequence
		{">Cmt\nAA\n\nGG\n", "Cmt", []base.Base{base.A, base.A}, nil},

		// Comment, sequence and folding
		{
			">Cmt\nACGUACGU\n.((..)).\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{-1, 6, 5, -1, -1, 2, 1, -1},
		},
		{
			">Cmt\nACGUACGU\n((...)).\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{6, 5, -1, -1, -1, 1, 0, -1},
		},
		{
			">Cmt\nACGUACGU\n.((...))\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{-1, 7, 6, -1, -1, -1, 2, 1},
		},
		// Folding with pseudoloop
		{
			">Cmt\nACGUACGUACGU\n.(([[..)).]]\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{-1, 8, 7, 11, 10, -1, -1, 2, 1, -1, 4, 3},
		},
		{
			">Cmt\nACGUACGUACGU\n.(([[())).]]\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{-1, 8, 7, 11, 10, 6, 5, 2, 1, -1, 4, 3},
		},
		{
			">Cmt\nACGUACGUACGU\n(.([[..)).]]\n",
			"Cmt",
			[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
			folding.FoldingPairs{8, -1, 7, 11, 10, -1, -1, 2, 0, -1, 4, 3},
		},
	}
	for _, tt := range tests {
		actual, err := ReadSequence(bufio.NewReader(strings.NewReader(tt.input)))
		if err != nil {
			t.Errorf("ReadSequence(%q): received error: %v", tt.input, err)
		}
		if actual == nil {
			t.Errorf("ReadSequence(%q): returned nil", tt.input)
		} else {
			if tt.comment != actual.Comment {
				t.Errorf("ReadSequence(%q): expected comment %#v, actual %#v", tt.input, tt.comment, actual.Comment)
			}
			if !reflect.DeepEqual(tt.bases, actual.Bases) {
				t.Errorf("ReadSequence(%q): expected bases %#v, actual %#v", tt.input, tt.bases, actual.Bases)
			}
			if !reflect.DeepEqual(tt.folding, actual.ReferenceFolding) {
				t.Errorf("ReadSequence(%q): expected folding %#v, actual %#v", tt.input, tt.folding, actual.ReferenceFolding)
			}
		}
	}
}

func TestReadMultiSequence(t *testing.T) {
	type seqs []struct {
		comment string
		bases   []base.Base
		folding folding.FoldingPairs
	}
	var tests = []struct {
		input string
		seqs  seqs
	}{
		{
			">Cmt1\nACGUACGU\n.((..)).\n\n" +
				">Cmt2\nACGUACGU\n((...)).\n\n" +
				">Cmt3\nACGUACGU\n.((...))\n",
			seqs{
				{
					"Cmt1",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{-1, 6, 5, -1, -1, 2, 1, -1},
				},
				{
					"Cmt2",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{6, 5, -1, -1, -1, 1, 0, -1},
				},
				{
					"Cmt3",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{-1, 7, 6, -1, -1, -1, 2, 1},
				},
			},
		},

		{
			">Cmt4\nACGUACGUACGU\n.(([[..)).]]\n\n" +
				">Cmt5\nACGUACGUACGU\n.(([[())).]]\n\n" +
				">Cmt6\nACGUACGUACGU\n(.([[..)).]]\n",
			seqs{
				{
					"Cmt4",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{-1, 8, 7, 11, 10, -1, -1, 2, 1, -1, 4, 3},
				},
				{
					"Cmt5",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{-1, 8, 7, 11, 10, 6, 5, 2, 1, -1, 4, 3},
				},
				{
					"Cmt6",
					[]base.Base{base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U, base.A, base.C, base.G, base.U},
					folding.FoldingPairs{8, -1, 7, 11, 10, -1, -1, 2, 0, -1, 4, 3},
				},
			},
		},
	}

	for _, tt := range tests {
		actuals, err := ReadMultiSequence(bufio.NewReader(strings.NewReader(tt.input)))
		if err != nil {
			t.Errorf("ReadSequence(%q): received error: %v", tt.input, err)
		}
		if len(actuals) != len(tt.seqs) {
			t.Errorf("ReadSequence(%q): expected %d sequences, actual %d", tt.input, len(tt.seqs), len(actuals))
		}
		for _, seq := range tt.seqs {
			actual := actuals[seq.comment]
			if actual == nil {
				t.Errorf("ReadSequence(%q): returned nil for sequence %q", tt.input, seq.comment)
			} else {
				if seq.comment != actual.Comment {
					t.Errorf("ReadSequence(%q): expected comment %#v, actual %#v", tt.input, seq.comment, actual.Comment)
				}
				if !reflect.DeepEqual(seq.bases, actual.Bases) {
					t.Errorf("ReadSequence(%q): expected bases %#v, actual %#v", tt.input, seq.bases, actual.Bases)
				}
				if !reflect.DeepEqual(seq.folding, actual.ReferenceFolding) {
					t.Errorf("ReadSequence(%q): expected folding %#v, actual %#v", tt.input, seq.folding, actual.ReferenceFolding)
				}
			}
		}
	}
}
