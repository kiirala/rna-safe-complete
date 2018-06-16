package fasta

import "testing"
import "strings"
import "bufio"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"

func TestReadComment(t *testing.T) {
	var tests = []struct {
		text     string
		expected string
	}{
		{">", ""},
		{";", ""},
		{">Foo", "Foo"},
		{";Foo", "Foo"},
		{"> Foo\nACGUT\n", "Foo"},
		{">sys|1234| foo [bar]\nACCGA", "sys|1234| foo [bar]"},
	}
	for _, tt := range tests {
		actual, err := readComment(bufio.NewReader(strings.NewReader(tt.text)))
		if err != nil {
			t.Errorf("readComment(%q): received error: %v", tt.text, err)
		}
		if !reflect.DeepEqual(tt.expected, actual) {
			t.Errorf("readComment(%q): expected %v, actual %v", tt.text, tt.expected, actual)
		}
	}
}

func TestReadCommentNoComment(t *testing.T) {
	input := "ACCGA\n>Comment here\nACCGA\n"
	actual, err := readComment(bufio.NewReader(strings.NewReader(input)))
	if actual != "" || err == nil || !strings.Contains(err.Error(), "no comment") {
		t.Errorf("readComment(%q): expected (\"\", error of no comment), actual (%v, %v)", input, actual, err)
	}
}

func TestReadBases(t *testing.T) {
	var tests = []struct {
		bases    string
		expected []base.Base
	}{
		{"", nil},
		{"\n", nil},
		{"\t \n", nil},
		{"ACGUT", []base.Base{base.A, base.C, base.G, base.U, base.U}},
		{"acgut", []base.Base{base.A, base.C, base.G, base.U, base.U}},
		{"A\nC\nG\nU\n", []base.Base{base.A, base.C, base.G, base.U}},
		{"ABCD", []base.Base{base.A, base.Other, base.C, base.Other}},
		// Comments inside a sequence are ignored
		{"AA\n;CC\n>UU\nGG", []base.Base{base.A, base.A, base.G, base.G}},
		// An asterisk at end of line  ends a sequence
		{"UUAA*\nGG\n", []base.Base{base.U, base.U, base.A, base.A}},
		{"UU\nAA\n*\nGG\n", []base.Base{base.U, base.U, base.A, base.A}},
		{"UU\nAA*GG", []base.Base{base.U, base.U, base.A, base.A, base.Other, base.G, base.G}},
		// An empty line ends a sequence
		{"AA\n\nGG\n", []base.Base{base.A, base.A}},
	}
	for _, tt := range tests {
		actual, err := readBases(bufio.NewReader(strings.NewReader(tt.bases)))
		if err != nil {
			t.Errorf("readBases(%q): received error: %v", tt.bases, err)
		}
		if !reflect.DeepEqual(tt.expected, actual) {
			t.Errorf("readBases(%q): expected %#v, actual %#v", tt.bases, tt.expected, actual)
		}
	}
}
