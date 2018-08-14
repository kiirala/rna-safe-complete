package fasta

import "testing"

func TestSplitLines(t *testing.T) {
	var tests = []struct {
		in       string
		count    int
		expected string
	}{
		{"abcdefg", 1, "a\nb\nc\nd\ne\nf\ng\n"},
		{"abcdefg", 2, "ab\ncd\nef\ng\n"},
		{"abcdef", 2, "ab\ncd\nef\n"},
		{"abcdefg", 80, "abcdefg\n"},
	}
	for _, tt := range tests {
		if actual := splitLines(tt.in, tt.count); actual != tt.expected {
			t.Errorf("splitlines(%q, %d): expected %q, got %q", tt.in, tt.count, tt.expected, actual)
		}
	}
}
