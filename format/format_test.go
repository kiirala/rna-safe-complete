package format /* import "keltainen.duckdns.org/rnafolding/format" */

import "testing"

import "keltainen.duckdns.org/rnafolding/folding"

func TestDotBracket(t *testing.T) {
	var tests = []struct {
		in       folding.FoldingPairs
		expected string
	}{
		{
			folding.FoldingPairs{},
			"",
		},
		{
			folding.FoldingPairs{1, 0},
			"()",
		},
		{
			folding.FoldingPairs{-1, -1, 3, 2},
			"..()",
		},
		{
			folding.FoldingPairs{3, 2, 1, 0},
			"(())",
		},
	}
	for _, tt := range tests {
		if actual := DotBracket(tt.in); actual != tt.expected {
			t.Errorf("DotBracket(%v): expected %v, actual %v", tt.in, tt.expected, actual)
		}
	}
}
