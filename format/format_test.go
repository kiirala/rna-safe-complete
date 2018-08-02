package format

import "testing"

import "keltainen.duckdns.org/rnafolding/types"

func TestDotBracket(t *testing.T) {
	var tests = []struct {
		in       types.FoldingPairs
		expected string
	}{
		{
			types.FoldingPairs{},
			"",
		},
		{
			types.FoldingPairs{1, 0},
			"()",
		},
		{
			types.FoldingPairs{-1, -1, 3, 2},
			"..()",
		},
		{
			types.FoldingPairs{3, 2, 1, 0},
			"(())",
		},
	}
	for _, tt := range tests {
		if actual := DotBracket(tt.in); actual != tt.expected {
			t.Errorf("DotBracket(%v): expected %v, actual %v", tt.in, tt.expected, actual)
		}
	}
}
