package trnadb /* import "keltainen.duckdns.org/rnafolding/trnadb" */

import "testing"

func TestBaseNames(t *testing.T) {
	var tests = []struct {
		name     string
		actual   baseName
		expected int
	}{
		{"b0", b0, 0},
		{"b10", b10, 10},
		{"b17", b17, 17},
		{"b17a", b17a, 18},
		{"b26", b26, 29},
		{"b45", b45, 48},
		{"b49", b49, 71},
		{"b66", b66, 88},
		{"b76", b76, 98},
	}
	for _, tt := range tests {
		if tt.expected != int(tt.actual) {
			t.Errorf("baseName %s=%d, expected %d", tt.name, tt.actual, tt.expected)
		}
	}
}
