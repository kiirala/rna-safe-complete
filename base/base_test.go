package base

import "testing"

func TestString(t *testing.T) {
	var tests = []struct {
		enum     Base
		expected string
	}{
		{Unknown, "Unknown"},
		{A, "Adenine"},
		{C, "Cytocine"},
		{G, "Guanine"},
		{U, "Uracil"},
		{Other, "Other"},
		{42, "Bad"},
	}
	for _, tt := range tests {
		if actual := tt.enum.String(); actual != tt.expected {
			t.Errorf("%v.String(): expected %v, actual %v", tt.enum, tt.expected, actual)
		}
	}
}

func TestFromCode(t *testing.T) {
	var tests = []struct {
		code     string
		expected Base
	}{
		{"A", A},
		{"a", A},
		{"C", C},
		{"c", C},
		{"G", G},
		{"g", G},
		{"U", U},
		{"u", U},
		{"T", U},
		{"t", U},
		{"X", Other},
		{".", Other},
	}
	for _, tt := range tests {
		if actual := FromCode(tt.code); actual != tt.expected {
			t.Errorf("Base.FromCode(\"%v\"): expected %v, actual %v", tt.code, tt.expected, actual)
		}
	}
}

func TestToCode(t *testing.T) {
	var tests = []struct {
		enum     Base
		expected string
	}{
		{Unknown, "-"},
		{A, "A"},
		{C, "C"},
		{G, "G"},
		{U, "U"},
		{Other, "N"},
		{42, "."},
	}
	for _, tt := range tests {
		if actual := tt.enum.ToCode(); actual != tt.expected {
			t.Errorf("%v.ToCode(): expected %v,actual %v", tt.enum, tt.expected, actual)
		}
	}
}
