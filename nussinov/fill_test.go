package nussinov

import "testing"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/format"

func TestFillArray(t *testing.T) {
	tests := []struct {
		bases    []base.Base
		expected [][]int
	}{
		{
			[]base.Base{base.A, base.U},
			[][]int{
				{0, 1},
				{0, 0},
			},
		},
		{
			[]base.Base{base.A, base.C},
			[][]int{
				{0, 0},
				{0, 0},
			},
		},
		{
			[]base.Base{base.G, base.G, base.G, base.A, base.A, base.A, base.U, base.C, base.C},
			[][]int{
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 2},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
			},
		},
	}

	for _, tt := range tests {
		if actual := FillArray(&base.Sequence{Bases: tt.bases}); !reflect.DeepEqual(tt.expected, actual) {
			t.Errorf("FillArray(%v):\nexpected:\n%v\nactual:\n%v", tt.bases, format.Matrix(tt.expected), format.Matrix(actual))
		}
	}
}

func TestFillComplementary(t *testing.T) {
	tests := []struct {
		bases    []base.Base
		v        [][]int
		expected [][]int
	}{
		{
			[]base.Base{base.A, base.U},
			[][]int{
				{0, 1},
				{0, 0},
			},
			[][]int{
				{0, 0},
				{0, 0},
			},
		},
		{
			[]base.Base{base.G, base.G, base.G, base.A, base.A, base.A, base.U, base.C, base.C},
			[][]int{
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 2},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 1, 1, 1},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
			},
			[][]int{
				{3, 2, 1, 1, 1, 0, 0, 0, 0},
				{0, 3, 2, 2, 2, 1, 1, 1, 0},
				{0, 0, 3, 3, 3, 2, 2, 1, 0},
				{0, 0, 0, 3, 3, 3, 2, 1, 0},
				{0, 0, 0, 0, 3, 3, 2, 1, 0},
				{0, 0, 0, 0, 0, 3, 2, 1, 0},
				{0, 0, 0, 0, 0, 0, 2, 1, 0},
				{0, 0, 0, 0, 0, 0, 0, 2, 1},
				{0, 0, 0, 0, 0, 0, 0, 0, 2},
			},
		},
	}

	for _, tt := range tests {
		if actual := FillComplementary(&base.Sequence{Bases: tt.bases}, tt.v); !reflect.DeepEqual(tt.expected, actual) {
			t.Errorf("FillComplementary(%v):\nexpected:\n%v\nactual:\n%v", tt.bases, format.Matrix(tt.expected), format.Matrix(actual))
		}
	}
}
