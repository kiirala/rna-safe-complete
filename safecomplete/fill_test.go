package safecomplete /* import "keltainen.duckdns.org/rnafolding/safecomplete" */

import "testing"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/format"

func TestFillArray(t *testing.T) {
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
				{0, 1},
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
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 3},
				{0, 0, 0, 0, 0, 0, 1, 2, 2},
				{0, 0, 0, 0, 0, 0, 1, 0, 0},
				{0, 0, 0, 0, 0, 0, 1, 0, 0},
				{0, 0, 0, 0, 0, 0, 1, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0},
			},
		},
	}

	for _, tt := range tests {
		p := &Predictor{Seq: &base.Sequence{Bases: tt.bases}, V: tt.v}
		p.FillArray()
		if !reflect.DeepEqual(tt.expected, p.W) {
			t.Errorf("FillArray(%v):\nexpected:\n%v\nactual:\n%v", tt.bases, format.Matrix(tt.expected), format.Matrix(p.W))
		}
	}
}
