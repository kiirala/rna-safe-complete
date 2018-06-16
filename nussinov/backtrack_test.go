package nussinov

import "testing"
import "reflect"

import "keltainen.duckdns.org/rnafolding/base"

func TestBacktrack(t *testing.T) {
	tests := []struct {
		sequence      *base.Sequence
		array         [][]int
		expectedCount int
		expectedPairs []int
	}{
		{
			base.SequenceFromString("CG"),
			[][]int{
				{0, 1},
				{0, 0},
			},
			1,
			[]int{1, 0},
		},
		{
			base.SequenceFromString("AC"),
			[][]int{
				{0, 0},
				{0, 0},
			},
			0,
			[]int{-1, -1},
		},
		{
			base.SequenceFromString("GGGAAAUCC"),
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
			3,
			[]int{8, 7, 6, -1, -1, -1, 2, 1, 0},
		},
	}

	for _, tt := range tests {
		actualCount, actualPairs := Backtrack(tt.sequence, tt.array)
		if tt.expectedCount != actualCount || !reflect.DeepEqual(tt.expectedPairs, actualPairs) {
			t.Errorf("Backtrack(%v): expected: %d pairs %v actual: %d pairs %v",
				tt.array, tt.expectedCount, tt.expectedPairs, actualCount, actualPairs)
		}
	}
}

func TestJoinedBacktrack(t *testing.T) {
	arr := [][]int{
		{0, 0, 0, 0, 0, 0, 1, 2, 3},
		{0, 0, 0, 0, 0, 0, 1, 2, 3},
		{0, 0, 0, 0, 0, 0, 1, 2, 2},
		{0, 0, 0, 0, 0, 0, 1, 1, 1},
		{0, 0, 0, 0, 0, 0, 1, 1, 1},
		{0, 0, 0, 0, 0, 0, 1, 1, 1},
		{0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0},
	}
	cmpl := [][]int{
		{3, 2, 1, 1, 1, 0, 0, 0, 0},
		{0, 3, 2, 2, 2, 1, 1, 1, 0},
		{0, 0, 3, 3, 3, 2, 2, 1, 0},
		{0, 0, 0, 3, 3, 3, 2, 1, 0},
		{0, 0, 0, 0, 3, 3, 2, 1, 0},
		{0, 0, 0, 0, 0, 3, 2, 1, 0},
		{0, 0, 0, 0, 0, 0, 2, 1, 0},
		{0, 0, 0, 0, 0, 0, 0, 2, 1},
		{0, 0, 0, 0, 0, 0, 0, 0, 2},
	}
	seq := base.SequenceFromString("GGGAAAUCC")
	tests := []struct {
		i        int
		j        int
		expected []int
	}{
		{0, 8, []int{8, 7, 6, -1, -1, -1, 2, 1, 0}},
		{1, 8, []int{-1, 8, 7, 6, -1, -1, 3, 2, 1}},
		{4, 6, []int{8, 7, -1, -1, 6, -1, 4, 1, 0}},
	}
	for _, tt := range tests {
		if _, actual := JoinedBacktrack(seq, arr, cmpl, tt.i, tt.j); !reflect.DeepEqual(tt.expected, actual) {
			t.Errorf("JoinedBacktrack(seq, arr, %v, %v): expected %v, actual %v", tt.i, tt.j, tt.expected, actual)
		}
	}
}
