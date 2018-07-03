package types

import "testing"
import "reflect"

func TestCollapseTreeAlternatives(t *testing.T) {
	f := &Folding{
		Pairs: []Pair{Pair{I: 0, J: 13}},
		JoinPrefix: &Folding{
			Pairs: []Pair{Pair{I: 1, J: 2}},
		},
		JoinSuffix: &Folding{
			JoinPrefix: &Folding{
				Pairs: []Pair{Pair{I: 3, J: 7}},
				Branches: []*Folding{
					&Folding{
						Free:  []int{4},
						Pairs: []Pair{Pair{I: 5, J: 6}},
					},
					&Folding{
						Free:  []int{6},
						Pairs: []Pair{Pair{I: 4, J: 5}},
					},
				},
			},
			JoinSuffix: &Folding{
				Pairs: []Pair{Pair{I: 8, J: 12}},
				Branches: []*Folding{
					&Folding{
						Free:  []int{9},
						Pairs: []Pair{Pair{I: 10, J: 11}},
					},
					&Folding{
						Free:  []int{11},
						Pairs: []Pair{Pair{I: 9, J: 10}},
					},
				},
			},
		},
	}

	expected := &Folding{
		Pairs: []Pair{Pair{I: 0, J: 13}, Pair{I: 1, J: 2}, Pair{I: 3, J: 7}, Pair{I: 8, J: 12}},
		JoinPrefix: &Folding{
			Branches: []*Folding{
				&Folding{
					Free:  []int{4},
					Pairs: []Pair{Pair{I: 5, J: 6}},
				},
				&Folding{
					Free:  []int{6},
					Pairs: []Pair{Pair{I: 4, J: 5}},
				},
			},
		},
		JoinSuffix: &Folding{
			Branches: []*Folding{
				&Folding{
					Free:  []int{9},
					Pairs: []Pair{Pair{I: 10, J: 11}},
				},
				&Folding{
					Free:  []int{11},
					Pairs: []Pair{Pair{I: 9, J: 10}},
				},
			},
		},
	}

	f.CollapseTree()
	if !reflect.DeepEqual(f, expected) {
		t.Errorf("CollapseTree(): got:\n%v\nexpected\n%v\n", f, expected)
	}
}
