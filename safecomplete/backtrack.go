package safecomplete

import "keltainen.duckdns.org/rnafolding/base"

type Pair struct {
	I int
	J int
}

type Folding struct {
	Pairs      []Pair
	Free       []int
	JoinPrefix *Folding
	JoinSuffix *Folding
	Branches   []*Folding
}

func BacktrackAll(seq *base.Sequence, v, w [][]int) *Folding {
	folding := &Folding{}
	recursiveBacktrack(seq, v, w, 0, len(v)-1, folding)
	collapseTree(folding)
	return folding
}

func recursiveBacktrack(seq *base.Sequence, v, w [][]int, i, j int, folding *Folding) {
	if i >= j {
		if i == j {
			folding.Free = append(folding.Free, i)
		}
		return
	}
	numFound := 0
	if seq.CanPair(i, j) && v[i][j] == v[i+1][j-1]+1 {
		numFound++
	}
	for k := i; k < j; k++ {
		if v[i][j] == v[i][k]+w[k+1][j] && w[k+1][j] > 0 {
			numFound++
		}
	}
	if seq.CanPair(i, j) && v[i][j] == v[i+1][j-1]+1 {
		branch := folding
		if numFound > 1 {
			branch = &Folding{}
			folding.Branches = append(folding.Branches, branch)
		}
		branch.Pairs = append(branch.Pairs, Pair{I: i, J: j})
		recursiveBacktrack(seq, v, w, i+1, j-1, branch)
	}
	for k := i; k < j; k++ {
		if v[i][j] == v[i][k]+w[k+1][j] && w[k+1][j] > 0 {
			branch := folding
			if numFound > 1 {
				branch = &Folding{}
				folding.Branches = append(folding.Branches, branch)
			}
			branch.JoinPrefix = &Folding{}
			recursiveBacktrack(seq, v, w, i, k, branch.JoinPrefix)

			branch.JoinSuffix = &Folding{}
			branch.JoinSuffix.Pairs = append(branch.JoinSuffix.Pairs, Pair{I: k + 1, J: j})
			recursiveBacktrack(seq, v, w, k+2, j-1, branch.JoinSuffix)
		}
	}
	if numFound == 0 {
		for k := i; k < j; k++ {
			if v[i][j] == v[i][k]+v[k+1][j] {
				folding.JoinPrefix = &Folding{}
				recursiveBacktrack(seq, v, w, i, k, folding.JoinPrefix)
				folding.JoinSuffix = &Folding{}
				recursiveBacktrack(seq, v, w, k+1, j, folding.JoinSuffix)
				break
			}
		}
	}
}

func branchless(f *Folding) bool {
	return f.JoinPrefix == nil && len(f.Branches) == 0
}

func collapseTree(f *Folding) {
	for _, b := range f.Branches {
		collapseTree(b)
	}
	if f.JoinPrefix != nil {
		collapseTree(f.JoinPrefix)
		collapseTree(f.JoinSuffix)

		if branchless(f.JoinPrefix) && branchless(f.JoinSuffix) {
			f.Pairs = append(f.Pairs, f.JoinPrefix.Pairs...)
			f.Free = append(f.Free, f.JoinPrefix.Free...)
			f.JoinPrefix = nil
			f.Pairs = append(f.Pairs, f.JoinSuffix.Pairs...)
			f.Free = append(f.Free, f.JoinSuffix.Free...)
			f.JoinSuffix = nil
		}
	}
}
