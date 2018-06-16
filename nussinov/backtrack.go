package nussinov

import "keltainen.duckdns.org/rnafolding/base"

func Backtrack(seq *base.Sequence, v [][]int) (int, []int) {
	paired := make([]int, len(v))
	for i := 0; i < len(v); i++ {
		paired[i] = -1
	}
	recursiveBacktrack(seq, v, 0, len(v)-1, paired)
	return v[0][len(v)-1], paired
}

func JoinedBacktrack(seq *base.Sequence, v, w [][]int, i, j int) (int, []int) {
	paired := make([]int, len(v))
	pairs := 0
	for i := 0; i < len(v); i++ {
		paired[i] = -1
	}
	if seq.CanPair(i, j) {
		paired[i] = j
		paired[j] = i
		pairs++
	}
	recursiveBacktrack(seq, v, i+1, j-1, paired)
	pairs += v[i+1][j-1]
	complementaryBacktrack(seq, v, w, i, j, paired)
	pairs += w[i][j]
	return pairs, paired
}

func recursiveBacktrack(seq *base.Sequence, v [][]int, i, j int, paired []int) {
	if i >= j {
		return
	} else if seq.CanPair(i, j) && v[i][j] == v[i+1][j-1]+1 {
		paired[i] = j
		paired[j] = i
		recursiveBacktrack(seq, v, i+1, j-1, paired)
	} else {
		for k := i; k < j; k++ {
			if v[i][j] == v[i][k]+v[k+1][j] {
				recursiveBacktrack(seq, v, i, k, paired)
				recursiveBacktrack(seq, v, k+1, j, paired)
				break
			}
		}
	}
}

func complementaryBacktrack(seq *base.Sequence, v, w [][]int, i, j int, paired []int) {
	if j-i >= len(v)-2 {
		return
	}
	if i == 0 {
		recursiveBacktrack(seq, v, j+1, len(v)-1, paired)
	} else if j == len(v)-1 {
		recursiveBacktrack(seq, v, 0, i-1, paired)
	} else if seq.CanPair(i-1, j+1) && w[i][j] == w[i-1][j+1]+1 {
		paired[i-1] = j + 1
		paired[j+1] = i - 1
		complementaryBacktrack(seq, v, w, i-1, j+1, paired)
	} else {
		for k := j + 1; k < i+len(v)-1; k++ {
			if k < len(v) && w[i][j] == v[j+1][k]+w[i][k] {
				recursiveBacktrack(seq, v, j+1, k, paired)
				complementaryBacktrack(seq, v, w, i, k, paired)
				break
			} else if k >= len(v) && w[i][j] == w[k+1-len(v)][j]+v[k+1-len(v)][i-1] {
				complementaryBacktrack(seq, v, w, k+1-len(v), j, paired)
				recursiveBacktrack(seq, v, k+1-len(v), i-1, paired)
				break
			}
		}
	}
}
