package safecomplete

import "fmt"

import "keltainen.duckdns.org/rnafolding/base"

func BacktrackAll(seq *base.Sequence, v, w [][]int) ([][]int, []string) {
	foldings := make([][]int, 1)
	foldings[0] = make([]int, len(v))
	for i := 0; i < len(v); i++ {
		foldings[0][i] = -1
	}
	makings := []string{""}
	foldings, makings = recursiveBacktrack(seq, v, w, 0, len(v)-1, foldings, makings, 0, false)
	return foldings, makings
}

func recursiveBacktrack(seq *base.Sequence, v, w [][]int, i, j int, foldings [][]int, makings []string, foldid int, onlyPairing bool) ([][]int, []string) {
	if i >= j {
		return foldings, makings
	}
	cachedFolding := make([]int, len(foldings[foldid]))
	copy(cachedFolding, foldings[foldid])
	cachedMaking := makings[foldid]
	numFound := 0
	if seq.CanPair(i, j) && v[i][j] == v[i+1][j-1]+1 {
		foldings[foldid][i] = j
		foldings[foldid][j] = i
		numFound++
		makings[foldid] += fmt.Sprintf(" p(%d,%d)", i, j)
		foldings, makings = recursiveBacktrack(seq, v, w, i+1, j-1, foldings, makings, foldid, false)
	}
	if onlyPairing {
		if numFound != 1 {
			makings[foldid] += fmt.Sprintf(" !!!(%d,%d)!!!", i, j)
		}
		return foldings, makings
	}
	for k := i; k < j; k++ {
		if v[i][j] == v[i][k]+w[k+1][j] && w[k+1][j] > 0 {
			numFound++
			newFoldID := foldid
			if numFound > 1 {
				newFolding := make([]int, len(foldings[0]))
				copy(newFolding, cachedFolding)
				foldings = append(foldings, newFolding)
				makings = append(makings, cachedMaking)
				newFoldID = len(foldings) - 1
			}
			makings[newFoldID] += fmt.Sprintf(" j(%d,%d)(%d,%d)", i, k, k+1, j)
			foldings, makings = recursiveBacktrack(seq, v, w, i, k, foldings, makings, newFoldID, false)
			foldings, makings = recursiveBacktrack(seq, v, w, k+1, j, foldings, makings, newFoldID, true)
		}
	}
	if numFound == 0 {
		for k := i; k < j; k++ {
			if v[i][j] == v[i][k]+v[k+1][j] {
				makings[foldid] += fmt.Sprintf(" x(%d,%d)(%d,%d)", i, k, k+1, j)
				foldings, makings = recursiveBacktrack(seq, v, w, i, k, foldings, makings, foldid, false)
				foldings, makings = recursiveBacktrack(seq, v, w, k+1, j, foldings, makings, foldid, false)
				break
			}
		}
	}
	return foldings, makings
}
