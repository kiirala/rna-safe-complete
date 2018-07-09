package safecomplete

func TrivialSafety(foldings [][]int) []bool {
	out := make([]bool, len(foldings[0]))
	for i := 0; i < len(out); i++ {
		out[i] = true
	}
	for _, f := range foldings {
		for i := 0; i < len(f); i++ {
			if f[i] != foldings[0][i] && (f[i] >= 0 || foldings[0][i] >= 0) {
				out[i] = false
			}
		}
	}
	return out
}
