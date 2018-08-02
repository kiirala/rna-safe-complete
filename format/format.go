package format

import "fmt"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"
import "keltainen.duckdns.org/rnafolding/folding"

func Matrix(m [][]int) string {
	strs := make([]string, len(m))
	for i, arr := range m {
		strs[i] = fmt.Sprint(arr)
	}
	return "  " + strings.Join(strs, "\n  ")
}

type canvas struct {
	c map[int]map[int]rune
}

func DotBracket(pairs folding.FoldingPairs) string {
	out := make([]rune, len(pairs))
	for i, p := range pairs {
		if p < 0 {
			out[i] = '.'
		} else {
			if p > i {
				out[i] = '('
			} else {
				out[i] = ')'
			}
		}
	}
	return string(out)
}

func Folding(seq *base.Sequence, pairs folding.FoldingPairs) string {
	safety := make([]bool, len(pairs))
	for i := 0; i < len(safety); i++ {
		safety[i] = true
	}
	return FoldingWithSafety(seq, pairs, safety)
}

func FoldingWithSafety(seq *base.Sequence, pairs folding.FoldingPairs, safety []bool) string {
	c := &canvas{c: map[int]map[int]rune{}}
	c.Set(0, -1, ">")
	c.Set(0, 1, "<")
	recursiveFolding(seq, pairs, safety, 1, 0, 0, len(seq.Bases)-1, 0, c)
	return c.Draw()
}

func directionToVectors(dir int) (int, int, int, int) {
	var x, y int
	if dir&1 == 0 {
		x = 1
		y = 0
	} else {
		x = 0
		y = 1
	}
	if dir&2 == 0 {
		return x, y, -y, x
	}
	return -x, -y, y, -x
}

func rotateCCW(dir int) int {
	return (dir - 1) % 4
}

func rotateCW(dir int) int {
	return (dir + 1) % 4
}

func recursiveFolding(seq *base.Sequence, pairs []int, safety []bool, x, y, i, j, dir int, c *canvas) {
	majorx, majory, minorx, minory := directionToVectors(dir)
	if i > j {
		//c.Set(x-minorx, y-minory, "\\")
		//c.Set(x, y, "|")
		//c.Set(x+minorx, y+minory, "/")
		return
	}
	if pairs[i] == j {
		c.Set(x-minorx, y-minory, seq.Bases[i].CodeAndSafety(safety[i]))
		c.Set(x, y, "#")
		c.Set(x+minorx, y+minory, seq.Bases[j].CodeAndSafety(safety[j]))
		recursiveFolding(seq, pairs, safety, x+majorx, y+majory, i+1, j-1, dir, c)
	} else if pairs[i] >= 0 && pairs[j] >= 0 {
		if pairs[i]+1 > pairs[j]-1 {
			if pairs[i]-i > j-pairs[j] {
				recursiveFolding(seq, pairs, safety, x+3*majorx, y+3*majory, i, pairs[i], dir, c)
				recursiveFolding(seq, pairs, safety, x+majorx+2*minorx, y+majory+2*minory, pairs[j], j, rotateCW(dir), c)
			} else {
				recursiveFolding(seq, pairs, safety, x+majorx-2*minorx, y+majory-2*minory, i, pairs[i], rotateCCW(dir), c)
				recursiveFolding(seq, pairs, safety, x+3*majorx, y+3*majory, pairs[j], j, dir, c)
			}
		} else {
			recursiveFolding(seq, pairs, safety, x+majorx-2*minorx, y+majory-2*minory, i, pairs[i], rotateCCW(dir), c)
			recursiveFolding(seq, pairs, safety, x+majorx+2*minorx, y+majory+2*minory, pairs[j], j, rotateCW(dir), c)
			recursiveFolding(seq, pairs, safety, x+3*majorx, y+3*majory, pairs[i]+1, pairs[j]-1, dir, c)
		}
	} else if pairs[i] < 0 && pairs[j] < 0 && i != j {
		c.Set(x-minorx, y-minory, seq.Bases[i].CodeAndSafety(safety[i]))
		c.Set(x+minorx, y+minory, seq.Bases[j].CodeAndSafety(safety[j]))
		recursiveFolding(seq, pairs, safety, x+majorx, y+majory, i+1, j-1, dir, c)
	} else if pairs[i] < 0 {
		c.Set(x-minorx, y-minory, seq.Bases[i].CodeAndSafety(safety[i]))
		c.Set(x+minorx, y+minory, "-")
		recursiveFolding(seq, pairs, safety, x+majorx, y+majory, i+1, j, dir, c)
	} else {
		c.Set(x-minorx, y-minory, "-")
		c.Set(x+minorx, y+minory, seq.Bases[j].CodeAndSafety(safety[j]))
		recursiveFolding(seq, pairs, safety, x+majorx, y+majory, i, j-1, dir, c)
	}
}

func (c *canvas) Set(x, y int, str string) {
	line, ok := c.c[y]
	if !ok {
		line = make(map[int]rune)
		c.c[y] = line
	}
	for i, ch := range str {
		line[x+i] = ch
	}
}

func keyFromMap(m map[int]interface{}) int {
	for k := range m {
		return k
	}
	return -1
}

func (c *canvas) Draw() string {
	var miny int
	for k := range c.c {
		miny = k
		break
	}
	maxy := miny
	var minx int
	for k := range c.c[miny] {
		miny = k
		break
	}
	maxx := minx
	for y := range c.c {
		if y < miny {
			miny = y
		}
		if y > maxy {
			maxy = y
		}
		for x := range c.c[y] {
			if x < minx {
				minx = x
			}
			if x > maxx {
				maxx = x
			}
		}
	}

	lines := make([][]rune, maxy-miny+1)
	for y := 0; y < len(lines); y++ {
		lines[y] = make([]rune, maxx-minx+1)
		for x := 0; x < len(lines[y]); x++ {
			lines[y][x] = ' '
		}
	}

	for y := range c.c {
		for x := range c.c[y] {
			lines[y-miny][x-minx] = c.c[y][x]
		}
	}

	var out string
	for y := 0; y < len(lines); y++ {
		out += string(lines[y]) + "\n"
	}
	return out
}
