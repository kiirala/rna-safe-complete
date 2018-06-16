package format

import "fmt"
import "strings"

import "keltainen.duckdns.org/rnafolding/base"

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

func Folding(seq *base.Sequence, pairs []int) string {
	c := &canvas{c: map[int]map[int]rune{}}
	recursiveFolding(seq, pairs, 0, 0, 0, len(seq.Bases)-1, c)
	return c.Draw()
}

func recursiveFolding(seq *base.Sequence, pairs []int, x, y, i, j int, c *canvas) {
	if i > j {
		c.Set(x, y+1, "\\")
		c.Set(x, y+2, "/")
		return
	}
	if pairs[i] == j {
		c.Set(x, y+1, seq.Bases[i].ToCode())
		c.Set(x, y+2, seq.Bases[j].ToCode())
		recursiveFolding(seq, pairs, x+1, y, i+1, j-1, c)
	} else if pairs[i] < 0 && pairs[j] < 0 && i != j {
		c.Set(x, y, seq.Bases[i].ToCode())
		c.Set(x, y+3, seq.Bases[j].ToCode())
		recursiveFolding(seq, pairs, x+1, y, i+1, j-1, c)
	} else if pairs[i] < 0 {
		c.Set(x, y, seq.Bases[i].ToCode())
		c.Set(x, y+3, "-")
		recursiveFolding(seq, pairs, x+1, y, i+1, j, c)
	} else {
		c.Set(x, y, "-")
		c.Set(x, y+3, seq.Bases[j].ToCode())
		recursiveFolding(seq, pairs, x+1, y, i, j-1, c)
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
