// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

package crimlab

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
	"testing"
)

var cut = regexp.MustCompile(`[\[\]{},:]`)

// same compares two JSON texts as cli/test.sh does: token by token, a
// number to 1e-12 relative to the larger (absolute when the reference
// is 0), since other C libraries and fused multiply-adds change the last
// bits.
func same(got, want string) bool {
	a, b := cut.Split(got, -1), cut.Split(want, -1)
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] == b[i] {
			continue
		}
		x, e1 := strconv.ParseFloat(a[i], 64)
		y, e2 := strconv.ParseFloat(b[i], 64)
		if e1 != nil || e2 != nil {
			return false
		}
		m := math.Max(math.Abs(x), math.Abs(y))
		if y == 0 {
			m = 1
		}
		if math.Abs(x-y) > 1e-12*m {
			return false
		}
	}
	return true
}

// The table of the command-line program (cli/test.tsv: a function, its
// arguments, the result the Python package gives), when the package sits
// in the repository.
func TestTable(t *testing.T) {
	f, err := os.Open("../../cli/test.tsv")
	if err != nil {
		t.Skip("no ../../cli/test.tsv (the module outside the repository)")
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	n := 0
	for sc.Scan() {
		p := strings.Split(sc.Text(), "\t")
		if len(p) != 3 {
			continue
		}
		got, err := Jrun(p[0], p[1])
		if err != nil {
			t.Errorf("%s %s: %v", p[0], p[1], err)
		} else if !same(got, p[2]) {
			t.Errorf("%s %s\n got:  %s\n want: %s", p[0], p[1], got, p[2])
		}
		n++
	}
	if n < 30 {
		t.Fatalf("only %d lines in cli/test.tsv", n)
	}
}

func TestResults(t *testing.T) {
	r, err := Orrisk(`{"or": 2, "rate": 0.1, "share": 0.3}`)
	if err != nil || !same(r, `{"exposed":0.14752049189933314,"unexposed":0.079634074900285817,"rr":1.8524795081006669,"or":2}`) {
		t.Fatalf("orrisk: %s %v", r, err)
	}
	r, err = Selftest("")
	if err != nil || !strings.Contains(r, `"nfail":0,"report":""`) {
		t.Fatalf("selftest: %s %v", r, err)
	}
	for _, n := range []string{"vodtst", "orrisk", "sntbnd", "kinidx", "selftest", "proof"} {
		if !strings.Contains(Jfns(), `"`+n+`":`) {
			t.Errorf("Jfns has no %s", n)
		}
	}
}

func TestErrors(t *testing.T) {
	_, err := Orrisk(`{"or": 2, "rate": 1, "share": 0.3}`)
	if err == nil || err.Error() != "orrisk: `rate` must be one number strictly between 0 and 1." {
		t.Fatalf("orrisk: %v", err)
	}
	if _, err := Jrun("nope", "{}"); err == nil || !strings.Contains(err.Error(), "no function nope") {
		t.Fatalf("nope: %v", err)
	}
	if _, err := Losbnd("[1"); err == nil {
		t.Fatal("losbnd [1: no error")
	}
	if _, err := Jrun("cllr", "{\"same\":\x00}"); err == nil {
		t.Fatal("NUL: no error")
	}
}

func TestLarge(t *testing.T) {
	// longer than the first buffer: the call is made again with room
	v := make([]string, 19999)
	for i := range v {
		v[i] = fmt.Sprint(float64(i+1) / 7)
	}
	s := strings.Join(v, ",")
	r, err := Lrcal(`{"same": [` + s + `], "diff": [` + s + `]}`)
	if err != nil || len(r) <= 1<<16 {
		t.Fatalf("lrcal: %d bytes, %v", len(r), err)
	}
}

func ExampleLosbnd() {
	r, err := Losbnd(`{"done": [30, 45, 60, 90, 120], "open": [100, 150, 200]}`)
	if err != nil {
		panic(err)
	}
	fmt.Println(r)
	// Output: {"mean":69,"age":150,"lower":99.375,"bias":30.375}
}
