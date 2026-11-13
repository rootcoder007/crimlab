// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

// Made by tools/vendor.sh from the table of cli/jrun.cpp; do not edit.

package crimlab

// Vodtst runs vodtst() of the R package (arguments: day night ref? shift?; ? marks an optional one).
func Vodtst(args string) (string, error) { return Jrun("vodtst", args) }

// Orrisk runs orrisk() of the R package (arguments: or rate share; ? marks an optional one).
func Orrisk(args string) (string, error) { return Jrun("orrisk", args) }

// Outtst runs outtst() of the R package (arguments: q t w?; ? marks an optional one).
func Outtst(args string) (string, error) { return Jrun("outtst", args) }

// Dspbm runs dspbm() of the R package (arguments: pop contact force ref?; ? marks an optional one).
func Dspbm(args string) (string, error) { return Jrun("dspbm", args) }

// Dyadrt runs dyadrt() of the R package (arguments: count pop; ? marks an optional one).
func Dyadrt(args string) (string, error) { return Jrun("dyadrt", args) }

// Gtmp runs gtmp() of the R package (arguments: p q? hom? theta?; ? marks an optional one).
func Gtmp(args string) (string, error) { return Jrun("gtmp", args) }

// Cllr runs cllr() of the R package (arguments: same diff; ? marks an optional one).
func Cllr(args string) (string, error) { return Jrun("cllr", args) }

// Lrpost runs lrpost() of the R package (arguments: lr prior; ? marks an optional one).
func Lrpost(args string) (string, error) { return Jrun("lrpost", args) }

// Mixlr runs mixlr() of the R package (arguments: p np xp nd xd; ? marks an optional one).
func Mixlr(args string) (string, error) { return Jrun("mixlr", args) }

// Dbsrch runs dbsrch() of the R package (arguments: mp n pop?; ? marks an optional one).
func Dbsrch(args string) (string, error) { return Jrun("dbsrch", args) }

// Cllrmin runs cllrmin() of the R package (arguments: same diff; ? marks an optional one).
func Cllrmin(args string) (string, error) { return Jrun("cllrmin", args) }

// Lrcal runs lrcal() of the R package (arguments: same diff; ? marks an optional one).
func Lrcal(args string) (string, error) { return Jrun("lrcal", args) }

// Tipplt runs tipplt() of the R package (arguments: same diff; ? marks an optional one).
func Tipplt(args string) (string, error) { return Jrun("tipplt", args) }

// Eceplt runs eceplt() of the R package (arguments: same diff lo?; ? marks an optional one).
func Eceplt(args string) (string, error) { return Jrun("eceplt", args) }

// Kinidx runs kinidx() of the R package (arguments: g1 g2 freq rel?; ? marks an optional one).
func Kinidx(args string) (string, error) { return Jrun("kinidx", args) }

// Patidx runs patidx() of the R package (arguments: child father freq mother?; ? marks an optional one).
func Patidx(args string) (string, error) { return Jrun("patidx", args) }

// Clrsmp runs clrsmp() of the R package (arguments: clra reca clrb recb; ? marks an optional one).
func Clrsmp(args string) (string, error) { return Jrun("clrsmp", args) }

// Expbnd runs expbnd() of the R package (arguments: y m gamma ref?; ? marks an optional one).
func Expbnd(args string) (string, error) { return Jrun("expbnd", args) }

// Misrd runs misrd() of the R package (arguments: r1 r0 prev se sp; ? marks an optional one).
func Misrd(args string) (string, error) { return Jrun("misrd", args) }

// Rtmhs runs rtmhs() of the R package (arguments: x1 x2 cut w?; ? marks an optional one).
func Rtmhs(args string) (string, error) { return Jrun("rtmhs", args) }

// Incap runs incap() of the R package (arguments: lam q s share?; ? marks an optional one).
func Incap(args string) (string, error) { return Jrun("incap", args) }

// Incapr runs incapr() of the R package (arguments: lam t0 s r?; ? marks an optional one).
func Incapr(args string) (string, error) { return Jrun("incapr", args) }

// Sellab runs sellab() of the R package (arguments: y rel rule w?; ? marks an optional one).
func Sellab(args string) (string, error) { return Jrun("sellab", args) }

// Losbnd runs losbnd() of the R package (arguments: done open; ? marks an optional one).
func Losbnd(args string) (string, error) { return Jrun("losbnd", args) }

// Sntbnd runs sntbnd() of the R package (arguments: y z treat? w?; ? marks an optional one).
func Sntbnd(args string) (string, error) { return Jrun("sntbnd", args) }

// Hzdsel runs hzdsel() of the R package (arguments: s h l surh surl; ? marks an optional one).
func Hzdsel(args string) (string, error) { return Jrun("hzdsel", args) }

// Selftest runs selftest() of the R package (no arguments; ? marks an optional one).
func Selftest(args string) (string, error) { return Jrun("selftest", args) }

// Proof runs proof() of the R package (arguments: fn?; ? marks an optional one).
func Proof(args string) (string, error) { return Jrun("proof", args) }
