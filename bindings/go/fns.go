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

// Impfpr runs impfpr() of the R package (arguments: p ppv fnr; ? marks an optional one).
func Impfpr(args string) (string, error) { return Jrun("impfpr", args) }

// Sqrtlw runs sqrtlw() of the R package (arguments: area units busy? metric? speed?; ? marks an optional one).
func Sqrtlw(args string) (string, error) { return Jrun("sqrtlw", args) }

// Fair runs fair() of the R package (arguments: score y grp cut? base? conf?; ? marks an optional one).
func Fair(args string) (string, error) { return Jrun("fair", args) }

// Fairdf runs fairdf() of the R package (arguments: pred grp? y? base? conf? nboot? seed?; ? marks an optional one).
func Fairdf(args string) (string, error) { return Jrun("fairdf", args) }

// Tmpaud runs tmpaud() of the R package (arguments: pred grp? period? city? base?; ? marks an optional one).
func Tmpaud(args string) (string, error) { return Jrun("tmpaud", args) }

// Lipaud runs lipaud() of the R package (arguments: h d? pairs? lip?; ? marks an optional one).
func Lipaud(args string) (string, error) { return Jrun("lipaud", args) }

// Nsbnd runs nsbnd() of the R package (arguments: pobs a b; ? marks an optional one).
func Nsbnd(args string) (string, error) { return Jrun("nsbnd", args) }

// Rnkres runs rnkres() of the R package (arguments: est hw; ? marks an optional one).
func Rnkres(args string) (string, error) { return Jrun("rnkres", args) }

// Nordet runs nordet() of the R package (arguments: x off radius p? draw? seed?; ? marks an optional one).
func Nordet(args string) (string, error) { return Jrun("nordet", args) }

// Fairsim runs fairsim() of the R package (arguments: n? groups? props? nareas? base? bias? seed?; ? marks an optional one).
func Fairsim(args string) (string, error) { return Jrun("fairsim", args) }

// Oaxbd runs oaxbd() of the R package (arguments: y x grp ref? conf?; ? marks an optional one).
func Oaxbd(args string) (string, error) { return Jrun("oaxbd", args) }

// Tlvc runs tlvc() of the R package (arguments: x item; ? marks an optional one).
func Tlvc(args string) (string, error) { return Jrun("tlvc", args) }

// Tlvlr runs tlvlr() of the R package (arguments: xa xb na nb bg; ? marks an optional one).
func Tlvlr(args string) (string, error) { return Jrun("tlvlr", args) }

// Aorwt runs aorwt() of the R package (arguments: from to bw? nb? origin? off? tt? toff? lo? hi? each?; ? marks an optional one).
func Aorwt(args string) (string, error) { return Jrun("aorwt", args) }

// Bdig runs bdig() of the R package (arguments: x nd? dp? sgn?; ? marks an optional one).
func Bdig(args string) (string, error) { return Jrun("bdig", args) }

// Bnft runs bnft() of the R package (arguments: dig nd? nsim? seed?; ? marks an optional one).
func Bnft(args string) (string, error) { return Jrun("bnft", args) }

// Skp runs skp() of the R package (arguments: x s2 alpha? nsim? seed?; ? marks an optional one).
func Skp(args string) (string, error) { return Jrun("skp", args) }

// Rptc runs rptc() of the R package (arguments: unit year rep y0 y1 minrep?; ? marks an optional one).
func Rptc(args string) (string, error) { return Jrun("rptc", args) }

// Rswp runs rswp() of the R package (arguments: d gam? alpha?; ? marks an optional one).
func Rswp(args string) (string, error) { return Jrun("rswp", args) }

// Rshl runs rshl() of the R package (arguments: d gam?; ? marks an optional one).
func Rshl(args string) (string, error) { return Jrun("rshl", args) }

// Lwss runs lwss() of the R package (arguments: x y f? iter? delta?; ? marks an optional one).
func Lwss(args string) (string, error) { return Jrun("lwss", args) }

// Chmb runs chmb() of the R package (arguments: a b wopt wval coarse?; ? marks an optional one).
func Chmb(args string) (string, error) { return Jrun("chmb", args) }

// Gstr runs gstr() of the R package (arguments: cnt ex? rad? nsim? seed?; ? marks an optional one).
func Gstr(args string) (string, error) { return Jrun("gstr", args) }

// Sptst runs sptst() of the R package (arguments: a b nboot? seed? alpha? scale?; ? marks an optional one).
func Sptst(args string) (string, error) { return Jrun("sptst", args) }

// Tfrp runs tfrp() of the R package (arguments: par deffect? tdist? ymax? nsim? seed? hl? nq?; ? marks an optional one).
func Tfrp(args string) (string, error) { return Jrun("tfrp", args) }

// Zgll runs zgll() of the R package (arguments: y x z? par prob off?; ? marks an optional one).
func Zgll(args string) (string, error) { return Jrun("zgll", args) }

// Zgem runs zgem() of the R package (arguments: y x z? par prob off? maxit? tol?; ? marks an optional one).
func Zgem(args string) (string, error) { return Jrun("zgem", args) }

// Mvtll runs mvtll() of the R package (arguments: y x beta sig prob; ? marks an optional one).
func Mvtll(args string) (string, error) { return Jrun("mvtll", args) }

// Mvtem runs mvtem() of the R package (arguments: y x beta sig prob maxit? tol?; ? marks an optional one).
func Mvtem(args string) (string, error) { return Jrun("mvtem", args) }

// Ntdeg runs ntdeg() of the R package (arguments: from to w? nodes? dir?; ? marks an optional one).
func Ntdeg(args string) (string, error) { return Jrun("ntdeg", args) }

// Ntbtw runs ntbtw() of the R package (arguments: from to nodes? dir? norm?; ? marks an optional one).
func Ntbtw(args string) (string, error) { return Jrun("ntbtw", args) }

// Ntcls runs ntcls() of the R package (arguments: from to nodes? dir?; ? marks an optional one).
func Ntcls(args string) (string, error) { return Jrun("ntcls", args) }

// Ntevc runs ntevc() of the R package (arguments: from to w? nodes? tol? maxit?; ? marks an optional one).
func Ntevc(args string) (string, error) { return Jrun("ntevc", args) }

// Ntcore runs ntcore() of the R package (arguments: from to nodes?; ? marks an optional one).
func Ntcore(args string) (string, error) { return Jrun("ntcore", args) }

// Ntcmp runs ntcmp() of the R package (arguments: from to nodes? dir?; ? marks an optional one).
func Ntcmp(args string) (string, error) { return Jrun("ntcmp", args) }

// Ntlvn runs ntlvn() of the R package (arguments: from to w? nodes? res? seed?; ? marks an optional one).
func Ntlvn(args string) (string, error) { return Jrun("ntlvn", args) }

// Ntmod runs ntmod() of the R package (arguments: from to memb w? nodes? res?; ? marks an optional one).
func Ntmod(args string) (string, error) { return Jrun("ntmod", args) }

// Ntbrk runs ntbrk() of the R package (arguments: from to grp nodes? dir?; ? marks an optional one).
func Ntbrk(args string) (string, error) { return Jrun("ntbrk", args) }

// Ntchg runs ntchg() of the R package (arguments: from to w? nodes? tol? maxit?; ? marks an optional one).
func Ntchg(args string) (string, error) { return Jrun("ntchg", args) }

// Coofg runs coofg() of the R package (arguments: inc off; ? marks an optional one).
func Coofg(args string) (string, error) { return Jrun("coofg", args) }

// Relchk runs relchk() of the R package (arguments: x value type? group? ref? event? digits? tol?; ? marks an optional one).
func Relchk(args string) (string, error) { return Jrun("relchk", args) }

// Catchk runs catchk() of the R package (arguments: x ref? nmin?; ? marks an optional one).
func Catchk(args string) (string, error) { return Jrun("catchk", args) }

// Hvdist runs hvdist() of the R package (arguments: from to r?; ? marks an optional one).
func Hvdist(args string) (string, error) { return Jrun("hvdist", args) }

// Eqproj runs eqproj() of the R package (arguments: x y? lon0? lat0? r?; ? marks an optional one).
func Eqproj(args string) (string, error) { return Jrun("eqproj", args) }

// Eqinv runs eqinv() of the R package (arguments: x y? lon0 lat0 r?; ? marks an optional one).
func Eqinv(args string) (string, error) { return Jrun("eqinv", args) }

// Cgstat runs cgstat() of the R package (arguments: x y? w? k? lonlat?; ? marks an optional one).
func Cgstat(args string) (string, error) { return Jrun("cgstat", args) }

// Cgmedn runs cgmedn() of the R package (arguments: x y? lonlat?; ? marks an optional one).
func Cgmedn(args string) (string, error) { return Jrun("cgmedn", args) }

// Webpt runs webpt() of the R package (arguments: x y? w? start? tol? maxit? lonlat?; ? marks an optional one).
func Webpt(args string) (string, error) { return Jrun("webpt", args) }

// Sdell runs sdell() of the R package (arguments: x y? w? k? c? lonlat?; ? marks an optional one).
func Sdell(args string) (string, error) { return Jrun("sdell", args) }

// Cvhull runs cvhull() of the R package (arguments: x y? lonlat?; ? marks an optional one).
func Cvhull(args string) (string, error) { return Jrun("cvhull", args) }

// Nnidx runs nnidx() of the R package (arguments: x y? win kmax? nsim? seed? dist? metric?; ? marks an optional one).
func Nnidx(args string) (string, error) { return Jrun("nnidx", args) }

// Ripk runs ripk() of the R package (arguments: x y? win r corr? metric?; ? marks an optional one).
func Ripk(args string) (string, error) { return Jrun("ripk", args) }

// Ripenv runs ripenv() of the R package (arguments: x y? win r corr? nsim? seed?; ? marks an optional one).
func Ripenv(args string) (string, error) { return Jrun("ripenv", args) }

// Fzmode runs fzmode() of the R package (arguments: x y? r count? metric? lonlat?; ? marks an optional one).
func Fzmode(args string) (string, error) { return Jrun("fzmode", args) }

// Kmpart runs kmpart() of the R package (arguments: x y? k init? cell? upd? maxit? seed? count? lonlat?; ? marks an optional one).
func Kmpart(args string) (string, error) { return Jrun("kmpart", args) }

// Nnhclu runs nnhclu() of the R package (arguments: x y? area? fixd? p? minpts? minhi? maxord? k? c? count? lonlat?; ? marks an optional one).
func Nnhclu(args string) (string, error) { return Jrun("nnhclu", args) }

// Stacel runs stacel() of the R package (arguments: x y? win nx? ny? r top? minpts? k? c? count?; ? marks an optional one).
func Stacel(args string) (string, error) { return Jrun("stacel", args) }

// Mantst runs mantst() of the R package (arguments: x y? t cs? ct? nsim? seed? metric? lonlat?; ? marks an optional one).
func Mantst(args string) (string, error) { return Jrun("mantst", args) }

// Stmavg runs stmavg() of the R package (arguments: x y? t span? lonlat?; ? marks an optional one).
func Stmavg(args string) (string, error) { return Jrun("stmavg", args) }

// Jtdecay runs jtdecay() of the R package (arguments: d type par tab?; ? marks an optional one).
func Jtdecay(args string) (string, error) { return Jrun("jtdecay", args) }

// Jtsurf runs jtsurf() of the R package (arguments: x y? win nx ny type par tab? metric? mode?; ? marks an optional one).
func Jtsurf(args string) (string, error) { return Jrun("jtsurf", args) }

// Jtcal runs jtcal() of the R package (arguments: d t h w? kern? refl?; ? marks an optional one).
func Jtcal(args string) (string, error) { return Jrun("jtcal", args) }

// Jtacc runs jtacc() of the R package (arguments: prob win base; ? marks an optional one).
func Jtacc(args string) (string, error) { return Jrun("jtacc", args) }

// Bjtc runs bjtc() of the R package (arguments: pj od dest; ? marks an optional one).
func Bjtc(args string) (string, error) { return Jrun("bjtc", args) }

// Spwt runs spwt() of the R package (arguments: x? y? type? k? d? nb? idw? style? metric? lonlat?; ? marks an optional one).
func Spwt(args string) (string, error) { return Jrun("spwt", args) }

// Gearc runs gearc() of the R package (arguments: x w nsim? seed?; ? marks an optional one).
func Gearc(args string) (string, error) { return Jrun("gearc", args) }

// Gorg runs gorg() of the R package (arguments: x w nsim? seed?; ? marks an optional one).
func Gorg(args string) (string, error) { return Jrun("gorg", args) }

// Lmorn runs lmorn() of the R package (arguments: x w nsim? seed?; ? marks an optional one).
func Lmorn(args string) (string, error) { return Jrun("lmorn", args) }

// Locgi runs locgi() of the R package (arguments: x w star?; ? marks an optional one).
func Locgi(args string) (string, error) { return Jrun("locgi", args) }

// Dirmn runs dirmn() of the R package (arguments: x? y? origin? w? angle? lonlat?; ? marks an optional one).
func Dirmn(args string) (string, error) { return Jrun("dirmn", args) }

// Hdbang runs hdbang() of the R package (arguments: x y? v w? nn? ang? ntrip? iter? lonlat?; ? marks an optional one).
func Hdbang(args string) (string, error) { return Jrun("hdbang", args) }

// Crwalk runs crwalk() of the R package (arguments: x y? h? nsim? seed? lonlat?; ? marks an optional one).
func Crwalk(args string) (string, error) { return Jrun("crwalk", args) }

// Rnnhcl runs rnnhcl() of the R package (arguments: x y? base bw? h p? minpts? k? c? lonlat?; ? marks an optional one).
func Rnnhcl(args string) (string, error) { return Jrun("rnnhcl", args) }

// Txcal runs txcal() of the R package (arguments: conc resp w? unk? k?; ? marks an optional one).
func Txcal(args string) (string, error) { return Jrun("txcal", args) }

// Lodimp runs lodimp() of the R package (arguments: x lod meth?; ? marks an optional one).
func Lodimp(args string) (string, error) { return Jrun("lodimp", args) }

// Pmrrat runs pmrrat() of the R package (arguments: cen per cut?; ? marks an optional one).
func Pmrrat(args string) (string, error) { return Jrun("pmrrat", args) }

// Etohlr runs etohlr() of the R package (arguments: eth prop? but? etg? thr?; ? marks an optional one).
func Etohlr(args string) (string, error) { return Jrun("etohlr", args) }

// Mtxrel runs mtxrel() of the R package (arguments: base dil days? stage? dmax? dcap? smax? scap?; ? marks an optional one).
func Mtxrel(args string) (string, error) { return Jrun("mtxrel", args) }

// Detrsp runs detrsp() of the R package (arguments: x ben sanc p; ? marks an optional one).
func Detrsp(args string) (string, error) { return Jrun("detrsp", args) }

// Dtrid runs dtrid() of the R package (arguments: p s c y? log? rho?; ? marks an optional one).
func Dtrid(args string) (string, error) { return Jrun("dtrid", args) }

// Inspg runs inspg() of the R package (arguments: a x steps? dt?; ? marks an optional one).
func Inspg(args string) (string, error) { return Jrun("inspg", args) }

// Insphs runs insphs() of the R package (arguments: temp? cost? steps? rate? x0?; ? marks an optional one).
func Insphs(args string) (string, error) { return Jrun("insphs", args) }

// Colarr runs colarr() of the R package (arguments: a e bg; ? marks an optional one).
func Colarr(args string) (string, error) { return Jrun("colarr", args) }

// Crsev runs crsev() of the R package (arguments: x off w unit? time? pop? per? base? sent?; ? marks an optional one).
func Crsev(args string) (string, error) { return Jrun("crsev", args) }

// Zcomp runs zcomp() of the R package (arguments: x unit cat w?; ? marks an optional one).
func Zcomp(args string) (string, error) { return Jrun("zcomp", args) }

// Cconc runs cconc() of the R package (arguments: x w? by? fracs? pop? nsim? seed?; ? marks an optional one).
func Cconc(args string) (string, error) { return Jrun("cconc", args) }

// Hilltl runs hilltl() of the R package (arguments: x xmin? discrete? approx? mintl?; ? marks an optional one).
func Hilltl(args string) (string, error) { return Jrun("hilltl", args) }

// Crcap runs crcap() of the R package (arguments: x kappa? int?; ? marks an optional one).
func Crcap(args string) (string, error) { return Jrun("crcap", args) }

// Dkfig runs dkfig() of the R package (arguments: v? r? amax? bmax? thr? va? ra? rb? lam? k?; ? marks an optional one).
func Dkfig(args string) (string, error) { return Jrun("dkfig", args) }

// Bndxw runs bndxw() of the R package (arguments: from to x unit w?; ? marks an optional one).
func Bndxw(args string) (string, error) { return Jrun("bndxw", args) }

// Synexp runs synexp() of the R package (arguments: y x ax pop; ? marks an optional one).
func Synexp(args string) (string, error) { return Jrun("synexp", args) }

// Clrilr runs clrilr() of the R package (arguments: x type? ps?; ? marks an optional one).
func Clrilr(args string) (string, error) { return Jrun("clrilr", args) }

// Preslr runs preslr() of the R package (arguments: x m1 s1 m2 s2; ? marks an optional one).
func Preslr(args string) (string, error) { return Jrun("preslr", args) }

// Tpdmc runs tpdmc() of the R package (arguments: pres? decay? mum? states? start? steps? seed?; ? marks an optional one).
func Tpdmc(args string) (string, error) { return Jrun("tpdmc", args) }

// Tpbhm runs tpbhm() of the R package (arguments: y x grp? pm? psd? level?; ? marks an optional one).
func Tpbhm(args string) (string, error) { return Jrun("tpbhm", args) }

// Enpath runs enpath() of the R package (arguments: x y fam? alpha? lambda? nlam? lmr? fold? tol? maxit? scr?; ? marks an optional one).
func Enpath(args string) (string, error) { return Jrun("enpath", args) }

// Lgtreg runs lgtreg() of the R package (arguments: x y tol? maxit?; ? marks an optional one).
func Lgtreg(args string) (string, error) { return Jrun("lgtreg", args) }

// Rfor runs rfor() of the R package (arguments: x y newx? ntree? mtry? minnode? frac? repl? hon? imb? seed?; ? marks an optional one).
func Rfor(args string) (string, error) { return Jrun("rfor", args) }

// Crchyp runs crchyp() of the R package (arguments: x y home?; ? marks an optional one).
func Crchyp(args string) (string, error) { return Jrun("crchyp", args) }

// Hitscr runs hitscr() of the R package (arguments: s home; ? marks an optional one).
func Hitscr(args string) (string, error) { return Jrun("hitscr", args) }

// Rtmfit runs rtmfit() of the R package (arguments: y lay cell? rad? kind? level?; ? marks an optional one).
func Rtmfit(args string) (string, error) { return Jrun("rtmfit", args) }

// Agemix runs agemix() of the R package (arguments: curves? shares? ages? panel? nperm? seed?; ? marks an optional one).
func Agemix(args string) (string, error) { return Jrun("agemix", args) }

// Cumrsk runs cumrsk() of the R package (arguments: age ev pop atrisk? dth? level?; ? marks an optional one).
func Cumrsk(args string) (string, error) { return Jrun("cumrsk", args) }

// Jpipe runs jpipe() of the R package (arguments: num den win? lag? grp?; ? marks an optional one).
func Jpipe(args string) (string, error) { return Jrun("jpipe", args) }

// Smdbal runs smdbal() of the R package (arguments: x d w? ps? thr?; ? marks an optional one).
func Smdbal(args string) (string, error) { return Jrun("smdbal", args) }

// Evlx runs evlx() of the R package (arguments: est lo? hi? type? rare? sd?; ? marks an optional one).
func Evlx(args string) (string, error) { return Jrun("evlx", args) }

// Jdgiv runs jdgiv() of the R package (arguments: z? d y w? dsh? deff?; ? marks an optional one).
func Jdgiv(args string) (string, error) { return Jrun("jdgiv", args) }

// Mmshr runs mmshr() of the R package (arguments: sent min grp?; ? marks an optional one).
func Mmshr(args string) (string, error) { return Jrun("mmshr", args) }

// Hwkres runs hwkres() of the R package (arguments: times hor? par? kernel? base? per?; ? marks an optional one).
func Hwkres(args string) (string, error) { return Jrun("hwkres", args) }

// Hwksim runs hwksim() of the R package (arguments: par hor kernel? base? per? win? sigma? seed? nmax?; ? marks an optional one).
func Hwksim(args string) (string, error) { return Jrun("hwksim", args) }

// Stpcp runs stpcp() of the R package (arguments: y x? minseg? nperm? seed?; ? marks an optional one).
func Stpcp(args string) (string, error) { return Jrun("stpcp", args) }

// Selftest runs selftest() of the R package (no arguments; ? marks an optional one).
func Selftest(args string) (string, error) { return Jrun("selftest", args) }

// Proof runs proof() of the R package (arguments: fn?; ? marks an optional one).
func Proof(args string) (string, error) { return Jrun("proof", args) }
