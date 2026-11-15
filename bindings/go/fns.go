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

// Cntreg runs cntreg() of the R package (arguments: y grp? adj? expo? clus? ref? gee? cadj? level?; ? marks an optional one).
func Cntreg(args string) (string, error) { return Jrun("cntreg", args) }

// Grprr runs grprr() of the R package (arguments: x n? grp? ref? type? per? level?; ? marks an optional one).
func Grprr(args string) (string, error) { return Jrun("grprr", args) }

// Dastd runs dastd() of the R package (arguments: x n strat grp? std ref? per? level?; ? marks an optional one).
func Dastd(args string) (string, error) { return Jrun("dastd", args) }

// Sirex runs sirex() of the R package (arguments: x n strat grp rates? level?; ? marks an optional one).
func Sirex(args string) (string, error) { return Jrun("sirex", args) }

// Ctabx runs ctabx() of the R package (arguments: row col? w? agg? correct? exact? level?; ? marks an optional one).
func Ctabx(args string) (string, error) { return Jrun("ctabx", args) }

// Calaud runs calaud() of the R package (arguments: area risk? outcome? group? population? nperm? seed?; ? marks an optional one).
func Calaud(args string) (string, error) { return Jrun("calaud", args) }

// Scrdsp runs scrdsp() of the R package (arguments: score group? ref? conf?; ? marks an optional one).
func Scrdsp(args string) (string, error) { return Jrun("scrdsp", args) }

// Jmono runs jmono() of the R package (arguments: judge d y? grp? w? lo? hi? minn? adj? alpha?; ? marks an optional one).
func Jmono(args string) (string, error) { return Jrun("jmono", args) }

// Mndla runs mndla() of the R package (arguments: days id? by? w? thr? hours? miss? alert? amin? elig? excl? level?; ? marks an optional one).
func Mndla(args string) (string, error) { return Jrun("mndla", args) }

// Yoypx runs yoypx() of the R package (arguments: x prd by? lag? units? minbase? level? rnd? kind? adj? alpha?; ? marks an optional one).
func Yoypx(args string) (string, error) { return Jrun("yoypx", args) }

// Lostay runs lostay() of the R package (arguments: days w? breaks? thr? flag? level?; ? marks an optional one).
func Lostay(args string) (string, error) { return Jrun("lostay", args) }

// Dktflw runs dktflw() of the R package (arguments: recv? disp? pend? prd? win? level? span? mtime? shares? floors? mix? scores?; ? marks an optional one).
func Dktflw(args string) (string, error) { return Jrun("dktflw", args) }

// Stkflw runs stkflw() of the R package (arguments: people? days? adp? alos? t? prd? from? to? expo? per? base? rel? early?; ? marks an optional one).
func Stkflw(args string) (string, error) { return Jrun("stkflw", args) }

// Bklog runs bklog() of the R package (arguments: arrive dispose horizon? target?; ? marks an optional one).
func Bklog(args string) (string, error) { return Jrun("bklog", args) }

// Trnmat runs trnmat() of the R package (arguments: state id time?; ? marks an optional one).
func Trnmat(args string) (string, error) { return Jrun("trnmat", args) }

// Lftab runs lftab() of the R package (arguments: breaks events cens? level? ci?; ? marks an optional one).
func Lftab(args string) (string, error) { return Jrun("lftab", args) }

// Bandx runs bandx() of the R package (arguments: x cnts? rule? cap? fac? floor? clsup? intsca? drop?; ? marks an optional one).
func Bandx(args string) (string, error) { return Jrun("bandx", args) }

// Pbbnd runs pbbnd() of the R package (arguments: x rnd? kind? suplim? supmar? prev? pop? prepop? per?; ? marks an optional one).
func Pbbnd(args string) (string, error) { return Jrun("pbbnd", args) }

// Hwkfit runs hwkfit() of the R package (arguments: times hor kernel? base? per? start? res? seed? maxit?; ? marks an optional one).
func Hwkfit(args string) (string, error) { return Jrun("hwkfit", args) }

// Trndrr runs trndrr() of the R package (arguments: y x? expo? seas? alt? level?; ? marks an optional one).
func Trndrr(args string) (string, error) { return Jrun("trndrr", args) }

// Cusrv runs cusrv() of the R package (arguments: y expd? start? k? h? trans? reset?; ? marks an optional one).
func Cusrv(args string) (string, error) { return Jrun("cusrv", args) }

// Xcser runs xcser() of the R package (arguments: x y lagmax? dif? level?; ? marks an optional one).
func Xcser(args string) (string, error) { return Jrun("xcser", args) }

// Kmrec runs kmrec() of the R package (arguments: time event? grp? unit? end? level?; ? marks an optional one).
func Kmrec(args string) (string, error) { return Jrun("kmrec", args) }

// Hypq runs hypq() of the R package (arguments: lam pref mu? queue? tt?; ? marks an optional one).
func Hypq(args string) (string, error) { return Jrun("hypq", args) }

// Sdbsim runs sdbsim() of the R package (arguments: n steps mode? a0? omega? eta? gam? theta? dt? b? rho? dx? dif? eps? z? sup? seed?; ? marks an optional one).
func Sdbsim(args string) (string, error) { return Jrun("sdbsim", args) }

// Sdbthr runs sdbthr() of the R package (arguments: a0 omega eps dif gam eta z? k?; ? marks an optional one).
func Sdbthr(args string) (string, error) { return Jrun("sdbthr", args) }

// Etasf runs etasf() of the R package (arguments: times hor mag? m0? start? maxit?; ? marks an optional one).
func Etasf(args string) (string, error) { return Jrun("etasf", args) }

// Mvhwk runs mvhwk() of the R package (arguments: times type hor beta? start? maxit?; ? marks an optional one).
func Mvhwk(args string) (string, error) { return Jrun("mvhwk", args) }

// Sthwk runs sthwk() of the R package (arguments: times x y hor win? par? at? start? maxit?; ? marks an optional one).
func Sthwk(args string) (string, error) { return Jrun("sthwk", args) }

// Didfe runs didfe() of the R package (arguments: y unit time trt adj? clus? w? fe? fam? expo? loo? bacon? mtch? exact? calip? level?; ? marks an optional one).
func Didfe(args string) (string, error) { return Jrun("didfe", args) }

// Evtst runs evtst() of the R package (arguments: y unit time first leads? lags? ref? bin? adj? clus? w? ma? norm? mbar? target? draws? seed? level?; ? marks an optional one).
func Evtst(args string) (string, error) { return Jrun("evtst", args) }

// Stgdid runs stgdid() of the R package (arguments: y unit time first? trt? meth? ctrl? adj? est? clus? reps? seed? level?; ? marks an optional one).
func Stgdid(args string) (string, error) { return Jrun("stgdid", args) }

// Synthc runs synthc() of the R package (arguments: y unit time first plac?; ? marks an optional one).
func Synthc(args string) (string, error) { return Jrun("synthc", args) }

// Sdidx runs sdidx() of the R package (arguments: y unit time first reps? seed? level?; ? marks an optional one).
func Sdidx(args string) (string, error) { return Jrun("sdidx", args) }

// Cbatt runs cbatt() of the R package (arguments: y d x est? tgt? ps? k? calip? nsub? reps? seed? level?; ? marks an optional one).
func Cbatt(args string) (string, error) { return Jrun("cbatt", args) }

// Svymx runs svymx() of the R package (arguments: y w? strata? psu? fpc? stat? den? by? ref? level?; ? marks an optional one).
func Svymx(args string) (string, error) { return Jrun("svymx", args) }

// Wcbx runs wcbx() of the R package (arguments: y unit time trt? adj? clus? meth? reps? wtype? ties? stat? ntrt? adopt? eff? alpha? power? seed? level?; ? marks an optional one).
func Wcbx(args string) (string, error) { return Jrun("wcbx", args) }

// Wystep runs wystep() of the R package (arguments: obs null; ? marks an optional one).
func Wystep(args string) (string, error) { return Jrun("wystep", args) }

// Dmlx runs dmlx() of the R package (arguments: y d x model? fold z? s? t? score? lrn? lrnd? alpha? nlam? nfin? tol? scr? hon? ntree? mtry? trim? clu? pred? stk? nstk? disc? seed?; ? marks an optional one).
func Dmlx(args string) (string, error) { return Jrun("dmlx", args) }

// Suplrn runs suplrn() of the R package (arguments: x y fam? lrns? nfold? fold? disc? newx? alpha? nlam? nfin? ntree? mtry? seed?; ? marks an optional one).
func Suplrn(args string) (string, error) { return Jrun("suplrn", args) }

// Knxnr runs knxnr() of the R package (arguments: x y? t ds dt nsim? seed? lonlat?; ? marks an optional one).
func Knxnr(args string) (string, error) { return Jrun("knxnr", args) }

// Kdsurf runs kdsurf() of the R package (arguments: x y? w? h? win? nx? ny? edge? at? cover?; ? marks an optional one).
func Kdsurf(args string) (string, error) { return Jrun("kdsurf", args) }

// Dbscn runs dbscn() of the R package (arguments: x y? eps minpts? border? lonlat?; ? marks an optional one).
func Dbscn(args string) (string, error) { return Jrun("dbscn", args) }

// Hsper runs hsper() of the R package (arguments: x y? prd cell top? org? lonlat?; ? marks an optional one).
func Hsper(args string) (string, error) { return Jrun("hsper", args) }

// Ebshr runs ebshr() of the R package (arguments: y v? w? b?; ? marks an optional one).
func Ebshr(args string) (string, error) { return Jrun("ebshr", args) }

// Kscan runs kscan() of the R package (arguments: x y? t r tw minc? top? nsim? seed? lonlat?; ? marks an optional one).
func Kscan(args string) (string, error) { return Jrun("kscan", args) }

// Gmoran runs gmoran() of the R package (arguments: x w nsim? seed?; ? marks an optional one).
func Gmoran(args string) (string, error) { return Jrun("gmoran", args) }

// Bvmor runs bvmor() of the R package (arguments: x y w nsim? seed?; ? marks an optional one).
func Bvmor(args string) (string, error) { return Jrun("bvmor", args) }

// Sarml runs sarml() of the R package (arguments: y x w model? durbin? method? int?; ? marks an optional one).
func Sarml(args string) (string, error) { return Jrun("sarml", args) }

// Spill runs spill() of the R package (arguments: fit? w? rho? b? th?; ? marks an optional one).
func Spill(args string) (string, error) { return Jrun("spill", args) }

// Crossk runs crossk() of the R package (arguments: x y? marks i j win r sigma? nsim? seed?; ? marks an optional one).
func Crossk(args string) (string, error) { return Jrun("crossk", args) }

// Qdrtst runs qdrtst() of the R package (arguments: x y? win nx ny?; ? marks an optional one).
func Qdrtst(args string) (string, error) { return Jrun("qdrtst", args) }

// Morsta runs morsta() of the R package (arguments: counts crit?; ? marks an optional one).
func Morsta(args string) (string, error) { return Jrun("morsta", args) }

// Dissix runs dissix() of the R package (arguments: x delta? w?; ? marks an optional one).
func Dissix(args string) (string, error) { return Jrun("dissix", args) }

// Pinpl runs pinpl() of the R package (arguments: x y? poly pop?; ? marks an optional one).
func Pinpl(args string) (string, error) { return Jrun("pinpl", args) }

// Arwt runs arwt() of the R package (arguments: x from to area? id? top?; ? marks an optional one).
func Arwt(args string) (string, error) { return Jrun("arwt", args) }

// Bufas runs bufas() of the R package (arguments: x y? site comp r buf? lonlat?; ? marks an optional one).
func Bufas(args string) (string, error) { return Jrun("bufas", args) }

// Thrtst runs thrtst() of the R package (arguments: stops srch hits grp area? off? ref? level? maxit? rank? fit? chains? iter? seed?; ? marks an optional one).
func Thrtst(args string) (string, error) { return Jrun("thrtst", args) }

// Ltrlrng runs ltrlrng() of the R package (arguments: x W model? m? navsd?; ? marks an optional one).
func Ltrlrng(args string) (string, error) { return Jrun("ltrlrng", args) }

// Swpwdth runs swpwdth() of the R package (arguments: x p; ? marks an optional one).
func Swpwdth(args string) (string, error) { return Jrun("swpwdth", args) }

// Srchpod runs srchpod() of the R package (arguments: coverage? model? m? W? effort? area?; ? marks an optional one).
func Srchpod(args string) (string, error) { return Jrun("srchpod", args) }

// Pswppod runs pswppod() of the R package (arguments: W S model? m? navsd? npts? ntrk?; ? marks an optional one).
func Pswppod(args string) (string, error) { return Jrun("pswppod", args) }

// Srchscs runs srchscs() of the R package (arguments: poc pod; ? marks an optional one).
func Srchscs(args string) (string, error) { return Jrun("srchscs", args) }

// Datumpmap runs datumpmap() of the R package (arguments: xedg yedg datum? center? sigma? end? nodes? wts?; ? marks an optional one).
func Datumpmap(args string) (string, error) { return Jrun("datumpmap", args) }

// Optcrcsrch runs optcrcsrch() of the R package (arguments: sigma W effort; ? marks an optional one).
func Optcrcsrch(args string) (string, error) { return Jrun("optcrcsrch", args) }

// Opteffort runs opteffort() of the R package (arguments: poc area? W effort tol?; ? marks an optional one).
func Opteffort(args string) (string, error) { return Jrun("opteffort", args) }

// Optsqr runs optsqr() of the R package (arguments: effort model? m? E?; ? marks an optional one).
func Optsqr(args string) (string, error) { return Jrun("optsqr", args) }

// Srchfig runs srchfig() of the R package (arguments: fig? num?; ? marks an optional one).
func Srchfig(args string) (string, error) { return Jrun("srchfig", args) }

// Pgsig runs pgsig() of the R package (arguments: x stab? set?; ? marks an optional one).
func Pgsig(args string) (string, error) { return Jrun("pgsig", args) }

// Briggs runs briggs() of the R package (arguments: x u diam exitvel stktmp ambtmp stab? dthdz? g?; ? marks an optional one).
func Briggs(args string) (string, error) { return Jrun("briggs", args) }

// Gplume runs gplume() of the R package (arguments: q u h rcpt stab? set? mixh? nimg?; ? marks an optional one).
func Gplume(args string) (string, error) { return Jrun("gplume", args) }

// Gspuff runs gspuff() of the R package (arguments: mass u h rcpt t stab? set? mixh? nimg?; ? marks an optional one).
func Gspuff(args string) (string, error) { return Jrun("gspuff", args) }

// Advdif2d runs advdif2d() of the R package (arguments: c0 u v kx ky dx dy dt nsteps source? maxwork?; ? marks an optional one).
func Advdif2d(args string) (string, error) { return Jrun("advdif2d", args) }

// Lagrprt runs lagrprt() of the R package (arguments: n x0 y0 u v kx ky dt nsteps seed? grid? maxwork?; ? marks an optional one).
func Lagrprt(args string) (string, error) { return Jrun("lagrprt", args) }

// Plmfld runs plmfld() of the R package (arguments: q u h x y z? stab? set? mixh? nimg?; ? marks an optional one).
func Plmfld(args string) (string, error) { return Jrun("plmfld", args) }

// Crfpm25 runs crfpm25() of the R package (arguments: expo y? refconc?; ? marks an optional one).
func Crfpm25(args string) (string, error) { return Jrun("crfpm25", args) }

// Crfno2 runs crfno2() of the R package (arguments: expo y? refconc? betaper10?; ? marks an optional one).
func Crfno2(args string) (string, error) { return Jrun("crfno2", args) }

// Attrfrac runs attrfrac() of the R package (arguments: rr expoprev rrlo? rrup?; ? marks an optional one).
func Attrfrac(args string) (string, error) { return Jrun("attrfrac", args) }

// Mortdisp runs mortdisp() of the R package (arguments: expdel pop brate betaperunit rare? hirate?; ? marks an optional one).
func Mortdisp(args string) (string, error) { return Jrun("mortdisp", args) }

// Polbrdn runs polbrdn() of the R package (arguments: expomean expoprev brate pop poll? y? refconc? rare? hirate?; ? marks an optional one).
func Polbrdn(args string) (string, error) { return Jrun("polbrdn", args) }

// Polbrdarea runs polbrdarea() of the R package (arguments: areatab areacol? expcol? popcol? basratecol? poll? y? rare? hirate?; ? marks an optional one).
func Polbrdarea(args string) (string, error) { return Jrun("polbrdarea", args) }

// Expconcidx runs expconcidx() of the R package (arguments: data expo income; ? marks an optional one).
func Expconcidx(args string) (string, error) { return Jrun("expconcidx", args) }

// Onlsmry runs onlsmry() of the R package (arguments: x?; ? marks an optional one).
func Onlsmry(args string) (string, error) { return Jrun("onlsmry", args) }

// Smryupd runs smryupd() of the R package (arguments: acc x; ? marks an optional one).
func Smryupd(args string) (string, error) { return Jrun("smryupd", args) }

// Smrymrg runs smrymrg() of the R package (arguments: a b; ? marks an optional one).
func Smrymrg(args string) (string, error) { return Jrun("smrymrg", args) }

// Smrystat runs smrystat() of the R package (arguments: acc; ? marks an optional one).
func Smrystat(args string) (string, error) { return Jrun("smrystat", args) }

// Rsvidx runs rsvidx() of the R package (arguments: n k seed?; ? marks an optional one).
func Rsvidx(args string) (string, error) { return Jrun("rsvidx", args) }

// Rsvsmp runs rsvsmp() of the R package (arguments: x k seed?; ? marks an optional one).
func Rsvsmp(args string) (string, error) { return Jrun("rsvsmp", args) }

// Dstnskch runs dstnskch() of the R package (arguments: x p? regstr?; ? marks an optional one).
func Dstnskch(args string) (string, error) { return Jrun("dstnskch", args) }

// Dstncnt runs dstncnt() of the R package (arguments: regstr; ? marks an optional one).
func Dstncnt(args string) (string, error) { return Jrun("dstncnt", args) }

// Skchmrg runs skchmrg() of the R package (arguments: a b; ? marks an optional one).
func Skchmrg(args string) (string, error) { return Jrun("skchmrg", args) }

// Co2int runs co2int() of the R package (arguments: loc?; ? marks an optional one).
func Co2int(args string) (string, error) { return Jrun("co2int", args) }

// Co2inttab runs co2inttab() of the R package (no arguments; ? marks an optional one).
func Co2inttab(args string) (string, error) { return Jrun("co2inttab", args) }

// Dtctloc runs dtctloc() of the R package (arguments: envloc? codloc? tz? locale?; ? marks an optional one).
func Dtctloc(args string) (string, error) { return Jrun("dtctloc", args) }

// Raplok runs raplok() of the R package (arguments: rapldir?; ? marks an optional one).
func Raplok(args string) (string, error) { return Jrun("raplok", args) }

// Fpeqv runs fpeqv() of the R package (arguments: co2eg; ? marks an optional one).
func Fpeqv(args string) (string, error) { return Jrun("fpeqv", args) }

// Exprspplr runs exprspplr() of the R package (arguments: data y expo confs nboot? rstate?; ? marks an optional one).
func Exprspplr(args string) (string, error) { return Jrun("exprspplr", args) }

// Inlhist runs inlhist() of the R package (arguments: x bins?; ? marks an optional one).
func Inlhist(args string) (string, error) { return Jrun("inlhist", args) }

// Freqtbl runs freqtbl() of the R package (arguments: data col? maxlvl? sort?; ? marks an optional one).
func Freqtbl(args string) (string, error) { return Jrun("freqtbl", args) }

// Cortbl runs cortbl() of the R package (arguments: data meth? minprs?; ? marks an optional one).
func Cortbl(args string) (string, error) { return Jrun("cortbl", args) }

// Topcors runs topcors() of the R package (arguments: data n? meth? minabs?; ? marks an optional one).
func Topcors(args string) (string, error) { return Jrun("topcors", args) }

// Drpempty runs drpempty() of the R package (arguments: data which? quiet?; ? marks an optional one).
func Drpempty(args string) (string, error) { return Jrun("drpempty", args) }

// Drpconst runs drpconst() of the R package (arguments: data naasval? quiet?; ? marks an optional one).
func Drpconst(args string) (string, error) { return Jrun("drpconst", args) }

// Clncolnms runs clncolnms() of the R package (arguments: data case? sep? quiet?; ? marks an optional one).
func Clncolnms(args string) (string, error) { return Jrun("clncolnms", args) }

// Mahaoutl runs mahaoutl() of the R package (arguments: data alpha? robust?; ? marks an optional one).
func Mahaoutl(args string) (string, error) { return Jrun("mahaoutl", args) }

// Missruns runs missruns() of the R package (arguments: data minrun?; ? marks an optional one).
func Missruns(args string) (string, error) { return Jrun("missruns", args) }

// Missgrid runs missgrid() of the R package (arguments: data height? width?; ? marks an optional one).
func Missgrid(args string) (string, error) { return Jrun("missgrid", args) }

// Mcartst runs mcartst() of the R package (arguments: data maxit? tol?; ? marks an optional one).
func Mcartst(args string) (string, error) { return Jrun("mcartst", args) }

// Misspat runs misspat() of the R package (arguments: data maxpat?; ? marks an optional one).
func Misspat(args string) (string, error) { return Jrun("misspat", args) }

// Duprows runs duprows() of the R package (arguments: data cols? quiet?; ? marks an optional one).
func Duprows(args string) (string, error) { return Jrun("duprows", args) }

// Misssmry runs misssmry() of the R package (arguments: data; ? marks an optional one).
func Misssmry(args string) (string, error) { return Jrun("misssmry", args) }

// Profcols runs profcols() of the R package (arguments: data qnts?; ? marks an optional one).
func Profcols(args string) (string, error) { return Jrun("profcols", args) }

// Dstntst runs dstntst() of the R package (arguments: data grp est lower upper qual? conflvl? pooled?; ? marks an optional one).
func Dstntst(args string) (string, error) { return Jrun("dstntst", args) }

// Dstfit runs dstfit() of the R package (arguments: x fams? nboot? conflvl? seed? k?; ? marks an optional one).
func Dstfit(args string) (string, error) { return Jrun("dstfit", args) }

// Ibench runs ibench() of the R package (arguments: y off x? minn? fdr? dr? local?; ? marks an optional one).
func Ibench(args string) (string, error) { return Jrun("ibench", args) }

// Lnroc runs lnroc() of the R package (arguments: tp pick conf proc? size? far? nboot? seed? level?; ? marks an optional one).
func Lnroc(args string) (string, error) { return Jrun("lnroc", args) }

// Lnsize runs lnsize() of the R package (arguments: choice k? susp? level?; ? marks an optional one).
func Lnsize(args string) (string, error) { return Jrun("lnsize", args) }

// Rptvic runs rptvic() of the R package (arguments: id t? pop? brk?; ? marks an optional one).
func Rptvic(args string) (string, error) { return Jrun("rptvic", args) }

// Wdqx runs wdqx() of the R package (arguments: a b c; ? marks an optional one).
func Wdqx(args string) (string, error) { return Jrun("wdqx", args) }

// Bayout runs bayout() of the R package (arguments: x freq? b? w? acty? alpha? at?; ? marks an optional one).
func Bayout(args string) (string, error) { return Jrun("bayout", args) }

// Bufexp runs bufexp() of the R package (arguments: x y? src r w? lonlat?; ? marks an optional one).
func Bufexp(args string) (string, error) { return Jrun("bufexp", args) }

// Kerexp runs kerexp() of the R package (arguments: x y? src h w? kern?; ? marks an optional one).
func Kerexp(args string) (string, error) { return Jrun("kerexp", args) }

// Lerprc runs lerprc() of the R package (arguments: w rho tau?; ? marks an optional one).
func Lerprc(args string) (string, error) { return Jrun("lerprc", args) }

// Bym2sc runs bym2sc() of the R package (arguments: w tau? phi? tol?; ? marks an optional one).
func Bym2sc(args string) (string, error) { return Jrun("bym2sc", args) }

// Ecolrg runs ecolrg() of the R package (arguments: y e x? family? z? tol? maxit?; ? marks an optional one).
func Ecolrg(args string) (string, error) { return Jrun("ecolrg", args) }

// Spbias runs spbias() of the R package (arguments: d b ub?; ? marks an optional one).
func Spbias(args string) (string, error) { return Jrun("spbias", args) }

// Spsir runs spsir() of the R package (arguments: i0 w? pop? beta? gamma? kappa? tmax? dt?; ? marks an optional one).
func Spsir(args string) (string, error) { return Jrun("spsir", args) }

// Sfrail runs sfrail() of the R package (arguments: time event x? reg w tau? tol?; ? marks an optional one).
func Sfrail(args string) (string, error) { return Jrun("sfrail", args) }

// Sfcure runs sfcure() of the R package (arguments: time event x? reg w tau? tol?; ? marks an optional one).
func Sfcure(args string) (string, error) { return Jrun("sfcure", args) }

// Evfilt runs evfilt() of the R package (arguments: v w x? method? cand? tol? maxv?; ? marks an optional one).
func Evfilt(args string) (string, error) { return Jrun("evfilt", args) }

// Vgemp runs vgemp() of the R package (arguments: x y? v nbins? cutoff?; ? marks an optional one).
func Vgemp(args string) (string, error) { return Jrun("vgemp", args) }

// Vgfit runs vgfit() of the R package (arguments: x y? v model?; ? marks an optional one).
func Vgfit(args string) (string, error) { return Jrun("vgfit", args) }

// Okrge runs okrge() of the R package (arguments: x y? v at vgm?; ? marks an optional one).
func Okrge(args string) (string, error) { return Jrun("okrge", args) }

// Bkrig runs bkrig() of the R package (arguments: x y? v blocks nq? nugget? sill? range?; ? marks an optional one).
func Bkrig(args string) (string, error) { return Jrun("bkrig", args) }

// Krigcv runs krigcv() of the R package (arguments: x y? v nugget? sill? range?; ? marks an optional one).
func Krigcv(args string) (string, error) { return Jrun("krigcv", args) }

// Trsurf runs trsurf() of the R package (arguments: x y? v order?; ? marks an optional one).
func Trsurf(args string) (string, error) { return Jrun("trsurf", args) }

// Sglmx runs sglmx() of the R package (arguments: v z? x y?; ? marks an optional one).
func Sglmx(args string) (string, error) { return Jrun("sglmx", args) }

// Smixm runs smixm() of the R package (arguments: v z? x y?; ? marks an optional one).
func Smixm(args string) (string, error) { return Jrun("smixm", args) }

// Spobj runs spobj() of the R package (arguments: data coords? crs?; ? marks an optional one).
func Spobj(args string) (string, error) { return Jrun("spobj", args) }

// Spxy runs spxy() of the R package (arguments: x; ? marks an optional one).
func Spxy(args string) (string, error) { return Jrun("spxy", args) }

// Crsof runs crsof() of the R package (arguments: x; ? marks an optional one).
func Crsof(args string) (string, error) { return Jrun("crsof", args) }

// Lfspec runs lfspec() of the R package (arguments: x; ? marks an optional one).
func Lfspec(args string) (string, error) { return Jrun("lfspec", args) }

// Spsamp runs spsamp() of the R package (arguments: x size dens? seed?; ? marks an optional one).
func Spsamp(args string) (string, error) { return Jrun("spsamp", args) }

// Prohs runs prohs() of the R package (arguments: x y? t at win cell bw? tw? tunit? top? ahead?; ? marks an optional one).
func Prohs(args string) (string, error) { return Jrun("prohs", args) }

// Lnkpr runs lnkpr() of the R package (arguments: x y? t series mo?; ? marks an optional one).
func Lnkpr(args string) (string, error) { return Jrun("lnkpr", args) }

// Spcure runs spcure() of the R package (arguments: time event x? dist? level? maxit?; ? marks an optional one).
func Spcure(args string) (string, error) { return Jrun("spcure", args) }

// Apcie runs apcie() of the R package (arguments: y expo level? maxit?; ? marks an optional one).
func Apcie(args string) (string, error) { return Jrun("apcie", args) }

// Hsacc runs hsacc() of the R package (arguments: pred obs area? top? prev?; ? marks an optional one).
func Hsacc(args string) (string, error) { return Jrun("hsacc", args) }

// Sctidy runs sctidy() of the R package (arguments: x; ? marks an optional one).
func Sctidy(args string) (string, error) { return Jrun("sctidy", args) }

// Qltyrcd runs qltyrcd() of the R package (arguments: x suppr? caution?; ? marks an optional one).
func Qltyrcd(args string) (string, error) { return Jrun("qltyrcd", args) }

// Limfrom runs limfrom() of the R package (arguments: est se? cv? moe? moelev? conflvl?; ? marks an optional one).
func Limfrom(args string) (string, error) { return Jrun("limfrom", args) }

// Indgap runs indgap() of the R package (arguments: data grp est lower upper ref? qual? conflvl? pooled?; ? marks an optional one).
func Indgap(args string) (string, error) { return Jrun("indgap", args) }

// Genexpo runs genexpo() of the R package (arguments: parent gparent partial?; ? marks an optional one).
func Genexpo(args string) (string, error) { return Jrun("genexpo", args) }

// Exptrnd runs exptrnd() of the R package (arguments: cases n dose? conflvl?; ? marks an optional one).
func Exptrnd(args string) (string, error) { return Jrun("exptrnd", args) }

// Trccmpl runs trccmpl() of the R package (no arguments; ? marks an optional one).
func Trccmpl(args string) (string, error) { return Jrun("trccmpl", args) }

// Trcprog runs trcprog() of the R package (arguments: done? yearend? total? start? conflvl?; ? marks an optional one).
func Trcprog(args string) (string, error) { return Jrun("trcprog", args) }

// Mhcont runs mhcont() of the R package (arguments: data grp ind est lower upper hope? belong? meaning? purp? qual? conflvl?; ? marks an optional one).
func Mhcont(args string) (string, error) { return Jrun("mhcont", args) }

// Cwbscr runs cwbscr() of the R package (arguments: incpercap highsch20p uni25p notcro nomajrep partic empl bounds?; ? marks an optional one).
func Cwbscr(args string) (string, error) { return Jrun("cwbscr", args) }

// Connassoc runs connassoc() of the R package (arguments: data place x y xlo? xup? ylo? yup? reps? seed? conflvl?; ? marks an optional one).
func Connassoc(args string) (string, error) { return Jrun("connassoc", args) }

// Healchg runs healchg() of the R package (arguments: data grp ind prd est lower upper from to well? qual? conflvl? distress?; ? marks an optional one).
func Healchg(args string) (string, error) { return Jrun("healchg", args) }

// Commidx runs commidx() of the R package (arguments: data grp ind est lower upper w qual? conflvl?; ? marks an optional one).
func Commidx(args string) (string, error) { return Jrun("commidx", args) }

// Frmaud runs frmaud() of the R package (arguments: x; ? marks an optional one).
func Frmaud(args string) (string, error) { return Jrun("frmaud", args) }

// Scests runs scests() of the R package (arguments: x; ? marks an optional one).
func Scests(args string) (string, error) { return Jrun("scests", args) }

// Ordcmp runs ordcmp() of the R package (arguments: data grp cat est lower upper lvls ref qual? conflvl? pooled?; ? marks an optional one).
func Ordcmp(args string) (string, error) { return Jrun("ordcmp", args) }

// Stdprev runs stdprev() of the R package (arguments: est lower upper std grp? qual? conflvl? scale?; ? marks an optional one).
func Stdprev(args string) (string, error) { return Jrun("stdprev", args) }

// Welprof runs welprof() of the R package (arguments: data grp ind est lower upper well distress domain? qual?; ? marks an optional one).
func Welprof(args string) (string, error) { return Jrun("welprof", args) }

// Holbal runs holbal() of the R package (arguments: data grp domain est lower upper w? qual? reps? seed? conflvl?; ? marks an optional one).
func Holbal(args string) (string, error) { return Jrun("holbal", args) }

// Dualcont runs dualcont() of the R package (arguments: x well distress conflvl?; ? marks an optional one).
func Dualcont(args string) (string, error) { return Jrun("dualcont", args) }

// Cnsprof runs cnsprof() of the R package (arguments: x; ? marks an optional one).
func Cnsprof(args string) (string, error) { return Jrun("cnsprof", args) }

// Profshr runs profshr() of the R package (arguments: x num den? smpfrac? deff? conflvl? roubase?; ? marks an optional one).
func Profshr(args string) (string, error) { return Jrun("profshr", args) }

// Cnsstd runs cnsstd() of the R package (arguments: x grp? geo? ages? resid?; ? marks an optional one).
func Cnsstd(args string) (string, error) { return Jrun("cnsstd", args) }

// Cnscwb runs cnscwb() of the R package (arguments: x dwell? groups? bounds? resid? minpop?; ? marks an optional one).
func Cnscwb(args string) (string, error) { return Jrun("cnscwb", args) }

// Trcind runs trcind() of the R package (no arguments; ? marks an optional one).
func Trcind(args string) (string, error) { return Jrun("trcind", args) }

// Trchlthgap runs trchlthgap() of the R package (arguments: health? lifeexp? infmort? suiide? ref? agetot?; ? marks an optional one).
func Trchlthgap(args string) (string, error) { return Jrun("trchlthgap", args) }

// Polex runs polex() of the R package (arguments: seed?; ? marks an optional one).
func Polex(args string) (string, error) { return Jrun("polex", args) }

// Votagr runs votagr() of the R package (arguments: x; ? marks an optional one).
func Votagr(args string) (string, error) { return Jrun("votagr", args) }

// Amscal runs amscal() of the R package (arguments: Z; ? marks an optional one).
func Amscal(args string) (string, error) { return Jrun("amscal", args) }

// Bbscal runs bbscal() of the R package (arguments: X ndim? minscale?; ? marks an optional one).
func Bbscal(args string) (string, error) { return Jrun("bbscal", args) }

// Optcls runs optcls() of the R package (arguments: votes ndim? maxit? nstart? seed?; ? marks an optional one).
func Optcls(args string) (string, error) { return Jrun("optcls", args) }

// Dblctr runs dblctr() of the R package (arguments: D; ? marks an optional one).
func Dblctr(args string) (string, error) { return Jrun("dblctr", args) }

// Cmdsf runs cmdsf() of the R package (arguments: D ndim?; ? marks an optional one).
func Cmdsf(args string) (string, error) { return Jrun("cmdsf", args) }

// Smacf runs smacf() of the R package (arguments: D ndim? maxit? tol? w? init?; ? marks an optional one).
func Smacf(args string) (string, error) { return Jrun("smacf", args) }

// Nmdsf runs nmdsf() of the R package (arguments: D ndim? maxit? tol?; ? marks an optional one).
func Nmdsf(args string) (string, error) { return Jrun("nmdsf", args) }

// Mdsfit runs mdsfit() of the R package (arguments: ev; ? marks an optional one).
func Mdsfit(args string) (string, error) { return Jrun("mdsfit", args) }

// Unfstr runs unfstr() of the R package (arguments: xr xs D w?; ? marks an optional one).
func Unfstr(args string) (string, error) { return Jrun("unfstr", args) }

// Mlsmu runs mlsmu() of the R package (arguments: D ndim? maxit? tol? nstart?; ? marks an optional one).
func Mlsmu(args string) (string, error) { return Jrun("mlsmu", args) }

// Smunf runs smunf() of the R package (arguments: D ndim? maxit? tol?; ? marks an optional one).
func Smunf(args string) (string, error) { return Jrun("smunf", args) }

// Ipsrec runs ipsrec() of the R package (arguments: xr xs?; ? marks an optional one).
func Ipsrec(args string) (string, error) { return Jrun("ipsrec", args) }

// Nomutl runs nomutl() of the R package (arguments: x zy zn beta? w?; ? marks an optional one).
func Nomutl(args string) (string, error) { return Jrun("nomutl", args) }

// Nompr runs nompr() of the R package (arguments: x zy zn beta? w?; ? marks an optional one).
func Nompr(args string) (string, error) { return Jrun("nompr", args) }

// Nomll runs nomll() of the R package (arguments: votes x zy zn beta? w?; ? marks an optional one).
func Nomll(args string) (string, error) { return Jrun("nomll", args) }

// Procrs runs procrs() of the R package (arguments: X target; ? marks an optional one).
func Procrs(args string) (string, error) { return Jrun("procrs", args) }

// Irtlik runs irtlik() of the R package (arguments: votes x alpha beta; ? marks an optional one).
func Irtlik(args string) (string, error) { return Jrun("irtlik", args) }

// Irtpst runs irtpst() of the R package (arguments: chain std?; ? marks an optional one).
func Irtpst(args string) (string, error) { return Jrun("irtpst", args) }

// Ordocl runs ordocl() of the R package (arguments: Y ndim? maxit?; ? marks an optional one).
func Ordocl(args string) (string, error) { return Jrun("ordocl", args) }

// Anchv runs anchv() of the R package (arguments: y V ncat?; ? marks an optional one).
func Anchv(args string) (string, error) { return Jrun("anchv", args) }

// Indscl runs indscl() of the R package (arguments: D ndim? maxit? tol?; ? marks an optional one).
func Indscl(args string) (string, error) { return Jrun("indscl", args) }

// Nrmvec runs nrmvec() of the R package (arguments: x y; ? marks an optional one).
func Nrmvec(args string) (string, error) { return Jrun("nrmvec", args) }

// Cutlns runs cutlns() of the R package (arguments: normals cuts xlim?; ? marks an optional one).
func Cutlns(args string) (string, error) { return Jrun("cutlns", args) }

// Dwnom runs dwnom() of the R package (arguments: votes ndim? maxit?; ? marks an optional one).
func Dwnom(args string) (string, error) { return Jrun("dwnom", args) }

// Nombt runs nombt() of the R package (arguments: votes x normals cuts nboot? seed?; ? marks an optional one).
func Nombt(args string) (string, error) { return Jrun("nombt", args) }

// Irtem runs irtem() of the R package (arguments: votes ndim? maxit? tol? xprior? bprior? start? conv?; ? marks an optional one).
func Irtem(args string) (string, error) { return Jrun("irtem", args) }

// Scboot runs scboot() of the R package (arguments: Z method? nboot? seed?; ? marks an optional one).
func Scboot(args string) (string, error) { return Jrun("scboot", args) }

// Wrdfsh runs wrdfsh() of the R package (arguments: dtm maxit? tol?; ? marks an optional one).
func Wrdfsh(args string) (string, error) { return Jrun("wrdfsh", args) }

// Bamscl runs bamscl() of the R package (arguments: Z nsamp? burn? polarity? seed?; ? marks an optional one).
func Bamscl(args string) (string, error) { return Jrun("bamscl", args) }

// Bmdsf runs bmdsf() of the R package (arguments: D ndim? nsamp? burn? sigma? seed?; ? marks an optional one).
func Bmdsf(args string) (string, error) { return Jrun("bmdsf", args) }

// Bunfld runs bunfld() of the R package (arguments: D ndim? nsamp? burn? seed?; ? marks an optional one).
func Bunfld(args string) (string, error) { return Jrun("bunfld", args) }

// Cjrirt runs cjrirt() of the R package (arguments: votes ndim? nsamp? burn? bvar? start? seed? chains?; ? marks an optional one).
func Cjrirt(args string) (string, error) { return Jrun("cjrirt", args) }

// Anomin runs anomin() of the R package (arguments: votes ndim? nsamp? burn? seed? thin? lop? minvotes? polarity? fixalpha?; ? marks an optional one).
func Anomin(args string) (string, error) { return Jrun("anomin", args) }

// Ordirt runs ordirt() of the R package (arguments: Y ndim? nsamp? burn? seed? lprec?; ? marks an optional one).
func Ordirt(args string) (string, error) { return Jrun("ordirt", args) }

// Dynirt runs dynirt() of the R package (arguments: votes period nsamp? burn? seed? thin? tau2? c0? d0? e0? E0? a0? A0? b0? B0? anchor?; ? marks an optional one).
func Dynirt(args string) (string, error) { return Jrun("dynirt", args) }

// Arimfc runs arimfc() of the R package (arguments: y h?; ? marks an optional one).
func Arimfc(args string) (string, error) { return Jrun("arimfc", args) }

// Sarfc runs sarfc() of the R package (arguments: y h? order? seasonal?; ? marks an optional one).
func Sarfc(args string) (string, error) { return Jrun("sarfc", args) }

// Catcor runs catcor() of the R package (arguments: units; ? marks an optional one).
func Catcor(args string) (string, error) { return Jrun("catcor", args) }

// Ousim runs ousim() of the R package (arguments: y npath? days? dt? seed?; ? marks an optional one).
func Ousim(args string) (string, error) { return Jrun("ousim", args) }

// Fpdens runs fpdens() of the R package (arguments: y ngrid? nstep?; ? marks an optional one).
func Fpdens(args string) (string, error) { return Jrun("fpdens", args) }

// Lvfit runs lvfit() of the R package (arguments: y year?; ? marks an optional one).
func Lvfit(args string) (string, error) { return Jrun("lvfit", args) }

// Mtrbnd runs mtrbnd() of the R package (arguments: y z ylo? yhi? assume? w?; ? marks an optional one).
func Mtrbnd(args string) (string, error) { return Jrun("mtrbnd", args) }

// Yoywr runs yoywr() of the R package (arguments: x file format? digits? na? meta? comment? pretty? align?; ? marks an optional one).
func Yoywr(args string) (string, error) { return Jrun("yoywr", args) }

// Yoymd runs yoymd() of the R package (arguments: x file? digits? align? meta?; ? marks an optional one).
func Yoymd(args string) (string, error) { return Jrun("yoymd", args) }

// Yoyjs runs yoyjs() of the R package (arguments: x file? digits? pretty? meta?; ? marks an optional one).
func Yoyjs(args string) (string, error) { return Jrun("yoyjs", args) }

// Yoytsv runs yoytsv() of the R package (arguments: x file? digits? na? meta? comment?; ? marks an optional one).
func Yoytsv(args string) (string, error) { return Jrun("yoytsv", args) }

// Yoycsv runs yoycsv() of the R package (arguments: x file? digits? na? meta? comment?; ? marks an optional one).
func Yoycsv(args string) (string, error) { return Jrun("yoycsv", args) }

// Yoysmry runs yoysmry() of the R package (arguments: x lag?; ? marks an optional one).
func Yoysmry(args string) (string, error) { return Jrun("yoysmry", args) }

// Yoylbl runs yoylbl() of the R package (arguments: x labels; ? marks an optional one).
func Yoylbl(args string) (string, error) { return Jrun("yoylbl", args) }

// Fylbl runs fylbl() of the R package (arguments: endyear sep? short?; ? marks an optional one).
func Fylbl(args string) (string, error) { return Jrun("fylbl", args) }

// Drfthomog runs drfthomog() of the R package (arguments: x y seed?; ? marks an optional one).
func Drfthomog(args string) (string, error) { return Jrun("drfthomog", args) }

// Drftchisq runs drftchisq() of the R package (arguments: obs expd; ? marks an optional one).
func Drftchisq(args string) (string, error) { return Jrun("drftchisq", args) }

// Drftpsi runs drftpsi() of the R package (arguments: x y bins? eps?; ? marks an optional one).
func Drftpsi(args string) (string, error) { return Jrun("drftpsi", args) }

// Drftks runs drftks() of the R package (arguments: x y; ? marks an optional one).
func Drftks(args string) (string, error) { return Jrun("drftks", args) }

// Rgnmapalt runs rgnmapalt() of the R package (arguments: map unit reg route known?; ? marks an optional one).
func Rgnmapalt(args string) (string, error) { return Jrun("rgnmapalt", args) }

// Rgnmapcmp runs rgnmapcmp() of the R package (arguments: pub obs unit cols?; ? marks an optional one).
func Rgnmapcmp(args string) (string, error) { return Jrun("rgnmapcmp", args) }

// Rgnmapintg runs rgnmapintg() of the R package (arguments: map unit reg regs?; ? marks an optional one).
func Rgnmapintg(args string) (string, error) { return Jrun("rgnmapintg", args) }

// Xfrvrf runs xfrvrf() of the R package (arguments: imported soucnts codebook? tol? strict?; ? marks an optional one).
func Xfrvrf(args string) (string, error) { return Jrun("xfrvrf", args) }

// Vrfmrg runs vrfmrg() of the R package (arguments: x pub tol? strict?; ? marks an optional one).
func Vrfmrg(args string) (string, error) { return Jrun("vrfmrg", args) }

// Vrfrcd runs vrfrcd() of the R package (arguments: orig recd decl; ? marks an optional one).
func Vrfrcd(args string) (string, error) { return Jrun("vrfrcd", args) }

// Grdbin runs grdbin() of the R package (arguments: x col; ? marks an optional one).
func Grdbin(args string) (string, error) { return Jrun("grdbin", args) }

// Declbld runs declbld() of the R package (arguments: x vlabs?; ? marks an optional one).
func Declbld(args string) (string, error) { return Jrun("declbld", args) }

// Grdlvls runs grdlvls() of the R package (arguments: x lvls ref?; ? marks an optional one).
func Grdlvls(args string) (string, error) { return Jrun("grdlvls", args) }

// Relbl runs relbl() of the R package (arguments: x map keep?; ? marks an optional one).
func Relbl(args string) (string, error) { return Jrun("relbl", args) }

// Grdrcd runs grdrcd() of the R package (arguments: x map keep?; ? marks an optional one).
func Grdrcd(args string) (string, error) { return Jrun("grdrcd", args) }

// Deccodes runs deccodes() of the R package (arguments: x vlabs keepna?; ? marks an optional one).
func Deccodes(args string) (string, error) { return Jrun("deccodes", args) }

// Shrs runs shrs() of the R package (arguments: x count? by? total? conflvl?; ? marks an optional one).
func Shrs(args string) (string, error) { return Jrun("shrs", args) }

// Wrtxtfb runs wrtxtfb() of the R package (arguments: text path; ? marks an optional one).
func Wrtxtfb(args string) (string, error) { return Jrun("wrtxtfb", args) }

// Asciifb runs asciifb() of the R package (arguments: x force?; ? marks an optional one).
func Asciifb(args string) (string, error) { return Jrun("asciifb", args) }

// Toascii runs toascii() of the R package (arguments: x; ? marks an optional one).
func Toascii(args string) (string, error) { return Jrun("toascii", args) }

// Cunif runs cunif() of the R package (arguments: n seed; ? marks an optional one).
func Cunif(args string) (string, error) { return Jrun("cunif", args) }

// Cgamcdf runs cgamcdf() of the R package (arguments: shape x; ? marks an optional one).
func Cgamcdf(args string) (string, error) { return Jrun("cgamcdf", args) }

// Cipw runs cipw() of the R package (arguments: treat ps trimlo? trimhi?; ? marks an optional one).
func Cipw(args string) (string, error) { return Jrun("cipw", args) }

// Cbootmean runs cbootmean() of the R package (arguments: x B? seed?; ? marks an optional one).
func Cbootmean(args string) (string, error) { return Jrun("cbootmean", args) }

// Ccovm runs ccovm() of the R package (arguments: x; ? marks an optional one).
func Ccovm(args string) (string, error) { return Jrun("ccovm", args) }

// Cmidrk runs cmidrk() of the R package (arguments: x; ? marks an optional one).
func Cmidrk(args string) (string, error) { return Jrun("cmidrk", args) }

// Cspr runs cspr() of the R package (arguments: x y; ? marks an optional one).
func Cspr(args string) (string, error) { return Jrun("cspr", args) }

// Cwtd runs cwtd() of the R package (arguments: x w; ? marks an optional one).
func Cwtd(args string) (string, error) { return Jrun("cwtd", args) }

// Cwinsmean runs cwinsmean() of the R package (arguments: x trim?; ? marks an optional one).
func Cwinsmean(args string) (string, error) { return Jrun("cwinsmean", args) }

// Ctrmmean runs ctrmmean() of the R package (arguments: x trim?; ? marks an optional one).
func Ctrmmean(args string) (string, error) { return Jrun("ctrmmean", args) }

// Ctkyfnc runs ctkyfnc() of the R package (arguments: x k?; ? marks an optional one).
func Ctkyfnc(args string) (string, error) { return Jrun("ctkyfnc", args) }

// Ciqr runs ciqr() of the R package (arguments: x; ? marks an optional one).
func Ciqr(args string) (string, error) { return Jrun("ciqr", args) }

// Cmad runs cmad() of the R package (arguments: x const?; ? marks an optional one).
func Cmad(args string) (string, error) { return Jrun("cmad", args) }

// Cmed runs cmed() of the R package (arguments: x; ? marks an optional one).
func Cmed(args string) (string, error) { return Jrun("cmed", args) }

// Cqnt runs cqnt() of the R package (arguments: x probs?; ? marks an optional one).
func Cqnt(args string) (string, error) { return Jrun("cqnt", args) }

// Cmom runs cmom() of the R package (arguments: x; ? marks an optional one).
func Cmom(args string) (string, error) { return Jrun("cmom", args) }

// Cnormlpdf runs cnormlpdf() of the R package (arguments: x mean? sd?; ? marks an optional one).
func Cnormlpdf(args string) (string, error) { return Jrun("cnormlpdf", args) }

// Cnormpdf runs cnormpdf() of the R package (arguments: x mean? sd?; ? marks an optional one).
func Cnormpdf(args string) (string, error) { return Jrun("cnormpdf", args) }

// Cdistm runs cdistm() of the R package (arguments: a b; ? marks an optional one).
func Cdistm(args string) (string, error) { return Jrun("cdistm", args) }

// Csdev runs csdev() of the R package (arguments: x ddof?; ? marks an optional one).
func Csdev(args string) (string, error) { return Jrun("csdev", args) }

// Ccorr runs ccorr() of the R package (arguments: x y; ? marks an optional one).
func Ccorr(args string) (string, error) { return Jrun("ccorr", args) }

// Cvari runs cvari() of the R package (arguments: x; ? marks an optional one).
func Cvari(args string) (string, error) { return Jrun("cvari", args) }

// Cmn runs cmn() of the R package (arguments: x; ? marks an optional one).
func Cmn(args string) (string, error) { return Jrun("cmn", args) }

// Clprep runs clprep() of the R package (arguments: data? preset?; ? marks an optional one).
func Clprep(args string) (string, error) { return Jrun("clprep", args) }

// Pedlr runs pedlr() of the R package (arguments: ped geno freq alt? mut? rate? range?; ? marks an optional one).
func Pedlr(args string) (string, error) { return Jrun("pedlr", args) }

// Gmixlr runs gmixlr() of the R package (arguments: peaks freq ref? hp? nc hd? at? pc? lambda? theta? stut? deg? par?; ? marks an optional one).
func Gmixlr(args string) (string, error) { return Jrun("gmixlr", args) }

// Cmpsc runs cmpsc() of the R package (arguments: a b len? npeaks? tx? full?; ? marks an optional one).
func Cmpsc(args string) (string, error) { return Jrun("cmpsc", args) }

// Ppmlfe runs ppmlfe() of the R package (arguments: data y x fe expo? clus? se? ssc? time? lag? nboot? seed? sep? level? tol? maxit? twc? psd?; ? marks an optional one).
func Ppmlfe(args string) (string, error) { return Jrun("ppmlfe", args) }

// Ivfe runs ivfe() of the R package (arguments: data y endo? inst? x? fe clus? se? ssc? time? lag? nboot? seed? level? tol? twc? psd?; ? marks an optional one).
func Ivfe(args string) (string, error) { return Jrun("ivfe", args) }

// Ppmlcf runs ppmlcf() of the R package (arguments: data y endo inst x? fe expo? clus? how? se? ssc? nboot? seed? level? tol? maxit?; ? marks an optional one).
func Ppmlcf(args string) (string, error) { return Jrun("ppmlcf", args) }

// Windiv runs windiv() of the R package (arguments: dir width? start? area? base? sect? edges? by? avg? share? spd?; ? marks an optional one).
func Windiv(args string) (string, error) { return Jrun("windiv", args) }

// Invsn runs invsn() of the R package (arguments: tlow tup sign? brk?; ? marks an optional one).
func Invsn(args string) (string, error) { return Jrun("invsn", args) }

// Updwn runs updwn() of the R package (arguments: data y ax ay day wind src area? tol? band? dist? twin? agg? norm? wx? se? ssc? clus? lag? nboot? seed? offs? level?; ? marks an optional one).
func Updwn(args string) (string, error) { return Jrun("updwn", args) }

// Expjn runs expjn() of the R package (arguments: expo outc by?; ? marks an optional one).
func Expjn(args string) (string, error) { return Jrun("expjn", args) }

// Smkiv runs smkiv() of the R package (arguments: data plume unit day pm? adj? nbr? z? seas? dens? cf? by?; ? marks an optional one).
func Smkiv(args string) (string, error) { return Jrun("smkiv", args) }

// Exprun runs exprun() of the R package (arguments: data var unit day thr? cap? step?; ? marks an optional one).
func Exprun(args string) (string, error) { return Jrun("exprun", args) }

// Selftest runs selftest() of the R package (no arguments; ? marks an optional one).
func Selftest(args string) (string, error) { return Jrun("selftest", args) }

// Proof runs proof() of the R package (arguments: fn?; ? marks an optional one).
func Proof(args string) (string, error) { return Jrun("proof", args) }
