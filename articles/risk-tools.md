# Risk tools and predictive policing

Fairness of a risk score across groups, its calibration, and audits of a
predictive-policing rule.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Differences between groups in a risk tool’s decisions

[`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md):
How the share of cases a tool flags (predicts positive) differs between
groups, as a difference and as a ratio against a base group, each with
an interval; with the outcomes, how its true and false positive rates
differ; and the bias amplification score, the largest parity gap times
the Gini coefficient of the groups’ shares.

``` r

set.seed(1)
g <- rep(c("a", "b"), each = 200)
y <- rbinom(400, 1, 0.3)
p <- rbinom(400, 1, ifelse(g == "a", 0.4, 0.25))
r <- fairdf(p, g, y, nboot = 199, seed = 1)
r$tab[, c("group", "ppos", "dp", "dplo", "dphi", "di", "dilo", "dihi")]
#>   group  ppos     dp       dplo        dphi        di      dilo      dihi
#> 1     a 0.435  0.000        NaN         NaN 1.0000000       NaN       NaN
#> 2     b 0.260 -0.175 -0.2640056 -0.08182849 0.5977011 0.4507581 0.7925463
r$sum
#>           di           dp         gini          bas       ginilo       ginihi 
#>  0.597701149 -0.175000000  0.125899281 -0.022032374  0.059494196  0.194227687 
#>        baslo        bashi          aod           eo 
#> -0.054404735 -0.004745647 -0.159044465 -0.202023149
```

## Disparity metrics over time and place

[`tmpaud()`](https://rootcoder007.github.io/crimlab/reference/tmpaud.md):
The disparate impact ratio, the parity gap, the Gini coefficient of the
groups’ shares flagged and the bias amplification score of
[`fairdf()`](https://rootcoder007.github.io/crimlab/reference/fairdf.md),
for every city and period, and how far each city’s ratio moves over the
periods (after Barman and Barman 2026).

``` r

set.seed(2)
n <- 600
g <- sample(c("a", "b"), n, TRUE)
per <- sample(c("2024-1", "2024-2", "2024-3"), n, TRUE)
cty <- sample(c("north", "south"), n, TRUE)
pr <- c(south = 0.2, north = 0.3)[cty]
pr[g == "a"] <- 0.4
p <- rbinom(n, 1, pr)
tmpaud(p, g, per, cty)$cities
#>    city ncell        di          dp       gini          bas     dimin     dimax
#> 1 north     3 0.7754133 -0.09992407 0.06547922 -0.008588778 0.6506747 0.9232143
#> 2 south     3 0.3603946 -0.28545866 0.23780500 -0.069788958 0.2545455 0.4545455
#>       dirng ngt1
#> 1 0.2725396    0
#> 2 0.2000000    0
```

## Calibration audit of a predictive-policing ranking

[`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md):
Does the order in which an algorithm ranks areas by predicted risk match
the order of what was later recorded there, and does the disagreement
follow the areas’ demographic make-up?
[`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md)
ranks the areas both ways, measures the agreement by Spearman’s
correlation, and gives each group the mean gap between its areas’ two
ranks, with a permutation test. It generalises the district-level
analysis of the SciencesPo predictive-policing study of Chicago
(Lacherade et al. 2021) to any place: the five columns of the audit are
`area`, `risk`, `outcome`, `group` and `population`.

``` r

set.seed(3)
d <- fairsim(4000, nareas = 24, bias = 0.5, seed = 1)
# the recorded outcome is `det`: the incidents that were detected
r <- calaud(d$area, d$risk, d$det, d$group, nperm = 499, seed = 1)
r$groups
#>   group  n gap     p
#> 1     A 12 -12 0.002
#> 2     B 12  12 0.002
r$spear
#>           rho             s           pex             t            pt 
#> -7.400000e-01  4.002000e+03  5.640811e-05 -5.160381e+00  3.575934e-05
```

## A risk score’s spread across groups

[`scrdsp()`](https://rootcoder007.github.io/crimlab/reference/scrdsp.md):
Each group’s number of cases, mean, median and standard deviation of a
score, each mean’s gap to a reference group with its interval, and the
one-way analysis of variance of the score on group, with Welch’s (1951)
test, which does not assume equal variances.

``` r

set.seed(7)
g <- rep(c("x", "y", "z"), c(40, 40, 20))
s <- rnorm(100, c(x = 0, y = 0.5, z = 1)[g])
scrdsp(s, g)
#> $tab
#>   group  n      mean    median        sd       gap        se          lo        hi
#> 1     x 40 0.2708131 0.1684252 1.0461883 0.0000000       NaN         NaN       NaN
#> 2     y 40 0.6225100 0.6625734 0.9275965 0.3516969 0.2144616 -0.07395004 0.7773438
#> 3     z 20 0.9068368 0.7723597 0.8266946 0.6360236 0.2626607  0.11471478 1.1573325
#> 
#> $anova
#>           f         df1         df2           p 
#>  3.18629701  2.00000000 97.00000000  0.04568425 
#> 
#> $welch
#>           f         df1         df2           p 
#>  3.31086371  2.00000000 54.14739082  0.04400091 
#> 
#> $spread
#> [1] 0.6360236
#> 
#> $ref
#> [1] "x"
```

## Individual fairness: a Lipschitz audit

[`lipaud()`](https://rootcoder007.github.io/crimlab/reference/lipaud.md):
Dwork et al. (2012) ask that similar individuals be treated similarly:
that the gap between two people’s scores be at most `lip` times the
distance between them, for a distance that says who is alike for the
task.
[`lipaud()`](https://rootcoder007.github.io/crimlab/reference/lipaud.md)
checks the given pairs against that rule and reports the smallest
constant the scores satisfy on them.

``` r

x <- c(1, 2, 2.5, 6)
h <- c(0.10, 0.20, 0.60, 0.65)
# scores on a scale of 0 to 1 against distances in years
lipaud(h, dist(x), lip = 0.2)
#> $lreq
#> [1] 0.8
#> 
#> $nviol
#> [1] 2
#> 
#> $rate
#> [1] 0.3333333
#> 
#> $ok
#> [1] FALSE
#> 
#> $maxgap
#> [1] 0.55
#> 
#> $worst
#> [1] 2 3
#> 
#> $pairs
#>   i j   d  gap      ratio  viol
#> 1 1 2 1.0 0.10 0.10000000 FALSE
#> 2 1 3 1.5 0.50 0.33333333  TRUE
#> 3 1 4 5.0 0.55 0.11000000 FALSE
#> 4 2 3 0.5 0.40 0.80000000  TRUE
#> 5 2 4 4.0 0.45 0.11250000 FALSE
#> 6 3 4 3.5 0.05 0.01428571 FALSE
```

## A true base rate behind a noisily recorded one

[`nsbnd()`](https://rootcoder007.github.io/crimlab/reference/nsbnd.md):
A recorded outcome is often a proxy: a rearrest for a reoffence, a
recorded crime for one committed. With false-positive noise `a`, the
chance a case without the outcome is recorded with it, and
false-negative noise `b`, the chance a case with it is not, a true base
rate `p` is recorded as pobs = p (1 - b) + (1 - p) a. Given only the
boxes the noise rates lie in,
[`nsbnd()`](https://rootcoder007.github.io/crimlab/reference/nsbnd.md)
gives the interval of true rates consistent with each recorded rate, and
says for each pair of groups whether the order of their true rates
follows from the records whatever the noise.

``` r

# rearrest rates of two groups; up to 10% false positives, 20% missed
nsbnd(c(a = 0.35, b = 0.55), a = 0.10, b = 0.20)
#> $bnd
#>   rate pobs        lo     hi     width
#> 1    a 0.35 0.2777778 0.4375 0.1597222
#> 2    b 0.55 0.5000000 0.6875 0.1875000
#> 
#> $pairs
#>   i j   ord  bd
#> 1 a b lower 0.3
# with up to 40% missed, the order is no longer identified
nsbnd(c(a = 0.35, b = 0.55), a = 0.05, b = 0.40)$pairs
#>   i j       ord    bd
#> 1 a b undecided 0.335
# known noise rates give the true rate itself
nsbnd(0.4 * 0.8 + 0.6 * 0.1, a = c(0.1, 0.1), b = c(0.2, 0.2))$bnd
#>   rate pobs  lo  hi width
#> 1    1 0.38 0.4 0.4     0
```

## The false positive rate the other rates imply

[`impfpr()`](https://rootcoder007.github.io/crimlab/reference/impfpr.md):
Chouldechova’s (2017) identity: in any confusion table, the false
positive rate is `p / (1 - p) (1 - ppv) / ppv` `(1 - fnr)`, from the
base rate `p`, the positive predictive value `ppv` and the false
negative rate `fnr`. So two groups whose scores are equally calibrated
(the same `ppv`) and miss equally (the same `fnr`) but whose base rates
differ must have different false positive rates: no score, however it is
built, can equalise all three. The identity is arithmetic, true of any
tool on any data; it says which equalities cannot hold together, not
which one should give way. See `proof("impfpr")`.

``` r

# the same ppv and fnr, different base rates: different fpr
impfpr(p = c(a = 0.3, b = 0.5), ppv = 0.6, fnr = 0.25)
#>         a         b 
#> 0.2142857 0.5000000
```

## Which pairs a noisy ranking orders

[`rnkres()`](https://rootcoder007.github.io/crimlab/reference/rnkres.md):
People (or areas, or courts) ranked by estimated scores are ordered
correctly whatever the noise only when their estimates differ by more
than their two interval half-widths together: below that, admissible
errors can reverse them, above it none can.
[`rnkres()`](https://rootcoder007.github.io/crimlab/reference/rnkres.md)
says which pairs that holds for, and the share it fails for; the
resolution of the ranking is twice the largest half-width. A league
table whose intervals are wide orders few of its pairs (Goldstein and
Spiegelhalter 1996). See `proof("rnkres")`.

``` r

rnkres(c(a = 0.2, b = 0.35, c = 0.8), hw = c(0.1, 0.1, 0.05))
#> $id
#>       a     b    c
#> a    NA FALSE TRUE
#> b FALSE    NA TRUE
#> c  TRUE  TRUE   NA
#> 
#> $npair
#> [1] 3
#> 
#> $nunid
#> [1] 1
#> 
#> $share
#> [1] 0.3333333
#> 
#> $res
#> [1] 0.2
```

## Oaxaca-Blinder decomposition of a gap between two groups

[`oaxbd()`](https://rootcoder007.github.io/crimlab/reference/oaxbd.md):
Splits the gap in mean outcome between two groups (a sentence length
between two groups of defendants, say) into the part their covariates’
means account for and the part they do not, by least squares in each
group (Oaxaca 1973; Blinder 1973), with delta-method standard errors
(Jann 2008).

``` r

set.seed(3)
g <- rep(c("a", "b"), each = 200)
sev <- rnorm(400, ifelse(g == "a", 5, 4))
pri <- rpois(400, ifelse(g == "a", 2, 1))
len <- 2 + 0.8 * sev + 0.5 * pri + (g == "a") + rnorm(400)
r <- oaxbd(len, data.frame(sev, pri), g, ref = "b")
r$dec
#>     part        est        se         lo        hi
#> 1    gap 2.23260343 0.1384465  1.9612533 2.5039535
#> 2    exp 1.40297787 0.1668697  1.0759193 1.7300365
#> 3  unexp 0.82962555 0.1680470  0.5002595 1.1589916
#> 4   expa 1.43218736 0.1341264  1.1693045 1.6950703
#> 5 unexpa 0.80041606 0.1344082  0.5369809 1.0638512
#> 6  endow 1.40297787 0.1668697  1.0759193 1.7300365
#> 7   coef 0.80041606 0.1344082  0.5369809 1.0638512
#> 8  inter 0.02920949 0.1641780 -0.2925734 0.3509924
```

## Simulated records with a known disparity

[`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md):
Incidents in segregated areas, whose true occurrence is the same for
every group but whose recording, and the risk an algorithm predicts for
them, differ by group by a set amount: data on which an audit’s answer
is known (after Barman and Barman 2026). The first group is the
reference.

``` r

d <- fairsim(2000, seed = 1)
head(d)
#>   area group truth det     risk
#> 1   10     B     0   0 379.1753
#> 2    2     B     0   1 324.0614
#> 3   10     B     0   1 280.4492
#> 4    1     A     0   0 312.1302
#> 5    6     B     1   0 326.8925
#> 6    8     B     0   1 304.4651
# the injected ratio is 1 - 0.5
fairdf(d$det, d$group, base = "A")$tab[, c("group", "di", "dilo", "dihi")]
#>   group        di      dilo      dihi
#> 1     A 1.0000000       NaN       NaN
#> 2     B 0.5751675 0.4849938 0.6821069
```

## Noisy-OR detection by patrols

[`nordet()`](https://rootcoder007.github.io/crimlab/reference/nordet.md):
The chance that an incident is seen by patrols when each officer within
`radius` of it sees it independently with probability `p`: with `k`
officers in range, `1 - (1 - p)^k` (the detection model of Barman and
Barman 2026). With `draw = TRUE`, one detection outcome is drawn for
each incident. A model of who is seen, under stated assumptions; how
patrols are actually placed and how much they see are matters of fact
the model does not settle.

``` r

set.seed(1)
inc <- matrix(runif(20), ncol = 2)
off <- matrix(runif(10), ncol = 2)
nordet(inc, off, radius = 0.3, p = 0.7)
#>    k prob
#> 1  2 0.91
#> 2  2 0.91
#> 3  0 0.00
#> 4  2 0.91
#> 5  1 0.70
#> 6  2 0.91
#> 7  0 0.00
#> 8  0 0.00
#> 9  1 0.70
#> 10 1 0.70
```

## Fairness metrics of a binary risk tool by group

[`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md):
The confusion matrix of a risk score cut at `cut` within each group, and
the usual rates computed from it: accuracy, sensitivity (true positive
rate), specificity, false positive and false negative rates, positive
and negative predictive values, the share and number predicted positive
and the Matthews correlation, each with its ratio to the base group and,
for the first eight, a Wilson score interval (Wilson 1927). A case is
predicted positive when its score exceeds `cut`. No single metric
defines fairness, and when base rates differ between groups the error
rates and the predictive values cannot all be equal at once
(Chouldechova 2017); the table shows which differ, not which should not.

``` r

set.seed(1)
g <- rep(c("a", "b"), each = 100)
y <- rbinom(200, 1, ifelse(g == "a", 0.3, 0.5))
s <- plogis(rnorm(200, ifelse(y == 1, 0.8, -0.8)))
fair(s, y, g)$met
#>   group acc       tpr       tnr       fpr       fnr       ppv       npv ppos npos
#> 1     a 0.8 0.8125000 0.7941176 0.2058824 0.1875000 0.6500000 0.9000000 0.40   40
#> 2     b 0.8 0.8148148 0.7826087 0.2173913 0.1851852 0.8148148 0.7826087 0.54   54
#>         mcc
#> 1 0.5776155
#> 2 0.5974235
```

## Sampling Shapley values of single predictions

[`xpshap()`](https://rootcoder007.github.io/crimlab/reference/xpshap.md):
Shares out one person’s prediction among the inputs: the Shapley value
of input `j` is its contribution averaged over every order in which the
inputs could be revealed, the inputs not yet revealed taking values from
the rows of `x` (the reference data). Strumbelj and Kononenko (2014)
estimate it by sampling: each sample draws an order of the inputs and a
row `w` of `x`, and takes the prediction at the person’s values on `j`
and the inputs before it (`w`’s on the rest) less the prediction with
`j` too taken from `w`. The value is the mean over `nsamp` samples, its
standard error their standard deviation over `sqrt(nsamp)`. Over all the
inputs the values add up, in expectation, to the person’s prediction
less the mean prediction over `x`.

``` r

set.seed(5)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpshap(x, f, newx = x[1:2, ], nsamp = 100)
#> xpshap: sampling Shapley values, 100 samples each; mean prediction 0.2695
#>        pred       prior         age
#>  0.09888651 -0.10661348 -0.06044138
#>  0.27926231  0.07142368 -0.05804795
```

## Accumulated local effects of one input

[`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md):
The effect of one input on a prediction, measured only where the data
are (Apley and Zhu 2020): the range of the column is cut at its
quantiles into intervals, each row’s prediction is taken at the upper
and at the lower edge of its own interval with its other inputs as they
are, the differences are averaged within each interval and accumulated
from the lowest edge. Unlike partial dependence
([`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md))
no row is moved far from where it lies, so an input correlated with
others (prior arrests and age, say) is not credited with their effects.

``` r

set.seed(4)
age <- runif(200, 18, 60)
prior <- rpois(200, (60 - age) / 10)
x <- cbind(prior, age)
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpale(x, f, "age", nbin = 5)
#> xpale: accumulated local effects of age over 5 interval(s)
#>         z         ale
#>  18.12276  0.11761138
#>  26.93682  0.06485629
#>  35.43904  0.01050556
#>  43.48715 -0.03293827
#>  52.64466 -0.06308333
#>  59.85774 -0.07813159
```

## Ceteris paribus profiles of single rows

[`xpcp()`](https://rootcoder007.github.io/crimlab/reference/xpcp.md):
The prediction for one person as one input moves over a grid while the
person’s other inputs stay as they are: what the risk tool would have
said had only that input been different. One profile per row of `newx`;
over all the rows of `x` the profiles are the individual conditional
expectation curves of Goldstein et al. (2015), and their mean is the
partial dependence of
[`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md).
A steep profile on a protected or proxy input is a decision that turns
on that input alone.

``` r

set.seed(3)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpcp(x, f, "prior", newx = x[1:2, ], grid = c(0, 2, 4, 6))
#> xpcp: ceteris paribus profiles of 2 row(s) over 4 values of prior
#>       pred         lo        hi
#>  0.1186177 0.07546754 0.6211458
#>  0.2891560 0.08321185 0.6457738
```

## Partial dependence of a prediction function on one input

[`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md):
The average prediction when one input is set to a value for every row
and the others keep their values (Friedman 2001, sec. 8.2): at each
value `g` of a grid, the mean of the predictions with column `v` set to
`g` in every row. It is the effect of the input averaged over the data,
as a risk tool’s user would read it; where the input is correlated with
the others it averages over combinations the data never hold, which
[`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md)
avoids.

``` r

set.seed(2)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xppdp(x, f, "age", grid = 5)
#> xppdp: partial dependence on age, averaged over 100 rows
#>      grid        pd
#>  18.57225 0.3701719
#>  28.72549 0.3092730
#>  38.87872 0.2541271
#>  49.03196 0.2055543
#>  59.18520 0.1638529
```

## Permutation importance of a prediction function

[`xpimp()`](https://rootcoder007.github.io/crimlab/reference/xpimp.md):
How much worse a risk tool predicts when one of its inputs is made
useless: the column is shuffled across the rows, which keeps its
distribution but breaks its link with the outcome and with the other
columns, and the loss of the predictions against `y` is measured again
(Breiman 2001; Fisher, Rudin and Dominici 2019, who call the ratio the
model reliance). Each column is shuffled `nrep` times and the mean and
the standard deviation of the increases are reported.

``` r

set.seed(1)
x <- cbind(prior = rpois(200, 2), age = runif(200, 18, 60),
           noise = rnorm(200))
y <- 0.3 * x[, "prior"] - 0.02 * x[, "age"] + rnorm(200, sd = 0.5)
f <- function(m) 0.3 * m[, "prior"] - 0.02 * m[, "age"]
xpimp(x, f, y, nrep = 5)
#> xpimp: permutation importance, 5 shuffles per column; loss se (increase), 0.2852 before shuffling
#>   term       imp         sd
#>  prior 0.3260732 0.02932846
#>    age 0.1047989 0.01306672
#>  noise 0.0000000 0.00000000
xpimp(x, "glm", y, nrep = 5, lctl = list(nfold = 2))
#> xpimp: permutation importance, 5 shuffles per column; loss se (increase), 0.284 before shuffling
#>   term          imp          sd
#>  prior 3.469121e-01 0.031564597
#>    age 1.071655e-01 0.013107077
#>  noise 6.448013e-06 0.000732987
```
