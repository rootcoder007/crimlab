# Energy and carbon footprint of a computation

Evaluates `expr` and reports the electricity it used and the
CO2-equivalent that electricity implies, by one of two methods, for the
grid it ran on, offline. Use it to report the energy and emissions of an
analysis next to its results (as journals and funders increasingly ask),
or to compare two ways of running the same computation. The result keeps
the value of `expr`, and lists every figure that was assumed rather than
measured.

## Usage

``` r
co2fp(
  expr,
  loc = NULL,
  meth = c("auto", "rapl", "model"),
  cpupwrw = NULL,
  usage = "process",
  loadcur = c("linear", "codecarbon"),
  cores = NULL,
  memgb = NULL,
  mempwrwpergb = 0.3725,
  pue = 1,
  rapldir = "/sys/class/powercap",
  procstat = "/proc/stat",
  devices = NULL
)

# S3 method for class 'co2fp'
print(x, ...)
```

## Arguments

- expr:

  The expression to evaluate. It is evaluated once, in the calling
  environment, when the function reaches it (lazily, as a promise).

- loc:

  The location, passed to
  [`co2int`](https://rootcoder007.github.io/crimlab/reference/co2int.md).
  It is a location code, alpha-3 code or name, a number in gCO2e/kWh, or
  `NULL` (the default) to detect it.

- meth:

  Character scalar: `"auto"` (the default; RAPL when readable, else the
  model), `"rapl"` (measured; an error when the counters are not
  readable) or `"model"` (always the Green Algorithms model).

- cpupwrw:

  The CPU's total thermal design power in watts (the model's
  \\P\_{cpu}\\); `NULL` uses codecarbon's unknown-CPU fallback.

- usage:

  `"process"` (measured from this process's CPU time), `"machine"`
  (measured from the whole machine's CPU accounting), or a number in
  between 0 and 1. `NULL` means `"process"`.

- loadcur:

  Character scalar: `"linear"` (the default; power proportional to
  utilisation, as Green Algorithms) or `"codecarbon"` (codecarbon's
  curves, see Details). Applies only when `cpupwrw` is given.

- cores:

  Number of logical cores the process utilisation is measured against;
  `NULL` detects it.

- memgb:

  Memory charged to the computation in GB; `NULL` detects the machine's
  total memory where it can, else assumes 8 GB.

- mempwrwpergb:

  Memory power per GB (default 0.3725 W/GB, Green Algorithms; codecarbon
  v2 used 0.375).

- pue:

  Power usage effectiveness: 1 for a laptop or desktop (default); 1.67
  is the Green Algorithms default for a data centre, 1.56 the Uptime
  Institute 2024 global average. Must be a single number of at least 1.

- rapldir:

  Character scalar, the powercap directory holding the RAPL domains;
  default `"/sys/class/powercap"` (Linux).

- procstat:

  Character scalar, the kernel CPU accounting file read for
  `usage = "machine"`; default `"/proc/stat"` (Linux).

- devices:

  `NULL` (the default), or a data frame of the other devices' measured
  energy, one row per device: `kind`, `vendor`, `name`, `watts` and/or
  `energy` (kWh), and `how` (see Details).

- x:

  A result of `co2fp()` or
  [`fpeqv()`](https://rootcoder007.github.io/crimlab/reference/fpeqv.md).

- ...:

  Not used.

## Value

A list of class `co2fp` with elements `value` (the result of `expr`),
`durs` (numeric, wall seconds), `cputms` (numeric, CPU seconds of this
process and its children), `usage` (numeric in \\\[0, 1\]\\, the
utilisation), `usagemode` (`"process"`, `"machine"` or `"fixed"`),
`cores` (integer), `memgb` (numeric), `enekwh` (a list of numerics
`cpu`, `memory` and `total`, in kWh, PUE included, and `devices` when
`devices` is given), `co2eg` (numeric, grams of CO2-equivalent),
`co2int` (the list returned by
[`co2int`](https://rootcoder007.github.io/crimlab/reference/co2int.md),
with its year and source), `method` (`"rapl"` or `"model"`, the method
actually used), `assumptions` (character vector of every default or
fallback that stood in for a measurement) and `citations` (character
vector). With `devices`, also `devices`: a data frame with `kind`,
`vendor`, `name`, `how`, `kwh` (the device's energy before PUE, `NA`
when not counted), `counted` (logical) and `note` (how the energy was
taken, or why it was not counted). The print method shows a short
summary.

## Details

**Measured** (`meth = "rapl"`): the CPU package energy counters (RAPL,
microjoules, read from `rapldir`) before and after, with wrap-around
handled through `max_energy_range_uj`. This is what codecarbon does on
Linux; the counters cover the CPU packages only, so memory is still
modelled. Available when
[`raplok()`](https://rootcoder007.github.io/crimlab/reference/raplok.md)
is `TRUE`.

**Modelled** (`meth = "model"`): the Green Algorithms formula
(Lannelongue, Grealey and Inouye 2021), \$\$E = t \\(P\_{cpu}\\ u + M\\
P\_{mem})\\ PUE\$\$ with \\t\\ the running time, \\P\_{cpu}\\ the CPU's
thermal design power in watts, \\u\\ its utilisation, \\M\\ the memory
in GB, \\P\_{mem} = 0.3725\\ W/GB and \\PUE\\ the data-centre overhead.
When the CPU's power is unknown (`cpupwrw = NULL`) the fallback is
codecarbon's: 50\\ running time whatever the measured utilisation (the
utilisation is still reported, but not applied).

**Utilisation** is measured, not assumed, in one of two ways.
`usage = "process"` (the default) is this R process's own CPU time,
children included, divided by the wall time and the core count:
codecarbon's process tracking mode, which charges only what the
computation itself used. `usage = "machine"` is the whole machine's busy
share over the run, from `/proc/stat` on Linux: codecarbon's machine
tracking mode, right when the computation is the only thing the machine
is doing or when the question is what the machine drew. Where
`/proc/stat` is absent the process figure is used and the assumptions
say so. A number between 0 and 1 fixes the utilisation, as the Green
Algorithms calculator asks for it. `loadcur` chooses how power scales
with that utilisation: linearly (Green Algorithms) or along codecarbon's
curves, a floor at 10\\ \\0.1 + 0.9 u^3\\ for a machine, so a figure can
be made comparable with either tool.

**Other devices** (`devices`): the CPU counters and the model cover the
CPU and memory only. The energy of a graphics card or another
accelerator, measured by its own tools, is added through a table with
one row per device and the columns `kind` (such as `"gpu"`), `vendor`,
`name`, `watts` (the device's average power over the run, in W) and
`energy` (its energy over the run, in kWh) (one of the two may be left
out, and a missing value means not measured) and `how` (how the figure
was obtained). A device's energy is its `energy`, or else `watts` times
the running time; that energy is added to the CPU and memory energy
before the PUE and the carbon intensity are applied. A row of kind
`"cpu"` or `"memory"` is not counted again, nor is a row with neither
figure; the result's `devices` table says which rows were counted and
why the others were not. Without `devices` the result is as described
below, with no `devices` element.

Emissions are \\E \times\\ the carbon intensity of the location
([`co2int`](https://rootcoder007.github.io/crimlab/reference/co2int.md):
a table compiled into the package of every country at its latest
published year and of sub-national zones, so no network access is needed
anywhere). With `loc = NULL` the location is detected from the
environment, the time zone or the locale
([`dtctloc`](https://rootcoder007.github.io/crimlab/reference/dtctloc.md))
and the detection is recorded. Every default that stood in for a
measurement is listed in `assumptions`, so a report can say what was
measured and what was assumed.

The energy figures include the PUE factor (CPU, memory and total alike).
The running time has a floor of the CPU time divided by the cores, so a
computation faster than the clock's resolution does not come out at
zero. Short computations (well under a second) give figures dominated by
timer resolution; measure a representative workload.

Every argument is checked, and the carbon intensity looked up, before
`expr` is evaluated, so a bad argument never costs the computation.
Errors: an unknown `loc`; `pue` below 1 or not a single finite number;
`cpupwrw` or `cores` not a single positive number (`cores` also whole);
`memgb` or `mempwrwpergb` negative; `usage` other than `"process"`,
`"machine"` or one number in \\\[0, 1\]\\; `meth = "rapl"` where the
counters are not readable; a `devices` table without its columns or with
a negative or infinite figure. An error inside `expr` propagates
unchanged. Nothing is written and nothing is fetched.

## References

Lannelongue, L., Grealey, J. and Inouye, M. (2021). Green Algorithms:
quantifying the carbon footprint of computation. Advanced Science 8(12),
2100707.
[doi:10.1002/advs.202100707](https://doi.org/10.1002/advs.202100707)

Courty, B. et al. (2024). CodeCarbon: estimate and track carbon
emissions from machine learning computing.
[doi:10.5281/zenodo.4658424](https://doi.org/10.5281/zenodo.4658424)

Lacoste, A., Luccioni, A., Schmidt, V. and Dandres, T. (2019).
Quantifying the carbon emissions of machine learning. arXiv:1910.09700.

## See also

[`co2int`](https://rootcoder007.github.io/crimlab/reference/co2int.md),
[`raplok`](https://rootcoder007.github.io/crimlab/reference/raplok.md),
[`fpeqv`](https://rootcoder007.github.io/crimlab/reference/fpeqv.md).

## Examples

``` r
fp <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON")
fp$enekwh$total
#> [1] 2.684137e-08
fp$co2eg
#> [1] 2.441759e-06
fp$co2int$year
#> [1] 2024
fp$assumptions
#> [1] "CPU power unknown: 85 W x 50% (codecarbon fallback)"
#> [2] "memory charged: the machine's total 15.6 GB"        
print(fp)
#> Computation footprint (modelled)
#>   Time: 0.002 s wall, 0.001 s CPU; utilisation 0.12 (process) on 4 cores
#>   Energy: 2.68e-08 kWh (cpu 2.36e-08, memory 3.23e-09)
#>   CO2e: 2.44e-06 g at 90.97 gCO2e/kWh (Canada, Ontario, 2024)
#>   Assumed: CPU power unknown: 85 W x 50% (codecarbon fallback); memory charged: the machine's total 15.6 GB
# wherever you are, offline: the location is detected and recorded
co2fp(1 + 1, cpupwrw = 45, memgb = 8)$assumptions
#> [1] "location WORLD via fallback, no location found in environment, time zone or locale"
# the whole machine's utilisation, on codecarbon's machine curve
co2fp(
  1 + 1,
  loc = "FR",
  usage = "machine",
  loadcur = "codecarbon",
  cpupwrw = 45,
  memgb = 8
)$usagemode
#> [1] "machine"
# a fixed utilisation and a data-centre PUE, with a carbon intensity given
# as a number
fp2 <- co2fp(
  sum(sqrt(seq_len(2e5))),
  loc = 30,
  usage = 0.5,
  cpupwrw = 65,
  memgb = 16,
  pue = 1.67,
  meth = "model"
)
fp2$enekwh
#> $cpu
#> [1] 3.015278e-08
#> 
#> $memory
#> [1] 5.529556e-09
#> 
#> $total
#> [1] 3.568233e-08
#> 
fpeqv(fp2$co2eg)
#> Footprint in everyday terms
#>   Car travel: 6.117e-09 km
#>   Tree uptake: 1.167e-09 tree-months
# a graphics card measured at 120 W on average (a simulated reading) and
# a CPU row, which is not counted twice
dev <- data.frame(
  kind = c("gpu", "cpu"), vendor = c("example", "example"),
  name = c("card 0", "package 0"), watts = c(120, 35),
  how = c("simulated reading", "simulated reading")
)
fp3 <- co2fp(sum(sqrt(seq_len(2e5))), loc = "CA-ON", cpupwrw = 65,
             memgb = 16, meth = "model", devices = dev)
fp3$devices
#>   kind  vendor      name               how   kwh counted
#> 1  gpu example    card 0 simulated reading 1e-07    TRUE
#> 2  cpu example package 0 simulated reading   NaN   FALSE
#>                                      note
#> 1              average power over the run
#> 2 already counted as CPU or memory energy
# a bad argument is refused before the expression runs
try(co2fp(Sys.sleep(10), pue = 0.9))
#> Error : co2fp: `pue` must be a single number >= 1.
```
