#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
"""Writes the guides in vignettes/ from the help pages, one guide per group
of families of _pkgdown.yml, so that every exported function is shown once,
with the example of its help page and, where the package draws its result,
the plot.

    Rscript dev/vigrd.R . rd.json                 the help pages
    Rscript dev/vigplot.R crimlab rd.json plot.tsv   what each example can plot
    python3 dev/vig.py . rd.json plot.tsv

The guides are generated: change a help page (or GROUPS below) and run this
again. vignettes/crimlab.Rmd, the introduction, is written by hand and left
alone.
"""
import json
import pathlib
import re
import sys

import yaml

repo = pathlib.Path(sys.argv[1])
pages = json.load(open(sys.argv[2]))
plot = {}
call = {}
slow = {}
for line in open(sys.argv[3]).read().split("\n"):
    if line:
        c = line.split("\t")
        f, ok, var, sec = c[:4]
        if ok == "ok" and var:
            plot[f] = var
        elif ok == "ok" and len(c) > 5 and c[5]:
            call[f] = c[5].replace("\\n", "\n")
        slow[f] = float(sec)

# guide file, title, one sentence, the families of _pkgdown.yml it holds
GROUPS = [
    ("disparities-and-checks", "Disparities, benchmarks and checks of a reported number",
     "Tests of disparity in stops and outcomes, and checks that a published ratio or table says what it is taken to say.",
     ["Disparities", "Checks", "Shares"]),
    ("forensics", "Forensic evidence: likelihood ratios and their calibration",
     "Likelihood ratios for forensic comparisons, how well a system that reports them is calibrated, and the plots that show it.",
     ["Forensics"]),
    ("eyewitness-and-toxicology", "Eyewitness identification, toxicology and taphonomy",
     "Lineup performance, calibration of a toxicology method, post-mortem redistribution and the time since death, and alcohol and drug use: a blood alcohol concentration from what was drunk and the scores of screening questionnaires.",
     ["Eyewitness identification", "Toxicology and taphonomy", "Alcohol and drug use"]),
    ("policing-and-deterrence", "Policing, police-reported crime and deterrence",
     "Clearance rates, exposure bounds, crime severity and reporting, and how inspection or arrest risk deters.",
     ["Policing", "Police-reported crime", "Deterrence"]),
    ("corrections-and-custody", "Corrections, custody and court flow",
     "Incapacitation, bounds on sentence effects, stock and flow in custody and court backlogs, and segregation against the Mandela Rules.",
     ["Corrections", "Custody and court flow", "Segregation and the Mandela Rules"]),
    ("placement-and-alerts", "Placement volatility and alert complexity",
     "How often a person in custody is moved within a year, how many alerts a record carries, and what follows from both.",
     ["Placement volatility and alert complexity"]),
    ("courts-and-life-course", "Courts, sentencing and the life course",
     "Judge-leniency designs, sentencing models, the age-crime curve, cumulative risk and the justice pipeline.",
     ["Courts and sentencing", "Life course and the justice pipeline"]),
    ("networks", "Networks of co-offending",
     "Building a co-offending network and measuring who is central in it.",
     ["Networks"]),
    ("spatial-geometry", "Spatial geometry and centrography",
     "Distances, projections, centres, ellipses and the other geometry that spatial analysis of crime starts from.",
     ["Spatial"]),
    ("spatial-crime-analysis", "Spatial crime analysis: hot spots, clusters and near repeats",
     "Hot spots, clustering, scan statistics, near-repeat tests and surfaces of crime density.",
     ["Spatial crime analysis"]),
    ("geostatistics-and-epidemiology", "Geostatistics and spatial epidemiology",
     "Variograms and kriging, and area-level models of rates with spatial dependence.",
     ["Geostatistics", "Spatial epidemiology"]),
    ("risk-tools", "Risk tools and predictive policing",
     "Fairness of a risk score across groups, its calibration, and audits of a predictive-policing rule.",
     ["Risk tools and predictive policing"]),
    ("counts-rates-and-tables", "Counts, rates, concentration and published tables",
     "Count regression, standardised rates, cross-tables, concentration of crime, the dark figure, and tables published in bands.",
     ["Counts, rates and tables", "Concentration", "Dark figure", "Published tables"]),
    ("change-tables-and-drift", "Change tables, drift between releases and region maps",
     "Year-over-year tables in several formats, tests that a new release of data still looks like the last, and checks of a region map.",
     ["Change tables out", "Drift between releases", "Region maps"]),
    ("describing-a-table", "Describing a table and keeping its categories whole",
     "Profiles of a data frame, frequency and correlation tables, and guards against categories that change silently.",
     ["Describing a table", "Categorical integrity", "Text"]),
    ("learning-and-presets", "Learners, presets and data larger than memory",
     "Penalised regression and forests written in the package, checks before a model is fitted, and summaries of data read in one pass.",
     ["Learning", "Presets", "Data larger than memory"]),
    ("causal-estimates", "Causal estimates and regression discontinuity",
     "Double machine learning, balance, sensitivity of an estimate to what was not measured, and regression discontinuity.",
     ["Causal estimates", "Regression discontinuity"]),
    ("experiments-and-policy", "Experiments and policy evaluation",
     "Randomisation, power, falsification tests, difference-in-differences, event studies and synthetic controls.",
     ["Experiments", "Policy evaluation"]),
    ("time-and-contagion", "Events in time: contagion, trends and operations",
     "Self-exciting processes, trends, change points, survival, forecasts and the other methods for events in time.",
     ["Time, contagion and operations", "Counts over time and across categories"]),
    ("surveys-and-summaries", "Surveys and summary statistics",
     "Survey estimates with their design, and the summary statistics the rest of the package is built on.",
     ["Surveys", "Summary statistics"]),
    ("indigenous-health", "Indigenous health and wellbeing",
     "Tables of First Nations, Inuit and M\\u00e9tis health and wellbeing: reading them, their quality flags, and gaps with their uncertainty.",
     ["Indigenous health and wellbeing"]),
    ("air-pollution", "Air pollution and the footprint of a computation",
     "Dispersion of a pollutant, the health burden attributed to it, and the energy and carbon cost of running an analysis.",
     ["Air pollution", "The footprint of a computation"]),
    ("exposure-and-crime", "Environmental exposure and crime",
     "Designs that estimate the effect of an environmental exposure on crime: fixed effects, instruments and their diagnostics.",
     ["Exposure and crime"]),
    ("environmental-crime-and-rescue", "Environmental crime, panels and search and rescue",
     "Panel models and tests used for environmental crime, and the search theory behind a search-and-rescue plan.",
     ["Environmental crime", "Search and rescue"]),
    ("scaling", "Scaling and ideal points",
     "Agreement scores and the scaling methods that place voters, judges or items on a line or in a plane.",
     ["Scaling and ideal points"]),
    ("oversight-reports", "Police oversight reports and open data",
     "Reading the text of police oversight reports into fields, and loading tables of the project's open data service.",
     ["Police oversight reports", "Open data"]),
    ("plots-and-proofs", "Plot helpers, proofs and the rest",
     "The plot helpers shared by the package, the theorems behind a function, and the functions no other guide holds.",
     ["Plots", "Proofs", "Package"]),
]


def quote(s):
    """The quoted code of a help page as code: backticks, or a code tag where
    knitr would take the backticks for R to run (text that begins with r and a space)."""
    def one(m):
        t = m.group(1)
        return "<code>%s</code>" % t if t == "r" or t.startswith("r ") else "`%s`" % t
    return re.sub("\u2018([^\u2018\u2019]*)\u2019", one, s)


HEAD = '''---
title: "{title}"
output: rmarkdown::html_vignette
vignette: >
  %\\VignetteIndexEntry{{{title}}}
  %\\VignetteEngine{{knitr::rmarkdown}}
  %\\VignetteEncoding{{UTF-8}}
---

```{{r, include = FALSE}}
# the vignette runs in full in GitHub CI and for the package tarball (NOT_CRAN),
# on r-universe (MY_UNIVERSE), on the package website (IN_PKGDOWN) and in an
# interactive session; when CRAN re-builds it in its check, none of these is set
# and only this chunk runs
if (!interactive() &&
      !any(nzchar(Sys.getenv(c("NOT_CRAN", "MY_UNIVERSE", "IN_PKGDOWN"))))) {{
  knitr::knit_exit()
}}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5,
                      fig.height = 4, fig.align = "center", dpi = 72)
options(width = 88)
set.seed(1)
library({pkg})
```

{lead}

Each section is one help page: what the function is for, the example of its
help page as it runs, and the plot of its result where the package draws
one. `?name` gives the arguments and the value in full.

'''

pkgname = re.search(r"^Package:\s*(\S+)", (repo / "DESCRIPTION").read_text(), re.M).group(1)
ref = yaml.safe_load((repo / "_pkgdown.yml").read_text())["reference"]
fam = {s["title"]: [c for c in (s.get("contents") or []) if isinstance(c, str)]
       for s in ref if "title" in s}
want = [f for g in GROUPS for f in g[3]]
assert sorted(want) == sorted(fam), ("families of _pkgdown.yml and of GROUPS differ",
                                    sorted(set(fam) ^ set(want)))

by_alias = {}
for p in pages:
    for a in p["aliases"]:
        by_alias.setdefault(a, p)

exports = set(re.findall(r"^export\((.*)\)$", (repo / "NAMESPACE").read_text(), re.M))
used = set()
out = repo / "vignettes"
written = []


def saved(code):
    """The call of an example that only printed its result, saved as res.

    Arguments aligned under the opening bracket move right with it; a call
    whose arguments hang two spaces in stays as it is. Where that makes a
    line longer than 80 characters the call goes on the line below the
    arrow, and None is returned if it does not fit there either."""
    ln = code.split("\n")
    hang = any(re.match(r" {3,}\S", x) for x in ln[1:])
    for first, pad in (("res <- " + ln[0], 7), ("res <-\n  " + ln[0], 2)):
        move = hang or pad == 2
        out = [first] + [(" " * pad + x if move and x.strip() else x)
                         for x in ln[1:]]
        text = "\n".join(out)
        if max(len(x) for x in text.split("\n")) <= 80:
            return text
    return None


def section(p):
    names = [a for a in p["aliases"] if a in exports] or p["aliases"][:1]
    calls = ", ".join("`%s()`" % n for n in names[:4]) + (" and others" if len(names) > 4 else "")
    desc = quote(p["desc"])
    s = ["## %s" % p["title"].rstrip("."), "", "%s: %s" % (calls, desc), ""]
    lab = re.sub(r"[^A-Za-z0-9]+", "-", p["file"][:-3]).strip("-")
    ex = "\n".join(l.rstrip() for l in p["ex"].split("\n"))
    s += ["```{r %s}" % lab, ex, "```", ""]
    v = plot.get(p["file"])
    if v:
        s += ["```{r %s-plot}" % lab, "plot(%s)" % v, "```", ""]
    elif p["file"] in call and saved(call[p["file"]]):
        s += ["```{r %s-plot}" % lab, saved(call[p["file"]]), "plot(res)",
              "```", ""]
    return "\n".join(s)


for k, (name, title, lead, fams) in enumerate(GROUPS):
    body = []
    n = 0
    for f in fams:
        secs = []
        for fn in fam[f]:
            p = by_alias.get(fn)
            if p is None or p["file"] in used:
                continue
            used.add(p["file"])
            secs.append(section(p))
        if not secs:
            continue
        if len(fams) > 1:
            body.append("# %s\n" % f)
        body += secs
        n += len(secs)
    if k == len(GROUPS) - 1:
        rest = [p for p in pages if p["file"] not in used
                and any(a in exports for a in p["aliases"])
                and not p["file"].endswith("-package.Rd")]
        if rest:
            body.append("# Other functions\n")
            for p in rest:
                used.add(p["file"])
                body.append(section(p))
                n += 1
    text = HEAD.format(title=title, lead=lead.encode().decode("unicode_escape"), pkg=pkgname) + "\n".join(body)
    (out / (name + ".Rmd")).write_text(text.rstrip("\n") + "\n", encoding="utf-8")
    written.append((name, n, sum(slow.get(pg, 0) for pg in used)))

# the website lists every guide
yml = (repo / "_pkgdown.yml").read_text()
m = re.search(r"(articles:\n  - title: Guides\n    navbar: ~\n    contents:\n)((?:      - .*\n)+)", yml)
assert m, "_pkgdown.yml: the Guides list was not found"
names = ["crimlab"] + [g[0] for g in GROUPS]
(repo / "_pkgdown.yml").write_text(yml[:m.start(2)] + "".join("      - %s\n" % n for n in names) + yml[m.end(2):])

covered = {a for p in pages if p["file"] in used for a in p["aliases"]}
missing = sorted(exports - covered)
for name, n, _ in written:
    print("%-34s %3d pages" % (name, n))
print(len(written), "guides,", len(used), "help pages,", len(exports & covered), "of", len(exports), "exports shown")
if missing:
    print("NOT SHOWN:", " ".join(missing))
    sys.exit(1)
