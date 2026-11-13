<!-- written by dev/bench.R from bench.csv; do not edit -->

Machine Kobe24: Apple M2 arm64, R version 4.6.1 (2026-06-24); crimlab built with R's default flags (-g -O2).

3 runs (start time, 1-minute load average at the start and at the end): 2026-11-12 21:38 EST, load 3.72 to 3.58; 2026-11-12 22:18 EST, load 3.54 to 3.24; 2026-11-12 22:58 EST, load 3.46 to 4.19. A row whose load rose by more than 1.5 while it was timed was timed again; the loads of each row are in bench.csv (l0, l1). Each row is from the run whose ratio is the median of the 3; lo and hi are the smallest and largest ratio over the runs. The machine also runs the interactive sessions that drive this work, so its load is well above 1 even when no benchmark runs, and R ran at scheduling priority nice 10, inherited from those sessions. Pairs with forensim were not run here: loading it starts Tk, which opens a window on the desktop.

This is one of three speed tables: this one covers the functions of the disparity, forensic, policing and corrections families (with the network functions in net.csv, from dev/bench-net.R); the spatial functions are in sp.csv (dev/bench-sp.R) and the learners in dml.csv.

Median seconds per call; ratio = other / crimlab (above 1: crimlab faster). Kind: third, another package; own, rmorie (the maintainer's own package); base, a base-R implementation of the formula where no package computes the quantity.

| function | counterpart | kind | n | of | crimlab | other | ratio | lo | hi | agree |
|---|---|---|--:|---|--:|--:|--:|--:|--:|---|
| vodtst | rmorie::morie_veil_of_darkness 1.5.3 | own |         5 | groups | 1.27e-06 | 0.000167 |  132 |  124 |  137 | yes |
| vodtst | rmorie::morie_veil_of_darkness 1.5.3 | own |     1,000 | groups | 1.1e-05 | 0.00154 |  139 |  135 |  146 | yes |
| vodtst | rmorie::morie_veil_of_darkness 1.5.3 | own |   100,000 | groups | 0.00134 | 0.277 |  206 |  177 |  210 | yes |
| orrisk | rmorie::morie_relative_risk_from_or 1.5.3 | own |         1 | calls | 5.33e-07 | 1.99e-05 | 37.3 | 37.2 | 37.9 | yes |
| outtst | rmorie::morie_outcome_test 1.5.3 | own |     1,000 | people | 6.64e-06 | 0.000119 | 17.9 | 17.8 | 18.6 | yes |
| outtst | rmorie::morie_outcome_test 1.5.3 | own |   100,000 | people | 0.000557 | 0.00271 | 4.87 | 4.83 | 4.98 | yes |
| outtst | rmorie::morie_outcome_test 1.5.3 | own | 1,000,000 | people | 0.00589 | 0.0336 |  5.7 | 5.52 | 7.16 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.5.3 | own |         5 | groups | 1.31e-06 | 0.00013 | 99.2 | 94.2 |  103 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.5.3 | own |     1,000 | groups | 1.15e-05 | 0.000198 | 17.2 | 16.7 | 17.4 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.5.3 | own |   100,000 | groups | 0.00165 | 0.0106 | 6.41 |  6.3 | 7.45 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.5.3 | own |         5 | groups | 1.35e-06 | 0.000161 |  119 |  118 |  123 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.5.3 | own |        50 | groups | 3.56e-05 | 0.00147 | 41.3 | 40.5 | 43.6 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.5.3 | own |       300 | groups | 0.00125 | 0.061 | 48.8 | 42.8 | 51.4 | yes |
| gtmp | Familias::FamiliasPosterior 2.6.4 | third |        20 | loci | 1.76e-06 | 0.00554 | 3.14e+03 | 2.96e+03 | 3.81e+03 | yes |
| gtmp | Familias::FamiliasPosterior 2.6.4 | third |       200 | loci | 3.4e-06 | 0.0445 | 1.31e+04 | 1.26e+04 | 1.88e+04 | yes |
| kinidx | Familias::FamiliasPosterior 2.6.4 | third |        20 | loci | 3.49e-06 | 0.00668 | 1.92e+03 | 1.45e+03 | 2.02e+03 | yes |
| kinidx | Familias::FamiliasPosterior 2.6.4 | third |       200 | loci | 1.72e-05 | 0.116 | 6.77e+03 | 1.95e+03 | 1.02e+04 | yes |
| kinidx | forrel::kinshipLR 1.10.0 | third |        20 | loci | 3.48e-06 | 0.0035 | 1.01e+03 |  961 | 1.02e+03 | yes |
| kinidx | forrel::kinshipLR 1.10.0 | third |       200 | loci | 1.69e-05 | 0.0276 | 1.63e+03 | 1.58e+03 | 1.7e+03 | yes |
| patidx | Familias::FamiliasPosterior 2.6.4 | third |        20 | loci | 4.18e-06 | 0.0049 | 1.17e+03 | 1.06e+03 | 1.37e+03 | yes |
| patidx | Familias::FamiliasPosterior 2.6.4 | third |       200 | loci | 2.32e-05 | 0.0415 | 1.79e+03 | 1.71e+03 | 3.96e+03 | yes |
| patidx | forrel::kinshipLR 1.10.0 | third |        20 | loci | 4.43e-06 | 0.00398 |  898 |  866 |  912 | yes |
| patidx | forrel::kinshipLR 1.10.0 | third |       200 | loci | 2.29e-05 | 0.0301 | 1.32e+03 | 1.3e+03 | 1.33e+03 | yes |
| kinidx | paramlink::LR 1.1.6 | third |        20 | loci | 3.4e-06 | 0.00476 | 1.4e+03 | 1.33e+03 | 1.41e+03 | yes |
| kinidx | paramlink::LR 1.1.6 | third |       200 | loci | 1.63e-05 | 0.0426 | 2.62e+03 | 2.57e+03 | 2.73e+03 | yes |
| patidx | paramlink::LR 1.1.6 | third |        20 | loci | 4.18e-06 | 0.00539 | 1.29e+03 | 1.22e+03 | 1.31e+03 | yes |
| patidx | paramlink::LR 1.1.6 | third |       200 | loci | 2.31e-05 | 0.051 | 2.21e+03 | 2.08e+03 | 2.66e+03 | yes |
| kinidx | relSim::lrSib 1.0.0 | third |        20 | loci | 3.4e-06 | 3.81e-06 | 1.12 | 1.07 | 1.18 | yes |
| kinidx | relSim::lrSib 1.0.0 | third |       200 | loci | 1.82e-05 | 2.68e-05 | 1.47 | 1.45 | 1.47 | yes |
| kinidx (parent) | relSim::lrPC 1.0.0 | third |        20 | loci | 3.61e-06 | 3.32e-06 | 0.92 | 0.901 | 0.951 | yes |
| kinidx (parent) | relSim::lrPC 1.0.0 | third |       200 | loci | 1.71e-05 | 1.85e-05 | 1.09 | 1.09 |  1.1 | yes |
| patidx (no mother) | relSim::lrPC 1.0.0 | third |        20 | loci | 3.73e-06 | 3.28e-06 | 0.879 | 0.842 | 0.907 | yes |
| patidx (no mother) | relSim::lrPC 1.0.0 | third |       200 | loci | 1.69e-05 | 1.83e-05 | 1.09 | 1.05 | 1.09 | yes |
| mixlr | forensim::LR | third |         4 | alleles |  |  |  |  |  | not run: not loaded on this machine: forensim starts Tk, which opens a window on its desktop |
| mixlr | forensim::LR | third |         6 | alleles |  |  |  |  |  | not run: not loaded on this machine: forensim starts Tk, which opens a window on its desktop |
| mixlr | forensim::LR | third |         8 | alleles |  |  |  |  |  | not run: not loaded on this machine: forensim starts Tk, which opens a window on its desktop |
| mixlr (two people) | relSim::lrMix 1.0.0 | third |         4 | alleles | 1.35e-06 | 4.02e-06 | 2.97 |  2.8 | 3.06 | yes |
| cllr | comparison::calc.ece 1.0.8 | third |     1,000 | ratios | 5.58e-06 | 0.0131 | 2.35e+03 | 2.3e+03 | 2.48e+03 | yes |
| cllr | comparison::calc.ece 1.0.8 | third |    10,000 | ratios | 7.45e-05 | 0.75 | 1.01e+04 | 9.57e+03 | 1.03e+04 | yes |
| cllr | base R | base |     1,000 | ratios | 5.62e-06 | 1.12e-05 |    2 | 1.97 | 2.05 | yes |
| cllr | base R | base |   100,000 | ratios | 0.000886 | 0.00088 | 0.994 | 0.972 | 1.01 | yes |
| cllr | base R | base | 1,000,000 | ratios | 0.00891 | 0.0109 | 1.23 |  1.2 | 1.26 | yes |
| cllrmin | comparison::calc.ece 1.0.8 | third |     1,000 | ratios | 0.000137 | 0.0113 | 82.8 | 82.6 | 84.6 | yes |
| cllrmin | comparison::calc.ece 1.0.8 | third |    10,000 | ratios | 0.000453 | 0.783 | 1.73e+03 | 1.72e+03 | 1.75e+03 | yes |
| cllr | idiolect:::Cllr 1.2.0 | third |     1,000 | ratios | 5.53e-06 | 2.69e-05 | 4.86 | 4.77 | 4.94 | yes |
| cllr | idiolect:::Cllr 1.2.0 | third |   100,000 | ratios | 0.000884 | 0.00183 | 2.07 | 2.06 | 2.09 | yes |
| cllr | idiolect:::Cllr 1.2.0 | third | 1,000,000 | ratios | 0.00887 | 0.0207 | 2.33 | 2.29 | 2.37 | yes |
| cllrmin | idiolect:::Cllr 1.2.0 | third |     1,000 | ratios | 0.000137 | 0.000157 | 1.15 | 1.15 | 1.18 | yes |
| cllrmin | idiolect:::Cllr 1.2.0 | third |   100,000 | ratios | 0.00388 | 0.00903 | 2.33 | 2.26 | 2.33 | yes |
| cllrmin | idiolect:::Cllr 1.2.0 | third | 1,000,000 | ratios | 0.0454 | 0.158 | 3.47 |  2.3 |  4.5 | yes |
| lrcal | idiolect:::opt.llr 1.2.0 | third |     1,000 | ratios | 0.000138 | 0.000134 | 0.968 | 0.964 | 0.97 | yes |
| lrcal | idiolect:::opt.llr 1.2.0 | third |   100,000 | ratios | 0.00288 | 0.00655 | 2.27 | 2.24 |  2.3 | yes |
| lrcal | idiolect:::opt.llr 1.2.0 | third | 1,000,000 | ratios | 0.0364 | 0.172 | 4.72 | 3.69 | 5.25 | yes |
| lrcal | comparison::calibrate.set 1.0.8 | third |     1,000 | ratios | 0.000136 | 0.0118 | 87.2 | 84.5 | 89.9 | yes |
| lrcal | comparison::calibrate.set 1.0.8 | third |    10,000 | ratios | 0.000368 | 0.744 | 2.02e+03 | 1.97e+03 | 2.17e+03 | yes |
| eceplt | comparison::calc.ece 1.0.8 | third |     1,000 | ratios | 0.00238 | 0.0151 | 6.33 | 5.51 | 6.34 | yes |
| eceplt | comparison::calc.ece 1.0.8 | third |    10,000 | ratios | 0.0174 | 0.737 | 42.3 | 40.9 | 50.5 | yes |
| lrpost | base R | base |     1,000 | ratios | 5.62e-06 | 3.08e-06 | 0.547 | 0.471 | 0.752 | yes |
| lrpost | base R | base |   100,000 | ratios | 0.000321 | 0.000123 | 0.384 | 0.382 | 0.409 | yes |
| lrpost | base R | base | 1,000,000 | ratios | 0.00348 | 0.00155 | 0.446 | 0.422 | 0.458 | yes |
| dbsrch | base R | base |         1 | calls | 5.33e-07 | 6.15e-07 | 1.15 |    1 | 1.15 | yes |
| tipplt | base R | base |     1,000 | ratios | 0.000686 | 0.000638 | 0.93 | 0.926 | 0.934 | yes |
| tipplt | base R | base |   100,000 | ratios | 0.00523 | 0.00769 | 1.47 | 1.47 | 1.47 | yes |
| tipplt | base R | base | 1,000,000 | ratios | 0.0562 | 0.0782 | 1.39 | 0.974 |  1.4 | yes |
| clrsmp | base R | base |         5 | offence types | 1.68e-06 | 2.79e-06 | 1.66 | 1.64 | 1.69 | yes |
| clrsmp | base R | base |     1,000 | offence types | 1.18e-05 | 2.25e-05 | 1.91 | 1.87 | 1.94 | yes |
| clrsmp | base R | base |   100,000 | offence types | 0.000988 | 0.00193 | 1.96 | 1.95 | 1.97 | yes |
| clrsmp | DasGuptR::dgnpop 2.2.0 | third |         5 | offence types | 1.72e-06 | 0.00253 | 1.47e+03 | 1.39e+03 | 1.54e+03 | yes |
| clrsmp | DasGuptR::dgnpop 2.2.0 | third |     1,000 | offence types | 1.2e-05 | 0.00813 |  679 |  678 |  717 | yes |
| clrsmp | DasGuptR::dgnpop 2.2.0 | third |   100,000 | offence types | 0.00101 | 0.697 |  689 |  683 |  690 | yes |
| clrsmp | DemoDecomp::horiuchi 1.14.1 | third |         5 | offence types | 1.72e-06 | 0.000465 |  270 |  258 |  282 | yes |
| clrsmp | DemoDecomp::horiuchi 1.14.1 | third |     1,000 | offence types | 1.18e-05 | 2.12 | 1.8e+05 | 1.77e+05 | 1.83e+05 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.5.3 | own |         5 | groups | 1.31e-06 | 0.000215 |  163 |  154 |  166 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.5.3 | own |     1,000 | groups | 1.12e-05 | 0.000296 | 26.5 | 26.5 | 26.7 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.5.3 | own |   100,000 | groups | 0.00128 | 0.00992 | 7.74 | 7.43 | 8.67 | yes |
| misrd | rmorie::morie_misclassified_rd 1.5.3 | own |         1 | calls | 6.56e-07 | 3.44e-06 | 5.25 | 5.06 | 5.31 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.5.3 | own |     1,000 | places | 7.75e-06 | 6.23e-05 | 8.04 | 7.83 | 8.08 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.5.3 | own |   100,000 | places | 0.000599 | 0.00611 | 10.2 | 10.1 | 10.2 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.5.3 | own | 1,000,000 | places | 0.00615 | 0.0813 | 13.2 | 12.9 | 13.7 | yes |
| incap (share) | rmorie::morie_incapacitation 1.5.3 | own |         5 | groups | 1.43e-06 | 0.000123 | 85.5 | 82.1 | 86.1 | yes |
| incap (share) | rmorie::morie_incapacitation 1.5.3 | own |     1,000 | groups | 2.89e-05 | 0.00015 | 5.18 | 5.13 |  5.5 | yes |
| incap (share) | rmorie::morie_incapacitation 1.5.3 | own |   100,000 | groups | 0.00276 | 0.00256 | 0.928 | 0.925 | 0.931 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.5.3 | own |         5 | groups | 1.27e-06 | 0.00012 | 94.2 | 93.9 | 98.7 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.5.3 | own |     1,000 | groups | 1.09e-05 | 0.000141 | 12.9 | 12.6 | 12.9 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.5.3 | own |   100,000 | groups | 0.000876 | 0.00186 | 2.13 | 2.12 | 2.16 | yes |
| incapr | rmorie::morie_incapacitation_career 1.5.3 | own |        10 | career periods | 6.56e-07 | 4.51e-06 | 6.87 | 6.53 | 7.33 | yes |
| incapr | rmorie::morie_incapacitation_career 1.5.3 | own |     1,000 | career periods | 2.87e-06 | 2.27e-05 |  7.9 | 7.73 | 8.03 | yes |
| incapr | rmorie::morie_incapacitation_career 1.5.3 | own |   100,000 | career periods | 0.000235 | 0.00172 | 7.29 | 7.29 | 7.31 | yes |
| sellab | rmorie::morie_selective_labels 1.5.3 | own |     1,000 | defendants | 6.48e-06 | 4.78e-05 | 7.37 | 7.34 | 7.67 | yes |
| sellab | rmorie::morie_selective_labels 1.5.3 | own |   100,000 | defendants | 0.00097 | 0.00537 | 5.53 | 5.53 | 5.61 | yes |
| sellab | rmorie::morie_selective_labels 1.5.3 | own | 1,000,000 | defendants | 0.01 | 0.0674 | 6.72 | 6.48 | 6.81 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.5.3 | own |     1,000 | stays | 2.09e-06 | 6.97e-06 | 3.33 | 3.29 | 3.41 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.5.3 | own |   100,000 | stays | 0.000171 | 0.000508 | 2.96 | 2.94 | 2.99 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.5.3 | own | 1,000,000 | stays | 0.00171 | 0.00501 | 2.93 | 2.93 | 2.99 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.5.3 | own |     1,000 | people | 2.09e-05 | 8.32e-05 | 3.98 |  3.7 | 4.08 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.5.3 | own |   100,000 | people | 0.00198 | 0.00708 | 3.58 | 3.54 | 3.67 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.5.3 | own | 1,000,000 | people | 0.0208 | 0.0828 | 3.98 | 3.77 | 4.23 | yes |
| hzdsel | rmorie::morie_hazard_selection 1.5.3 | own |         1 | calls | 1.15e-06 | 2.75e-06 | 2.39 | 2.34 | 2.45 | yes |
| proof | base R | base |         1 | calls | 3.94e-06 | 1.93e-05 |  4.9 |  4.6 | 4.98 | yes |
