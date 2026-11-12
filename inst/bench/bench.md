<!-- written by dev/bench.R from bench.csv; do not edit -->

R 4.5.0, Cortex-A76 aarch64, 2026-11-12. Median seconds per call; ratio = other / crimlab (above 1: crimlab faster).

| function | counterpart | n | of | crimlab | other | ratio | agree |
|---|---|--:|---|--:|--:|--:|---|
| vodtst | rmorie::morie_veil_of_darkness 1.4.3 |         5 | groups | 3.09e-06 | 0.00104 |  335 | yes |
| vodtst | rmorie::morie_veil_of_darkness 1.4.3 |     1,000 | groups | 2.45e-05 | 0.00503 |  205 | yes |
| vodtst | rmorie::morie_veil_of_darkness 1.4.3 |   100,000 | groups | 0.00324 | 0.951 |  294 | yes |
| orrisk | rmorie::morie_relative_risk_from_or 1.4.3 |         1 | calls | 1.93e-06 | 8.3e-05 | 43.1 | yes |
| outtst | rmorie::morie_outcome_test 1.4.3 |     1,000 | people | 1.13e-05 | 0.000704 | 62.5 | yes |
| outtst | rmorie::morie_outcome_test 1.4.3 |   100,000 | people | 0.00109 | 0.015 | 13.7 | yes |
| outtst | rmorie::morie_outcome_test 1.4.3 | 1,000,000 | people | 0.0235 | 0.144 | 6.13 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.4.3 |         5 | groups | 3.13e-06 | 0.000803 |  257 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.4.3 |     1,000 | groups | 2.58e-05 | 0.000797 | 30.9 | yes |
| dspbm | rmorie::morie_disparity_benchmark 1.4.3 |   100,000 | groups | 0.00259 | 0.0533 | 20.6 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.4.3 |         5 | groups | 5.63e-06 | 0.000718 |  128 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.4.3 |        50 | groups | 0.000119 | 0.00501 |   42 | yes |
| dyadrt | rmorie::morie_interracial_rates 1.4.3 |       300 | groups | 0.00514 | 0.205 | 39.9 | yes |
| gtmp | Familias::FamiliasPosterior 2.6.4 |        20 | loci | 5.61e-06 | 0.0216 | 3.85e+03 | yes |
| gtmp | Familias::FamiliasPosterior 2.6.4 |       200 | loci | 1.28e-05 | 0.194 | 1.52e+04 | yes |
| kinidx | Familias::FamiliasPosterior 2.6.4 |        20 | loci | 1.13e-05 | 0.016 | 1.41e+03 | yes |
| kinidx | Familias::FamiliasPosterior 2.6.4 |       200 | loci | 6.17e-05 | 0.227 | 3.68e+03 | yes |
| kinidx | forrel::kinshipLR 1.10.0 |        20 | loci | 1.13e-05 | 0.0163 | 1.44e+03 | yes |
| kinidx | forrel::kinshipLR 1.10.0 |       200 | loci | 6.11e-05 | 0.142 | 2.32e+03 | yes |
| patidx | Familias::FamiliasPosterior 2.6.4 |        20 | loci | 1.35e-05 | 0.0158 | 1.17e+03 | yes |
| patidx | Familias::FamiliasPosterior 2.6.4 |       200 | loci | 8.24e-05 | 0.178 | 2.15e+03 | yes |
| patidx | forrel::kinshipLR 1.10.0 |        20 | loci | 1.38e-05 | 0.0507 | 3.68e+03 | yes |
| patidx | forrel::kinshipLR 1.10.0 |       200 | loci | 8.43e-05 | 0.283 | 3.36e+03 | yes |
| mixlr | forensim::LR 4.3.3 |         4 | alleles | 8.11e-06 | 9.94e-05 | 12.3 | yes |
| mixlr | forensim::LR 4.3.3 |         6 | alleles | 1.14e-05 | 0.00033 | 29.1 | yes |
| mixlr | forensim::LR 4.3.3 |         8 | alleles | 1.9e-05 | 0.011 |  579 | yes |
| cllr | comparison::calc.ece 1.0.8 |     1,000 | ratios | 2.27e-05 | 0.0507 | 2.24e+03 | yes |
| cllr | comparison::calc.ece 1.0.8 |    10,000 | ratios | 0.000226 | 2.78 | 1.23e+04 | yes |
| cllr | base R |     1,000 | ratios | 2.24e-05 | 9.29e-05 | 4.14 | yes |
| cllr | base R |   100,000 | ratios | 0.00232 | 0.00961 | 4.15 | yes |
| cllr | base R | 1,000,000 | ratios | 0.027 | 0.098 | 3.63 | yes |
| cllrmin | comparison::calc.ece 1.0.8 |     1,000 | ratios | 0.00044 | 0.0504 |  115 | yes |
| cllrmin | comparison::calc.ece 1.0.8 |    10,000 | ratios | 0.00153 | 2.52 | 1.65e+03 | yes |
| lrcal | comparison::calibrate.set 1.0.8 |     1,000 | ratios | 0.000384 | 0.0459 |  120 | yes |
| lrcal | comparison::calibrate.set 1.0.8 |    10,000 | ratios | 0.00124 |    3 | 2.42e+03 | yes |
| eceplt | comparison::calc.ece 1.0.8 |     1,000 | ratios | 0.00966 | 0.068 | 7.04 | yes |
| eceplt | comparison::calc.ece 1.0.8 |    10,000 | ratios | 0.078 | 3.44 | 44.1 | yes |
| lrpost | base R |     1,000 | ratios | 8.13e-06 | 8.52e-06 | 1.05 | yes |
| lrpost | base R |   100,000 | ratios | 0.000753 | 0.00103 | 1.37 | yes |
| lrpost | base R | 1,000,000 | ratios | 0.00758 | 0.0184 | 2.43 | yes |
| dbsrch | base R |         1 | calls | 6e-06 | 2.76e-06 | 0.46 | yes |
| tipplt | base R |     1,000 | ratios | 0.00358 | 0.00417 | 1.16 | yes |
| tipplt | base R |   100,000 | ratios | 0.0453 | 0.0577 | 1.27 | yes |
| tipplt | base R | 1,000,000 | ratios | 0.355 | 0.435 | 1.22 | yes |
| clrsmp | base R |         5 | offence types | 7.72e-06 | 1.16e-05 |  1.5 | yes |
| clrsmp | base R |     1,000 | offence types | 4.88e-05 | 0.000176 | 3.61 | yes |
| clrsmp | base R |   100,000 | offence types | 0.00489 | 0.0219 | 4.48 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.4.3 |         5 | groups | 4.69e-06 | 0.00143 |  306 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.4.3 |     1,000 | groups | 3.41e-05 | 0.00144 | 42.3 | yes |
| expbnd | rmorie::morie_disparity_exposure_bounds 1.4.3 |   100,000 | groups | 0.00312 | 0.107 | 34.4 | yes |
| misrd | rmorie::morie_misclassified_rd 1.4.3 |         1 | calls | 4.26e-06 | 1.69e-05 | 3.98 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.4.3 |     1,000 | places | 1.56e-05 | 0.000275 | 17.6 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.4.3 |   100,000 | places | 0.00159 | 0.038 |   24 | yes |
| rtmhs | rmorie::morie_regression_to_mean 1.4.3 | 1,000,000 | places | 0.0228 | 0.316 | 13.9 | yes |
| incap (share) | rmorie::morie_incapacitation 1.4.3 |         5 | groups | 6.13e-06 | 0.000678 |  111 | yes |
| incap (share) | rmorie::morie_incapacitation 1.4.3 |     1,000 | groups | 0.000181 | 0.000973 | 5.39 | yes |
| incap (share) | rmorie::morie_incapacitation 1.4.3 |   100,000 | groups | 0.0352 | 0.0207 | 0.588 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.4.3 |         5 | groups | 1.02e-05 | 0.000803 | 78.7 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.4.3 |     1,000 | groups | 4.77e-05 | 0.000911 | 19.1 | yes |
| incap (no share) | rmorie::morie_incapacitation 1.4.3 |   100,000 | groups | 0.00817 | 0.0154 | 1.88 | yes |
| incapr | rmorie::morie_incapacitation_career 1.4.3 |        10 | career periods | 3.32e-06 | 3.41e-05 | 10.3 | yes |
| incapr | rmorie::morie_incapacitation_career 1.4.3 |     1,000 | career periods | 8.41e-06 | 9.77e-05 | 11.6 | yes |
| incapr | rmorie::morie_incapacitation_career 1.4.3 |   100,000 | career periods | 0.000537 | 0.00678 | 12.6 | yes |
| sellab | rmorie::morie_selective_labels 1.4.3 |     1,000 | defendants | 1.63e-05 | 0.000272 | 16.6 | yes |
| sellab | rmorie::morie_selective_labels 1.4.3 |   100,000 | defendants | 0.00192 | 0.0313 | 16.3 | yes |
| sellab | rmorie::morie_selective_labels 1.4.3 | 1,000,000 | defendants | 0.0299 | 0.445 | 14.9 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.4.3 |     1,000 | stays | 6.13e-06 | 4.8e-05 | 7.84 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.4.3 |   100,000 | stays | 0.000348 | 0.00388 | 11.2 | yes |
| losbnd | rmorie::morie_backlog_censoring 1.4.3 | 1,000,000 | stays | 0.00345 | 0.0372 | 10.8 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.4.3 |     1,000 | people | 4.99e-05 | 0.000451 | 9.05 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.4.3 |   100,000 | people | 0.00507 | 0.0327 | 6.45 | yes |
| sntbnd | rmorie::morie_sentence_effect_bounds 1.4.3 | 1,000,000 | people | 0.0467 |  0.3 | 6.43 | yes |
| hzdsel | rmorie::morie_hazard_selection 1.4.3 |         1 | calls | 4e-06 | 1.41e-05 | 3.53 | yes |
| proof | base R |         1 | calls | 1.1e-05 | 0.000126 | 11.4 | yes |
