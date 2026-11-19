// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// rOpenSci statistical software standards for spatial software, read by
// the srr package (srr::srr_report()), as src/srr.h holds the General and
// EDA ones. Standards met in one place are tagged there (spin.cpp,
// bsp.cpp, test-sp.R); the rest are here. Nothing here is compiled: no
// source file includes this header.
#ifndef CLAB_SRRSP_H
#define CLAB_SRRSP_H

//' srr_stats
//'
//' Spatial standards the package meets.
//'
//' @srrstats {SP1.0} Each spatial help page says whether it takes planar
//'   coordinates, longitude and latitude, or either through `lonlat`; the
//'   planar kernels assume Cartesian coordinates in one length unit, and
//'   longitude and latitude reach them through the equal-area projection of
//'   eqproj() (README, "Spatial input").
//' @srrstats {SP1.1} Every spatial function is two-dimensional; sf points
//'   with a third coordinate are refused.
//' @srrstats {SP2.2} sf objects are taken as input, and plain results
//'   (data frames of coordinates) go back to sf with sf::st_as_sf().
//' @srrstats {SP2.2a} README "Spatial input" shows sf points going in and
//'   a result going back to sf.
//' @srrstats {SP2.3} test-sp.R writes points to a GeoPackage, reads them
//'   back with sf::st_read() and passes them in.
//' @srrstats {SP2.6} Each spatial help page documents the forms its points
//'   may take (vectors, a matrix, a data frame or sf points) and their
//'   units.
//' @srrstats {SP2.9} Points keep their input order through every function,
//'   and the CRS of sf input decides whether they are degrees.
//' @srrstats {SP4.0} Results are base R data frames, named vectors and
//'   lists, as for the rest of the package (SP4.0b).
//' @srrstats {SP4.0b} Each help page documents the result's form, which
//'   is the same whatever form the input took.
//' @srrstats {SP4.1} Lengths come back in the unit of the coordinates
//'   (metres for longitude and latitude, and the CRS's unit for projected
//'   sf input) and areas in its square, as each help page states; a units
//'   class on an input column is not carried to the plain numbers returned.
//' @srrstats {SP4.2} Each help page gives the type and class of the
//'   result; test-sp.R checks them.
//' @srrstats {SP3.0} Neighbourhoods are the user's: a radius (fzmode(),
//'   stacel()), a threshold distance fixed or from the random expectation
//'   at a chosen probability (nnhclu()), a neighbour order (nnidx()), and
//'   for the autocorrelation statistics k nearest neighbours, a distance
//'   band or a contiguity list (spwt()).
//' @srrstats {SP3.0a} On grids, gstr()'s radius 1 gives the rook and 1.5
//'   the queen neighbourhood (cell included); spwt() takes any contiguity
//'   list, rook or queen.
//' @srrstats {SP3.0b} In irregular space neighbours are set by an integer
//'   number (nnidx()'s `kmax`) or a distance (fzmode()'s `r`, nnhclu()'s
//'   `fixd`, stacel()'s `r`).
//' @srrstats {SP3.1} jtsurf() weights each incident by a decay function of
//'   its distance (six forms, or a calibrated table from jtcal()), and
//'   spwt() gives inverse-distance weights of any power, not only a uniform
//'   cut-off.
//' @srrstats {SP3.4} nnhclu() and stacel() cluster by explicitly spatial
//'   rules (distance thresholds, search circles); kmpart() is K-means on
//'   the coordinates alone, where proximity is the only variable, not an
//'   added weight.
//' @srrstats {SP3.3} sarml() estimates the spatial autoregressive
//'   coefficient of the outcome (`rho`) and that of the errors (`lambda`)
//'   apart from the coefficients of the covariates, and spill() separates
//'   each covariate's direct effect from its spillover.
//' @srrstats {SP5.0} The spatial results that are drawn have plot
//'   methods (help page `plotm`): kernel and prospective surfaces, scan
//'   clusters, Knox tables, quadrat counts, hot spots over periods, the
//'   cross K function, buffers, thin-plate surfaces and DBSCAN clusters.
//' @srrstats {SP5.1} Those plots put the first coordinate on the
//'   horizontal axis and the second on the vertical one.
//' @srrstats {SP6.3} test-sp.R checks each way of defining neighbours
//'   against brute force: radii 0 and positive (fzmode()), fixed and
//'   random-expectation thresholds (nnhclu()), orders 1 to 3 (nnidx()),
//'   search circles (stacel()).
//' @srrstats {SP6.4} test-sp.R checks jtsurf() under the decay forms and
//'   both distances (Euclidean, Manhattan) against the surfaces written out
//'   in R, and each decay form against its closed form.
//' @noRd
#define CLAB_SRR_SP 1

//' NA_standards
//'
//' Spatial standards that do not apply, with the reason for each.
//'
//' @srrstatsNA {SP2.0a} No new class for spatial data is implemented:
//'   points come in as vectors, a matrix, a data frame or sf points. The
//'   classes that results carry are for print() and plot() only.
//' @srrstatsNA {SP2.5} No new class for spatial data is implemented
//'   (SP2.0a); sf input brings its own CRS, which is read (SP2.4).
//' @srrstatsNA {SP2.5a} No new class for spatial data is implemented
//'   (SP2.0a).
//' @srrstatsNA {SP3.2} Simulations draw complete spatial randomness, the
//'   null hypothesis being tested; sampling by local density would change
//'   the null, not estimate it.
//' @srrstatsNA {SP3.5} No spatial machine learning: the learners
//'   (enpath(), rfor(), dmlx()) take no coordinates.
//' @srrstatsNA {SP3.6} No spatial machine learning (SP3.5).
//' @srrstatsNA {SP4.0a} Results are base data frames and vectors whatever
//'   the input class (SP4.0b).
//' @srrstatsNA {SP5.2} Planar coordinates come in the caller's unit,
//'   which no argument states, so the axes of a plot are labelled x and
//'   y, without a unit.
//' @srrstatsNA {SP5.3} The plots use base graphics; the results are
//'   plain tables and matrices, which an interactive library takes as
//'   they are.
//' @srrstatsNA {SP6.6} No spatial machine learning (SP3.5).
//' @noRd
#endif
