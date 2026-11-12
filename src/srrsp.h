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
//'   at a chosen probability (nnhclu()), a neighbour order (nnidx()).
//' @srrstats {SP3.0b} In irregular space neighbours are set by an integer
//'   number (nnidx()'s `kmax`) or a distance (fzmode()'s `r`, nnhclu()'s
//'   `fixd`, stacel()'s `r`).
//' @srrstats {SP3.1} jtsurf() weights each incident by a decay function of
//'   its distance (six forms, or a calibrated table from jtcal()), not by a
//'   uniform cut-off.
//' @srrstats {SP3.4} nnhclu() and stacel() cluster by explicitly spatial
//'   rules (distance thresholds, search circles); kmpart() is K-means on
//'   the coordinates alone, where proximity is the only variable, not an
//'   added weight.
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
//' @srrstatsNA {SP2.0a} No new class is implemented.
//' @srrstatsNA {SP2.5} No new class is implemented; sf input brings its
//'   own CRS, which is read (SP2.4).
//' @srrstatsNA {SP2.5a} No new class is implemented.
//' @srrstatsNA {SP3.0a} No function considers neighbours on a grid.
//' @srrstatsNA {SP3.2} Simulations draw complete spatial randomness, the
//'   null hypothesis being tested; sampling by local density would change
//'   the null, not estimate it.
//' @srrstatsNA {SP3.3} No spatial regression.
//' @srrstatsNA {SP3.5} No machine learning.
//' @srrstatsNA {SP3.6} No machine learning.
//' @srrstatsNA {SP4.0a} Results are base data frames and vectors whatever
//'   the input class (SP4.0b).
//' @srrstatsNA {SP5.0} No class is implemented, so no plot method.
//' @srrstatsNA {SP5.1} No spatial function draws a plot.
//' @srrstatsNA {SP5.2} No spatial function draws a plot.
//' @srrstatsNA {SP5.3} No spatial function draws a plot.
//' @srrstatsNA {SP6.6} No machine learning.
//' @noRd
#endif
