/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_cauchy_h
#define GUARD_cauchy_h

#include <Rcpp.h>
// [[Rcpp::depends(BART)]]
#include "BARTfiles/bart.h"
#include "BARTfiles/heterbart.h"
#include "BARTfiles/latent.h"
#include "BARTfiles/rtnorm.h"
// #include <bart.h>
// #include <heterbart.h>
// #include <latent.h>
// #include <rtnorm.h>




//data
class dinfo_s : public dinfo{
public:
  dinfo_s() {using_cauchy=false;cauchy_factor=1;}
  bool using_cauchy;
  double cauchy_factor;
};

class bart_s : public bart {
public:
  bart_s():bart() {}
  bart_s(size_t m):bart(m) {}
  void set_cauchy(double s) { di.using_cauchy =true; di.cauchy_factor=s;}
  dinfo_s di;
  double* get_allfit() {return this->allfit;}
  void set_allfit(double* s) {this->allfit = s;}
};

class heterbart_s : public heterbart
{
public:
  heterbart_s():heterbart() {}
  heterbart_s(size_t m):heterbart(m) {}
  void set_cauchy(double s) { di.using_cauchy =true; di.cauchy_factor=s;}
  dinfo_s di;
  double* get_allfit() {return this->allfit;}
  void set_allfit(double* s) {this->allfit = s;}
};

void drmu(tree& t, xinfo& xi, dinfo_s& di, pinfo& pi, double sigma, rn& gen);
void heterdrmu(tree& t, xinfo& xi, dinfo_s& di, pinfo& pi, double *sigma, rn& gen);


#endif
