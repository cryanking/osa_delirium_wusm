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


#include "cauchy_addon.h"


void heterdrmu(tree& t, xinfo& xi, dinfo_s& di, pinfo& pi, double *sigma, rn& gen)
{
  tree::npv bnv;
  std::vector<double> bv;
  std::vector<double> Mv;
  heterallsuff(t,xi,di,bnv,bv,Mv,sigma);
  if(di.using_cauchy) {
    for(tree::npv::size_type i=0;i!=bnv.size();i++)
      bnv[i]->settheta(heterdrawnodemu(bv[i]/di.cauchy_factor,Mv[i]/di.cauchy_factor,pi.tau,gen)*di.cauchy_factor );
  } else {
    for(tree::npv::size_type i=0;i!=bnv.size();i++)
      bnv[i]->settheta(heterdrawnodemu(bv[i],Mv[i],pi.tau,gen));
  }
}

void drmu(tree& t, xinfo& xi, dinfo_s& di, pinfo& pi, double sigma, rn& gen)
{
  tree::npv bnv;
  std::vector<size_t> nv;
  std::vector<double> syv;
  allsuff(t,xi,di,bnv,nv,syv);
  
  if(di.using_cauchy) {
    for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      bnv[i]->settheta(drawnodemu(nv[i],syv[i]/di.cauchy_factor,pi.tau,sigma/di.cauchy_factor,gen)*di.cauchy_factor);}
  } else {
    for(tree::npv::size_type i=0;i!=bnv.size();i++) 
    {bnv[i]->settheta(drawnodemu(nv[i],syv[i] ,pi.tau,sigma,gen));}
  }
}
