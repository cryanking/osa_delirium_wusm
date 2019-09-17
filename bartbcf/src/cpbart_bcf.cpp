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

#include <ctime>

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "rtnorm.h"
#include "rand_draws.h"
#include "heterbart.h"

// #include <tree.h>
// #include <treefuns.h>
// #include <info.h>
// #include <bartfuns.h>
// #include <bd.h>
// #include <bart.h>
// #include <rtnorm.h>
// #include <rand_draws.h>
// #include <heterbart.h>

#include "cauchy_addon.h"


#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define TRDRAW_treated(a, b) trdraw_treated(a, b)
#define TEDRAW_treated(a, b) tedraw_treated(a, b)

RcppExport SEXP cpbart_causal_forest(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _inp,		//number of observations in test data
   SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,		//y, train,  nx1
   SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
  SEXP _num_treated, // first num_treated are exposed, rest unexposed
  SEXP _extra_var, //draw sigma independently in treated group?
  SEXP _use_cauchy, //use half-cauchy overparameterization?
   SEXP _im,		//number of trees
  SEXP _im_treated, //number of trees in treated function
   SEXP _inc,		//number of cut points
   SEXP _inc_treated,		//number of cut points
   SEXP _ind,		//number of kept draws (except for thinnning ..)
   SEXP _iburn,		//number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _binaryOffset,
   SEXP _itau,
  //   SEXP _iM, // number of shards for Modified LISA MCMC
   SEXP _idart,
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,
   SEXP _ib,
   SEXP _irho,
   SEXP _iaug,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
//   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo,
   SEXP _Xinfo_treated
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   size_t num_treated  = Rcpp::as<int>(_num_treated );
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::IntegerVector  yv(_iy); // binary
   int *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   size_t m_treated = Rcpp::as<int>(_im_treated);
   
   //size_t nc = Rcpp::as<int>(_inc);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   Rcpp::IntegerVector _nc_treated(_inc_treated);
   int *numcut_treated = &_nc_treated[0];
   
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double binaryOffset = Rcpp::as<double>(_binaryOffset);
   double tau = Rcpp::as<double>(_itau);
   bool extra_var = Rcpp::as<bool>(_extra_var);
   bool use_cauchy = Rcpp::as<bool>(_use_cauchy);
//   double rootM = sqrt(Rcpp::as<double>(_iM));
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
//   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   Rcpp::NumericMatrix Xinfo_treated(_Xinfo_treated);
   //Rcpp::List Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
/*
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
*/
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   
   Rcpp::NumericMatrix trdraw_treated(nkeeptrain,n); // hold the treatment effect estimates
   Rcpp::NumericMatrix tedraw_treated(nkeeptest,np);
   Rcpp::NumericMatrix varprb_treated(nkeeptreedraws,p); // hold the treatment effect counts
   Rcpp::IntegerMatrix varcnt_treated(nkeeptreedraws,p);
   
   //random number generation
   arn gen;

//    if(use_cauchy) {
     bart_s bm(m);
     bart_s bm_treated(m_treated);
//    } else {
//      bart bm(m);
//      bart bm_treated(m_treated);
//    }

     if(Xinfo.size()>0) {
       xinfo _xi;
       _xi.resize(p);
       for(size_t i=0;i<p;i++) {
         _xi[i].resize(numcut[i]);
         //Rcpp::IntegerVector cutpts(Xinfo[i]);
         for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
       }
       bm.setxinfo(_xi);
     }
     
     if(Xinfo_treated.size()>0) {
       xinfo _xi_treated;
       _xi_treated.resize(p);
       for(size_t i=0;i<p;i++) {
         _xi_treated[i].resize(numcut_treated[i]);
         for(size_t j=0;j<numcut_treated[i];j++) _xi_treated[i][j]=Xinfo_treated(i, j);
       }
       bm_treated.setxinfo(_xi_treated);
     }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
#define TRDRAW_treated(a, b) trdraw_treated[a][b]
#define TEDRAW_treated(a, b) tedraw_treated[a][b]

void cpbart_causal_forest(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   int* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
 size_t num_treated,		//number of records in the treated group
 bool extra_var, // heterogenous variance?
 bool use_cauchy,
 size_t m,		//number of trees
 size_t m_treated,		//number of trees
   int *numcut, //size_t nc,		//number of cut points
   int* numcut_treated,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double binaryOffset,
   double tau,
   bool dart,
   double theta,
   double omega,
   int *grp,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   //size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
/*
   double* trmean,
   double* temean,
*/
   double* _trdraw,
   double* _tedraw
)
{

   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;
   
   std::vector<double*> trdraw_treated(nkeeptrain);
   std::vector<double*> tedraw_treated(nkeeptest);
   
   
   for(size_t i=0; i<nkeeptrain; ++i) trdraw_treated[i]=&_trdraw_treated[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw_treated[i]=&_tedraw_treated[i*np];
   
   std::vector< std::vector<size_t> > varcnt_treated;
   std::vector< std::vector<double> > varprb_treated;   
   
   //random number generation
   arn gen(n1, n2);

//    if(use_cauchy) {
     bart_s bm(m);
     bart_s bm_treated(m_treated);
//    } else {
//      bart bm(m);
//      bart bm_treated(m_treated);
//    }
#endif

/*
   for(size_t i=0; i<n; ++i) trmean[i]=0.0;
   for(size_t i=0; i<np; ++i) temean[i]=0.0;
*/

   double* iz = new double[n];
   double* psuedo_xp = ix + (num_treated*p);
   int num_untreated = n - num_treated;
   
   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::stringstream treess_treated;  //string stream to write trees to  
   treess_treated.precision(10);
   treess_treated << nkeeptreedraws << " " << m_treated << " " << p << endl;
   
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   printf("*****Into main of pbart with causal forest mod\n");

   size_t skiptr,skipte,skiptreedraws;
   //size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
/*
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
*/
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %d, %d\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:mybeta,alpha,tau: %lf,%lf,%lf\n",
                   mybeta,alpha,tau);
   printf("*****binaryOffset: %lf\n",binaryOffset);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest,nkeeptreedraws: %zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptreedraws);
/*
   printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
*/
   printf("*****printevery: %zu\n",printevery);
   //printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);
   printf("*****skiptr,skipte,skiptreedraws: %zu,%zu,%zu\n",skiptr,skipte,skiptreedraws);

   //--------------------------------------------------
   //bart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,iz,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);
   
   bm_treated.setprior(alpha/3.0,mybeta+1.,tau);
   bm_treated.setdata(p,num_treated,ix,iz,numcut_treated); 
   bm_treated.setdart(a,b,rho,aug,dart,theta,omega);
   
   double *allfit_holder = bm_treated.get_allfit();
   bm_treated.set_allfit(bm.get_allfit() );
   
   //--------------------------------------------------
//init
  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(0., binaryOffset, 1., gen);
    else iz[k]=rtnorm(0., -binaryOffset, 1., gen);
  }

   //--------------------------------------------------

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   if(np) { fhattest = new double[np]; }
   double* fhatwork = new double[n];
   double restemp=0.0,rss=0.0;
   
   if(use_cauchy) {
     bm.set_cauchy(1.0);
   }
   
   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
//   size_t temecnt=0; //count test draws into posterior mean
//   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,/*keeptestme,*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();
   xinfo& xi_treated = bm_treated.getxinfo();
   
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(1., gen);
      //bm.draw(rootM, gen);

  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(-bm.f(k), binaryOffset, 1., gen);
    else iz[k]=rtnorm(bm.f(k), -binaryOffset, 1., gen);
  }

      if(i>=burn) {
        bm.predict(p,num_treated,ix,fhatwork); 
        
         //for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
           for(size_t k=num_treated;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
           for(size_t k=0;k<num_treated;k++) TRDRAW(trcnt,k)=fhatwork[k];
//             trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         //keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         //if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
//             tecnt+=1;
         }
/*
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
*/
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
//	   #ifndef NoRcpp
//	   Rcpp::List lists(m*treesaslists);
//	   #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
/*
	      #ifndef NoRcpp
	      varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
	      #endif
*/
	    }
	    #ifndef NoRcpp
//	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      //varcnt(i-burn,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	      //varprb(i-burn,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif

            //treedrawscnt +=1;
         }
      }
      // usual args for drawnodemu are sum(1/sigma^2), sum(y/sigma^2)
      // data are (y-untreated.predict)/ite which has sd^2 sigma^2/ite^2
      // so want to use sum [(y-untreated) * ite / sigma^2] and  sum(ite^2/sigma^2)
      // sigma is set to 1 for this update
      // Hahn actually recommends not using the cauchy update for the ite, just for the confounding variables, 
      
      if(use_cauchy) {
        //     if( i < burn)  bm.predict(p,num_treated,ix,fhatwork); 
        restemp =0. ; rss=0.;
        for(size_t k=0;k<n;k++) {
          restemp+= pow( (bm.f(k)  ) ,2 ) ;
          rss += (iz[k]  ) * (bm.f(k)  ) ;
          bm.set_cauchy(heterdrawnodemu(restemp,rss, tau, gen  ) ); // this is correct - the transformation induces heteroscekacitity in this update
        }
      }

      bm_treated.draw(1.0,gen);
      
      if(i>=burn) {

        bm_treated.predict(p,num_untreated,psuedo_xp,fhatwork);
        
    
        //for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
          //index = trcnt*n;;
          //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
          for(size_t k=0;k<num_treated;k++) TRDRAW_treated(trcnt,k)=bm_treated.f(k);
          for(size_t k=num_treated;k<n;k++) TRDRAW_treated(trcnt,k)=fhatwork[k-num_treated]+TRDRAW(trcnt,k);
          trcnt+=1;          
        }
        keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
        //keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
        //if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        if(keeptest) {
          bm_treated.predict(p,np,ixp,fhattest);
          //index=tecnt*np;
          //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
          for(size_t k=0;k<np;k++) TEDRAW_treated(tecnt,k)=fhattest[k]+TEDRAW(trcnt,k);
          tecnt+=1;
        }
        /*
         *  if(keeptestme) {
         *     for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
         *     temecnt+=1;
      }
      */
        keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
        if(keeptreedraw) {
          //	   #ifndef NoRcpp
          //	   Rcpp::List lists(m*treesaslists);
          //	   #endif
          
          for(size_t j=0;j<m_treated;j++) {
            treess_treated << bm_treated.gettree(j);
            /*
             * #ifndef NoRcpp
             * varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
             * if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
             * #endif
             */
          }
          #ifndef NoRcpp
          //	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
          ivarcnt=bm_treated.getnv();
          ivarprb=bm_treated.getpv();
          size_t k=(i-burn)/skiptreedraws;
          for(size_t j=0;j<p;j++){
            varcnt_treated(k,j)=ivarcnt[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb_treated(k,j)=ivarprb[j];
            //varprb(i-burn,j)=ivarprb[j];
          }
          #else
          varcnt_treated.push_back(bm.getnv());
          varprb_treated.push_back(bm.getpv());
          #endif
          
          //treedrawscnt +=1;
        }
      }
        
        
        
        
        
      
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
/*
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
*/
   printf("check counts\n");
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
   //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   delete[] iz;
   if(fhatwork) delete[] fhatwork;
   
  for(size_t j=0;j<nkeeptrain;j++) {
  for(size_t k=0; k<n; k++) {
      trdraw[j,k] += binaryOffset ;
      trdraw_treated[j,k] += binaryOffset ;
  }
  }
  for(size_t j=0;j<nkeeptest;j++) {
  for(size_t k=0; k<np; k++) {
      tedraw[j,k] += binaryOffset ;
      tedraw_treated[j,k] += binaryOffset ;
  }
  }
  

   bm_treated.set_allfit(allfit_holder );
   allfit_holder =0;

#ifndef NoRcpp
   //--------------------------------------------------
   //return
   Rcpp::List ret;
   //ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   //ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;
   ret["yhat.train.treated"]=trdraw_treated;
   ret["yhat.test.treated"]=tedraw_treated;
   ret["varcount.treated"]=varcnt_treated;
   ret["varprob.treated"]=varprb_treated;
   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   //if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   Rcpp::List xiret_treated(xi_treated.size());
   for(size_t i=0;i<xi_treated.size();i++) {
     Rcpp::NumericVector vtemp(xi_treated[i].size());
     std::copy(xi_treated[i].begin(),xi_treated[i].end(),vtemp.begin());
     xiret_treated[i] = Rcpp::NumericVector(vtemp);
   }
   
   Rcpp::List treesL_treated;
   treesL_treated["cutpoints"] = xiret_treated;
   treesL_treated["trees"]=Rcpp::CharacterVector(treess_treated.str());
   ret["treedraws.treated"] = treesL_treated;
   
   return ret;
#else

#endif

}
