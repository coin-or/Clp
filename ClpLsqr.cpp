// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpLsqr.hpp"
#include "ClpPdco.hpp"


void ClpLsqr::do_lsqr( CoinDenseVector &b,  
		real damp, real atol, real btol, real conlim, int itnlim, 
		bool show, Info info, CoinDenseVector &x , int *istop,
		int *itn, Outfo *outfo, bool precon, CoinDenseVector &Pr){

 /**
 Special version of LSQR for use with pdco.m.
 It continues with a reduced atol if a pdco-specific test isn't
 satisfied with the input atol.
*/

//     Initialize.
  static char term_msg[8][80] = 
  {
    "The exact solution is x = 0",
    "The residual Ax - b is small enough, given ATOL and BTOL",
    "The least squares error is small enough, given ATOL",
    "The estimated condition number has exceeded CONLIM",
    "The residual Ax - b is small enough, given machine precision",
    "The least squares error is small enough, given machine precision",
    "The estimated condition number has exceeded machine precision",
    "The iteration limit has been reached"
  };

  //  printf("***************** Entering LSQR *************\n");
  assert (model_);
  
  char str1[100], str2[100], str3[100], str4[100], head1[100], head2[100];

  int m = nrows_;  int n = ncols_;     // set  m,n from lsqr object
  
  *itn     = 0;            *istop   = 0;       int nstop  = 0;
  real ctol   = 0;         
  if (conlim > 0) ctol = 1/conlim;

  real anorm  = 0;         real acond  = 0;
  real dampsq = damp*damp; real ddnorm = 0;      real res2   = 0;
  real xnorm  = 0;         real xxnorm = 0;      real z      = 0;
  real cs2    = -1;        real sn2    = 0;

  // Set up the first vectors u and v for the bidiagonalization.
  // These satisfy  beta*u = b,  alfa*v = A'u.

  CoinDenseVector u(b);
  CoinDenseVector v(n, 0.0);
  x.clear();
  real alfa   = 0;      real beta = u.twoNorm();
  if (beta > 0){
    u = (1/beta) * u;
    matVecMult( 2, v, u );
    if (precon)
      v = v*Pr;
    alfa = v.twoNorm();
  }
  if (alfa > 0){
    v.scale(1/alfa);
  }
  CoinDenseVector w(v);

  real arnorm = alfa * beta; 
  if (arnorm == 0){
     printf("  %s\n\n", term_msg[0]);
     return;
  }

  real rhobar = alfa;  real phibar = beta;   real bnorm  = beta;
  real rnorm  = beta;
  sprintf(head1,"   Itn      x(1)      Function");
  sprintf(head2," Compatible   LS      Norm A   Cond A");

  if (show){
    printf(" %s%s\n",head1, head2);
    real test1  = 1;  real test2  = alfa / beta;
    sprintf(str1,"%6d %12.5e %10.3e",   *itn, x[0], rnorm );
    sprintf(str2,"  %8.1e  %8.1e",       test1, test2 );
    printf("%s%s\n",str1, str2);
  }

  //----------------------------------------------------------------
  // Main iteration loop.
  //----------------------------------------------------------------
  while (*itn < itnlim){
    *itn += 1;
    // Perform the next step of the bidiagonalization to obtain the
    // next beta, u, alfa, v.  These satisfy the relations
    // beta*u  =  a*v   -  alfa*u,
    // alfa*v  =  A'*u  -  beta*v.

    u.scale((-alfa));
    if (precon){
      CoinDenseVector pv(v*Pr);
      matVecMult( 1, u, pv);
    }else{
      matVecMult( 1, u, v);
    }
    beta = u.twoNorm();
    if (beta > 0){
      u.scale((1/beta));
      anorm = sqrt(anorm*anorm + alfa*alfa + beta*beta +  damp*damp);
      v.scale((-beta));
      CoinDenseVector vv(n);
      vv.clear();
      matVecMult( 2, vv, u );
      if (precon)
	vv = vv*Pr;
      v = v + vv;
      alfa  = v.twoNorm();
      if (alfa > 0)
	v.scale((1/alfa));
    }

    // Use a plane rotation to eliminate the damping parameter.
    // This alters the diagonal (rhobar) of the lower-bidiagonal matrix.

    real rhobar1 = sqrt(rhobar*rhobar + damp*damp);
    real cs1     = rhobar / rhobar1;
    real sn1     = damp   / rhobar1;
    real psi     = sn1 * phibar;
    phibar       *= cs1;

    // Use a plane rotation to eliminate the subdiagonal element (beta)
    // of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

    real rho     = sqrt(rhobar1*rhobar1 +  beta*beta);
    real cs      =   rhobar1/ rho;
    real sn      =   beta   / rho;
    real theta   =   sn * alfa;
    rhobar  = - cs * alfa;
    real phi     =   cs * phibar;
    phibar  =   sn * phibar;
    real tau     =   sn * phi;

    // Update x and w.

    real t1      =   phi  /rho;
    real t2      = - theta/rho;
    //    dk           =   ((1/rho)*w);

    real w_norm = w.twoNorm();
    x       = x      +  t1*w;
    w       = v      +  t2*w;
    ddnorm  = ddnorm +  (w_norm/rho)*(w_norm/rho);
    // if wantvar, var = var  +  dk.*dk; end

    // Use a plane rotation on the right to eliminate the
    // super-diagonal element (theta) of the upper-bidiagonal matrix.
    // Then use the result to estimate  norm(x).

    real delta   =   sn2 * rho;
    real gambar  = - cs2 * rho;
    real rhs     =   phi  -  delta * z;
    real zbar    =   rhs / gambar;
    xnorm   =   sqrt(xxnorm + zbar*zbar);
    real gamma   =   sqrt(gambar*gambar + theta*theta);
    cs2     =   gambar / gamma;
    sn2     =   theta  / gamma;
    z       =   rhs    / gamma;
    xxnorm  =   xxnorm  +  z*z;

    // Test for convergence.
    // First, estimate the condition of the matrix  Abar,
    // and the norms of  rbar  and  Abar'rbar.

    acond   =   anorm * sqrt( ddnorm );
    real res1    =   phibar*phibar;
    real res2    =   res1  +  psi*psi;
    rnorm   =   sqrt( res1 + res2 );
    arnorm  =   alfa * fabs( tau );

    // Now use these norms to estimate certain other quantities,
    // some of which will be small near a solution.

    real test1   =   rnorm / bnorm;
    real test2   =   arnorm/( anorm * rnorm );
    real test3   =       1 / acond;
    t1      =   test1 / (1    +  anorm * xnorm / bnorm);
    real rtol    =   btol  +  atol *  anorm * xnorm / bnorm;

    // The following tests guard against extremely small values of
    // atol, btol  or  ctol.  (The user may have set any or all of
    // the parameters  atol, btol, conlim  to 0.)
    // The effect is equivalent to the normal tests using
    // atol = eps,  btol = eps,  conlim = 1/eps.

    if (*itn >= itnlim)
      *istop = 7; 
    if (1 + test3  <= 1)
      *istop = 6;
    if (1 + test2  <= 1)
      *istop = 5;
    if (1 + t1     <= 1)
      *istop = 4;

    // Allow for tolerances set by the user.

    if  (test3 <= ctol)
      *istop = 3;
    if  (test2 <= atol)
      *istop = 2;
    if  (test1 <= rtol)
      *istop = 1;

    //-------------------------------------------------------------------
    // SPECIAL TEST THAT DEPENDS ON pdco.m.
    // Aname in pdco   is  iw in lsqr.
    // dy              is  x
    // Other stuff     is in info.

    // We allow for diagonal preconditioning in pdDDD3.
    //-------------------------------------------------------------------
    if (*istop > 0){
      real r3new     = arnorm;
      real r3ratio   = r3new / info.r3norm;
      real atolold   = atol;
      real atolnew   = atol;
         
      if (atol > info.atolmin){
        if     (r3ratio <= 0.1){   // dy seems good
          // Relax
        }else if (r3ratio <= 0.5){ // Accept dy but make next one more accurate.
          atolnew = atolnew * 0.1;
        }else{                     // Recompute dy more accurately
	  if (show){
	    printf("\n                                ");
	    printf("                                \n");
	    printf(" %5.1f%7d%7.3f", log10(atolold), *itn, r3ratio);
	  }
          atol    = atol * 0.1;
          atolnew = atol;
          *istop   = 0;
        }
    
      outfo->atolold = atolold;
      outfo->atolnew = atolnew;
      outfo->r3ratio = r3ratio;
    }
      
    //-------------------------------------------------------------------
    // See if it is time to print something.
    //-------------------------------------------------------------------
    int prnt = 0;
    if (n     <= 40       ) prnt = 1;
    if (*itn   <= 10       ) prnt = 1;
    if (*itn   >= itnlim-10) prnt = 1;
    if (fmod(*itn,10) == 0  ) prnt = 1;
    if (test3 <=  2*ctol  ) prnt = 1;
    if (test2 <= 10*atol  ) prnt = 1;
    if (test1 <= 10*rtol  ) prnt = 1;
    if (*istop !=  0       ) prnt = 1;

    if (prnt == 1){
      if (show){
        sprintf(str1, "   %6d %12.5e %10.3e",   *itn, x[0], rnorm );
        sprintf(str2, "  %8.1e %8.1e",       test1, test2 );
        sprintf(str3, " %8.1e %8.1e",        anorm, acond );
        printf("%s%s%s\n", str1, str2, str3);
      }
    }
    if (*istop > 0)
      break;
    }
  }
  // End of iteration loop.
  // Print the stopping condition.

  if (show){
    printf("\n LSQR finished\n");
  //    disp(msg(istop+1,:))
  //    disp(' ')
    printf("%s\n", term_msg[*istop]);
    sprintf(str1, "istop  =%8d     itn    =%8d",      *istop, *itn    );
    sprintf(str2, "anorm  =%8.1e   acond  =%8.1e",  anorm, acond  );
    sprintf(str3, "rnorm  =%8.1e   arnorm =%8.1e",  rnorm, arnorm );
    sprintf(str4, "bnorm  =%8.1e   xnorm  =%8.1e",  bnorm, xnorm  );
    printf("%s %s\n", str1, str2);
    printf("%s %s\n", str3, str4);
  }
}




void ClpLsqr::matVecMult( int mode, CoinDenseVector *x, CoinDenseVector *y){
  int n = model_->numberColumns();
  int m = model_->numberRows();
  CoinDenseVector *temp = new CoinDenseVector(n, 0.0);
  real *t_elts = temp->getElements();
  real *x_elts = x->getElements();
  real *y_elts = y->getElements();
  ClpPdco * pdcoModel = (ClpPdco *) model_;
  if (mode == 1){
    pdcoModel->matVecMult( 2, temp, y);
    for (int k=0; k<n; k++)
      x_elts[k] += (diag1_[k] * t_elts[k]);
    for (int k=0; k<m; k++)
      x_elts[n+k] += (diag2_ * y_elts[k]);
  }else{
    for (int k=0; k<n; k++)
      t_elts[k] = diag1_[k] * y_elts[k];
    pdcoModel->matVecMult( 1, x, temp);
    for (int k=0; k<m; k++)
      x_elts[k] += diag2_ * y_elts[n+k];
  }
  delete temp;
  return;
}

void ClpLsqr::matVecMult( int mode, CoinDenseVector &x, CoinDenseVector &y){
  matVecMult( mode, &x, &y);
  return;
}
/* Default constructor */
ClpLsqr::ClpLsqr() :
  nrows_(0),
  ncols_(0),
  model_(NULL),
  diag1_(NULL),
  diag2_(0.0)
{}

/* Constructor for use with Pdco model (note modified for pdco!!!!) */
ClpLsqr::ClpLsqr(ClpInterior *model) :
  diag1_(NULL),
  diag2_(0.0)
{
  model_ = model;
  nrows_ = model->numberRows() + model->numberColumns();
  ncols_ = model->numberRows();
}
/** Destructor */
ClpLsqr::~ClpLsqr()
{
  // delete [] diag1_; no as we just borrowed it
}
bool 
ClpLsqr::setParam(char *parmName, int parmValue){
  std::cout << "Set lsqr integer parameter " << parmName << "to " << parmValue
	    << std::endl;
  if ( strcmp(parmName, "nrows") == 0) {
    nrows_ = parmValue;
    return 1;
  }else if ( strcmp(parmName, "ncols") == 0) {
    ncols_ = parmValue;
    return 1;
  }
  std::cout <<"Attempt to set unknown integer parameter name "<<parmName<<std::endl;
  return 0;
}
ClpLsqr::ClpLsqr(const ClpLsqr &rhs) :
  nrows_(rhs.nrows_),
  ncols_(rhs.ncols_),
  model_(rhs.model_),
  diag2_(rhs.diag2_)
{
  diag1_ = ClpCopyOfArray(rhs.diag1_,nrows_);
}
// Assignment operator. This copies the data
ClpLsqr & 
ClpLsqr::operator=(const ClpLsqr & rhs)
{
  if (this != &rhs) {
    delete [] diag1_;
    diag1_ = ClpCopyOfArray(rhs.diag1_,nrows_);
    nrows_ = rhs.nrows_;
    ncols_ = rhs.ncols_;
    model_ = rhs.model_;
    diag2_ = rhs.diag2_;
  }
  return *this;
}
