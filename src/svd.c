#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void lcr(double* theta, int* L, int* n, double* lr)
{
  double lc[*L]; 
  lc[0] = 0.0;
  double ls[*L];
  double tmp;
  int l,k;
  lr[0] = 0.0;
  
  for (l = 1 ; l <= *L ; ++l) 
  {
    tmp=0.0; 

    for(k = 0; k < *n; ++k) 
    { 
      tmp = tmp + pow(theta[k] / theta[0], l); 
    }

    ls[l] = l * log(theta[0]) + log(tmp);
    tmp=0.0; 
    for (k = 0 ; k < l ; ++k)
    {
      tmp = tmp + exp( lc[k] + ls[l - k] - ls[l] -log(2.0 * l) ) ; 
    }
    lc[l] = ls[l] + log(tmp);
    lr[l] = lc[l] + lgamma(l + 1) + lgamma(*n / 2.0) - lgamma(*n / 2.0 + l) ;
  }
}

void ln2moment(double* mu, double* s2 , int* mmax, double* l2mom)
{
  int i;
  double lmom[*mmax * 2 + 1];
  l2mom[0] = 0;
  lmom[0] = 0;
  lmom[1] = log(*mu);
  
  for (i = 2; i <= *mmax * 2; ++i) {
    lmom[i] = lmom[i - 1] + 
      log( (i - 1) * *s2 * exp(lmom[i - 2] - lmom[i - 1]) + *mu);
    if(i % 2 == 0) 
    { 
      l2mom[i / 2] = lmom[i]; 
    }
  }
}

void rW(double* kap, int* m, double* W)
{
  
  GetRNGstate();
  
  double b = (-2.0 * *kap + sqrt(4 * pow(*kap, 2) + pow(*m - 1.0, 2)) ) /
    (*m - 1.0);
  
  double x0 = (1.0 - b) / (1.0 + b);
  double c = *kap * x0 + (*m - 1.0) * log(1.0 - pow(x0, 2));
  
  double Z, U;
  int done = 0;

  while (done == 0) {
    
    Z = rbeta((*m - 1.0) / 2.0, (*m - 1.0) / 2.0);
    *W = (1 - (1 + b) * Z) / (1.0 - (1.0 - b) * Z);
    U = runif(0, 1);
    if(*kap * *W + (*m - 1) * log(1 - x0 * *W) - c > log(U)) 
    { 
      done=1;
    }
  }

  PutRNGstate();
  
}

