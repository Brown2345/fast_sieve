#include <stdlib.h>
#include <stdio.h>
#include <gmpxx.h>
#include <math.h>
#include <sys/time.h>
#include "diophappr.h"
#include "simplesieve.h"
#include <omp.h>
#include <algorithm>    // std::min

/*  Compile using
g++ sievesimp.cpp -osievesimp -O3 -frounding-math -mfpmath=387 -finline-functions -lgmp -fopenmp */
double tzero = 0, tone = 0, ttwo = 0;
long counter = 0;
void OrigSegSiev(short *s, unsigned long n, long D, long K, long quot)
/* ensure: s[j] = (n-D+j is prime) for 0<=j<=2 D */
/* parameter kappa = 1/quot */
{
  mpz_class sqt, M2kap, sqtR, kmpz;
  mpz_class npD(n+D);
  long M, Mp, R, m, m0, x;
  unsigned long np;
  rat alpha0, alpha1;
  mpq_class eta, etaq;
  long c,cainv,a,ainv,q,k,r0,j;
  
  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  x = sqt.get_si();                     /* x = (int) sqrt(n+D) */
  		  
  SubSegSiev(s,n-D,2*D,K*D);
  for(M = K*D+1; M <=x; M+=(2*R+1)) {
     M2kap = ((((((mpz_class) M)*M)*D)/quot)/n);  /* (M^2)*D/(quot*n) */
    mpz_sqrt(sqtR.get_mpz_t(), M2kap.get_mpz_t());
    R = sqtR.get_si();  /* (int) (M*sqrt(kappa*D/n)) */
    
    m0 = M+R; Mp = M+2*R;
    alpha0.num = n%m0; alpha0.den = m0;
    alpha1.den = m0*m0; alpha1.num = (alpha1.den-n%alpha1.den);
    eta = (((mpq_class) D)/((mpq_class) M))*(1+1/((mpq_class) quot)); /* eta = 3D/2M */

    diophappr(alpha1,2*R,&a,&ainv,&q);
    etaq = eta*q;
    kmpz = etaq.get_num()/etaq.get_den(); k=kmpz.get_si();
    c = (alpha0.num*q+m0/2)/m0;
    cainv  = (ainv*c)%q;


    for(r0= -cainv, j=0; j <= k+1; j++, r0 -= ainv) {
      if(r0 <= -q)
	r0 += q;
      for(m= m0+r0; m>=M; m-=q)
	if(m%2) {
	  np = ((n+D)/m)*m;
	  if(np>=n-D && np<=n+D && np>m) 
	    s[np-(n-D)]=0;		
	}
      for(m=m0+r0+q; m<=Mp; m+=q)
	if(m%2) {
	  np = ((n+D)/m)*m;
	  if(np>=n-D && np<=n+D && np>m) 
	    s[np-(n-D)]=0;	
	}
    }

    for(r0= -cainv+ainv, j = -1; j >= -(k+1); j--, r0+=ainv) {
      if(r0>0)
	r0 -=q;
      for(m= m0+r0; m>=M; m-=q)
	if(m%2) {
	  np = ((n+D)/m)*m;
	  if(np>=n-D && np<=n+D && np>m) 
	    s[np-(n-D)]=0;
	}
      for(m=m0+r0+q; m<=Mp; m+=q)
	if(m%2) {
	  np = ((n+D)/m)*m;
	  if(np>=n-D && np<=n+D && np>m) 
	    s[np-(n-D)]=0;
	}
    }
  }
}


void NewSegSiev(short *s, unsigned long n, long D, long K, long quot)
/* ensure: s[j] = (n-D+j is prime) for 0<=j<=2 D */
/* parameter kappa = 1/quot */
{
  mpz_class sqt, M2kap, sqtR, kmpz;
  mpz_class npD(n+D);
  long M, Mp, R, m, m0, x;
  unsigned long np;
  rat alpha0, alpha1;
  mpq_class eta, etaq;
  long c,cainv,a,ainv,q,k,r0;
  
  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  x = sqt.get_si();                     /* x = (int) sqrt(n+D) */
  		  
  SubSegSiev(s,n-D,2*D,K*D);
  for(M = K*D+1; M <=x; M+=(2*R+1)) { //unabh: M, K, D, R
     M2kap = ((((((mpz_class) M)*M)*D)/quot)/n);  /* (M^2)*D/(quot*n) */
    mpz_sqrt(sqtR.get_mpz_t(), M2kap.get_mpz_t());
    R = sqtR.get_si();  /* (int) (M*sqrt(kappa*D/n)) */
    
    m0 = M+R; Mp = M+2*R;
    alpha0.num = n%m0; alpha0.den = m0;
    alpha1.den = m0*m0; alpha1.num = (alpha1.den-n%alpha1.den);
    eta = (((mpq_class) D)/((mpq_class) M))*(1+1/((mpq_class) quot)); /* eta = 3D/2M */

    diophappr(alpha1,2*R,&a,&ainv,&q);
    etaq = eta*q;
    kmpz = etaq.get_num()/etaq.get_den(); k=kmpz.get_si();
    c = (alpha0.num*q+m0/2)/m0;
    cainv  = (ainv*c)%q;

		long chunksize = ceil((1.0*(k+2)) / omp_get_max_threads());
		#pragma omp parallel // private (np, r0)
		{

		int tid = omp_get_thread_num();
		r0 = -cainv;
		r0 -= tid * chunksize * ainv;
		r0 += (abs(r0)/q) * q;

		long start = tid * chunksize;	
		long stop  = (tid+1) * chunksize;
		if (stop > k+2) stop = k+2;

		//#pragma omp for nowait
    for(long j=start; j < stop; j++, r0-=ainv) {
      //r0 = -cainv - j * ainv;	
      //r0 += (abs(r0)/q) * q;
      if (r0+q <= 0)
				r0 += q;
      
			for(m= m0+r0; m>=M; m-=q)
				if(m%2) {
					np = ((n+D)/m)*m;
					if(np>=n-D && np<=n+D && np>m) 
			  		s[np-(n-D)]=0;		
				}
      for(m=m0+r0+q; m<=Mp; m+=q)
				if(m%2) {
	  			np = ((n+D)/m)*m;
	  			if(np>=n-D && np<=n+D && np>m) 
	    			s[np-(n-D)]=0;	
				}
    }

		r0 = -cainv + ainv;
		r0 += tid * chunksize * ainv;
    r0 -= (abs(r0)/q) * q;

		start = tid * chunksize;	
		stop  = (tid+1) * chunksize;
		if (stop > k+1) stop = k+1;
		
    //#pragma omp for nowait
		for(long j = start; j < stop; j++, r0+=ainv) {
			//r0 = -cainv + (j+1) * ainv;
      //r0 -= (abs(r0)/q) * q;
      if (r0 > 0)
    		r0 -= q;
      for(m= m0+r0; m>=M; m-=q)
				if(m%2) {
					np = ((n+D)/m)*m;
					if(np>=n-D && np<=n+D && np>m) 
	    			s[np-(n-D)]=0;
				}
      for(m=m0+r0+q; m<=Mp; m+=q)
				if(m%2) {
	  			np = ((n+D)/m)*m;
	  			if(np>=n-D && np<=n+D && np>m) 
	    			s[np-(n-D)]=0;
				}
		}


		} // parallel part end
  }
}

int main(int argc, char *argv[])
{
  int print = 0;
  int test = 0;
  int errors = 0;
  int nthreads = 4;
  unsigned long n, D, j, n1, D1, n2, D2;
  short *s1, *s2; //manu
  timeval tstart, tend;

//./sievesimp 5000000000000000000 10000000 4 

  
  if(argc<5) {
    n = 5000000000000000000;
    D = 40000000;
    nthreads = 4;
    test = 0;
  } else {
    n = atol(argv[1]);
    D = atol(argv[2]);
    nthreads = atol(argv[3]);
    test = atol(argv[4]);
  }

  s2 = (short *) calloc(2*D+1,sizeof(short));
  n1 = n2 = n;
  D1 = D2 = D;
  if (test) {
    s1 = (short *) calloc(2*D+1,sizeof(short));
    omp_set_num_threads(1);
    OrigSegSiev(s1,n1,D1,8,4);
  }
  omp_set_num_threads(nthreads);

  gettimeofday(&tstart,0);
  NewSegSiev(s2,n2,D2,8,4);
  gettimeofday(&tend,0);

  if (print)
    for(j=0; j<=2*D; j++)
	if (s2[j])
	    printf("%lu\n",n-D+j);
  if (test)
    for(j=0; j<=2*D; j++)
     if(s1[j] != s2[j])
  	errors++;
  if (test)
    printf("#errors: %i\n", errors);

  long seconds = tend.tv_sec-tstart.tv_sec;
  long useconds = tend.tv_usec-tstart.tv_usec;

  fprintf(stderr,"Seconds: %lf\n",seconds + useconds/1000000.0);
  //fprintf(stderr,"tzero: %f tone: %f ttwoL %f \n",tzero, tone, ttwo);
}
