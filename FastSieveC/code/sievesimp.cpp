#include <stdlib.h>
#include <stdio.h>
#include <gmpxx.h>
#include <math.h>
#include <sys/time.h>
#include "diophappr.h"
#include "simplesieve.h"
#include <omp.h>

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
  timeval t1, t2, t4;
    gettimeofday(&t1,0);


    SubSegSiev(s,n-D,2*D,K*D);

    gettimeofday(&t2,0);
    long seconds = t2.tv_sec-t1.tv_sec;
    long useconds= t2.tv_usec-t1.tv_usec;
    fprintf(stderr,"Seconds SubSegSiev: %lf\n",seconds + useconds/1000000.0);

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
  gettimeofday(&t4,0);
  seconds = t4.tv_sec-t2.tv_sec;
  useconds= t4.tv_usec-t2.tv_usec;
  fprintf(stderr,"Seconds 3: %lf\n",seconds + useconds/1000000.0);
}

int main(int argc, char *argv[])
{
  int print = 0;
  int test;
  int errors = 0;
  unsigned long n, D, j;
  short *s; //manu
  int nthreads;
  timeval tstart, tend;

//./sievesimp 5000000000000000000 10000000 nthreads


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

  s = (short *) calloc(2*D+1,sizeof(short));

  gettimeofday(&tstart,0);
  OrigSegSiev(s,n,D,8,4);
  gettimeofday(&tend,0);

  if (print)
    for(j=0; j<=2*D; j++)
			if (s[j])
	    	printf("%lu\n",n-D+j);

  long seconds = tend.tv_sec-tstart.tv_sec;
  long useconds = tend.tv_usec-tstart.tv_usec;

  fprintf(stderr,"Seconds: %lf\n",seconds + useconds/1000000.0);
  //fprintf(stderr,"tzero: %f tone: %f ttwoL %f \n",tzero, tone, ttwo);
}
