#include <stdlib.h>
#include <stdio.h>
#include <gmpxx.h>
#include <math.h>
#include <sys/time.h>
#include "diophappr.h"
#include "simplesieve.h"
#include <omp.h>
#include <algorithm>    // std::min
#include <vector>

/*  Compile using
g++ sievesimp.cpp -osievesimp -O3 -frounding-math -mfpmath=387 -finline-functions -lgmp -fopenmp */
double tzero = 0, tone = 0, ttwo = 0;
long counter = 0, counter1 = 0;
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

void DoSieving(short *s, long k, long q, long M, long m0, long Mp, long D, unsigned long n, long ainv, long cainv, long chunksize) {

    long j,m, r0;
    unsigned long np;
    
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
} // DoSieving end

long calcEnd (long D, long quot, long n, long K, mpz_class npD) {
  mpz_class sqt, M2kap, sqtR, kmpz;
  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  long x = sqt.get_si();                     /* x = (int) sqrt(n+D) */
  long R;
  long end = 0;
  for(long M = K*D+1; M <=x; M+=(2*R+1)) { //unabh: M, K, D, R
    M2kap = ((((((mpz_class) M)*M)*D)/quot)/n);  /* (M^2)*D/(quot*n) */
    mpz_sqrt(sqtR.get_mpz_t(), M2kap.get_mpz_t());
    R = sqtR.get_si();  /* (int) (M*sqrt(kappa*D/n)) */
    end++;
  }
  return end;
}


std::vector<long> calcStartVals (long D, long quot, long n, long K, mpz_class npD, long chunksize) {
  mpz_class sqt, M2kap, sqtR, kmpz;
  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  long x = sqt.get_si();                     /* x = (int) sqrt(n+D) */
  long R;
  std::vector<long> MStartVal;
  long end2 = 0;
  for(long M = K*D+1; M <=x; M+=(2*R+1)) { 
    if (end2++ % chunksize == 0) MStartVal.push_back(M);
    M2kap = ((((((mpz_class) M)*M)*D)/quot)/n);  // (M^2)*D/(quot*n) 
    mpz_sqrt(sqtR.get_mpz_t(), M2kap.get_mpz_t());
    R = sqtR.get_si();  // (int) (M*sqrt(kappa*D/n)) 
    //end++;
  }
  printf("MVals: %li and %li", MStartVal[0], MStartVal[1]);
  return MStartVal;
}

void NewSegSiev(short *s, unsigned long n, long D, long K, long quot, int nthreads)
/* ensure: s[j] = (n-D+j is prime) for 0<=j<=2 D */
/* parameter kappa = 1/quot */
{
  timeval t1, t2, t3, t4;
  mpz_class npDtmp(n+D), sqttmp;
  long x;

  mpz_sqrt(sqttmp.get_mpz_t(), npDtmp.get_mpz_t());
  x = sqttmp.get_si();                     /* x = (int) sqrt(n+D) */
  		  
  SubSegSiev(s,n-D,2*D,K*D);

  gettimeofday(&t1,0);


  long end = calcEnd(D, quot, n, K, npDtmp); 


  gettimeofday(&t2,0);
  long seconds = t2.tv_sec-t1.tv_sec;
  long useconds= t2.tv_usec-t1.tv_usec;
  fprintf(stderr,"Seconds 1: %lf\n",seconds + useconds/1000000.0);


  long chunksize = ceil((1.0*end)/nthreads);
  if (chunksize*nthreads < end) {printf("error"); return;}

  std::vector<long> MStartVal = calcStartVals(D, quot, n, K, npDtmp, chunksize);

  gettimeofday(&t3,0);
  seconds = t3.tv_sec-t2.tv_sec;
  useconds= t3.tv_usec-t2.tv_usec;
  fprintf(stderr,"Seconds 2: %lf\n",seconds + useconds/1000000.0);

  fprintf(stderr,"end: %li chunksize: %li MStartValLength: %li\n",end, chunksize, MStartVal.size());

  
  #pragma omp barrier

  #pragma omp parallel num_threads(nthreads)
   {
    mpz_class npD(n+D), sqt;
    mpz_class M2kap, sqtR, kmpz;
    long M, Mp, m, m0, R;
    rat alpha0, alpha1;
    mpq_class eta, etaq;
    long c,cainv,a,ainv,q,k;
    int tid = omp_get_thread_num();
    //printf("\ntid: %i\n", tid);
    long start = tid * chunksize;
    long stop = (tid+1) * chunksize; //> end ? end : (tid+1) * chunksize ;
    if (stop > end) stop = end;
    M = MStartVal[tid];
    //#pragma omp parallel for /*private(sqt, M2kap, sqtR, kmpz, R, M)*/ schedule(static,chunksize) num_threads(nthreads)
    for (long i = start; i < stop; i++/*, M += (2*R+1)*/) {
    //for (long i = 0; i < end; i++/*, M += (2*R+1)*/) {
    //for(M = K*D+1; M <=x; M+=(2*R+1)) { //unabh: M, K, D, R
      if (tid ==0)counter++; else counter1++;
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

      DoSieving(s, k, q, M, m0, Mp, D, n, ainv, cainv, chunksize);
      M += 2*R+1;   
    }
  }

  gettimeofday(&t4,0);
  seconds = t4.tv_sec-t3.tv_sec;
  useconds= t4.tv_usec-t3.tv_usec;
  fprintf(stderr,"Seconds 3: %lf\n",seconds + useconds/1000000.0);
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

  //fprintf(stderr, "Device count: %i/n",omp_get_num_devices());  

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
  NewSegSiev(s2,n2,D2,8,4, nthreads);
  gettimeofday(&tend,0);
  long seconds = tend.tv_sec-tstart.tv_sec;
  long useconds = tend.tv_usec-tstart.tv_usec;

  fprintf(stderr,"Seconds total: %lf\n",seconds + useconds/1000000.0);
  //fprintf(stderr,"tzero: %f tone: %f ttwoL %f \n",tzero, tone, ttwo);

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
  printf("counter: %li\n",counter+counter1);

}
