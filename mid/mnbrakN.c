#include <math.h>
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


void NmnbrakN(int Nvec, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(int Nvec,double,int *locator,unsigned long Nspikes,unsigned long Ntrials), double *p, double *xi, unsigned long n,double *stimuli,int *locator,unsigned long Nspikes,unsigned long Ntrials)
{
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(Nvec,*ax,locator,Nspikes,Ntrials);
  *fb=(*func)(Nvec,*bx,locator,Nspikes,Ntrials);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
        SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(Nvec,*cx,locator,Nspikes,Ntrials);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*SIGN(DMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(Nvec,u,locator,Nspikes,Ntrials);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(Nvec,u,locator,Nspikes,Ntrials);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(Nvec,u,locator,Nspikes,Ntrials);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
            SHFT(*fb,*fc,fu,(*func)(Nvec,u,locator,Nspikes,Ntrials))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(Nvec,u,locator,Nspikes,Ntrials);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(Nvec,u,locator,Nspikes,Ntrials);
    }
    SHFT(*ax,*bx,*cx,u)
        SHFT(*fa,*fb,*fc,fu)
  }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
//#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #.3. */
