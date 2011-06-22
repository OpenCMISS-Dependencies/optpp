/*----------------------------------------------------------------------
/
/     These subroutines were taken from the book "NUMERICAL RECIPES 
/     IN C - The Art of Scientific Computing" by W.H. Press, 
/     B.P. Flannery, S.A. Teukolsky and W.T. Vetterling, published 
/     by Cambridge Univ. Press, New York, 1989.
/
/     All routines have been converted to double precision.
/
/
/----------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define DEBUG 0
#define EPS 1.0e-8
#define NR_END 1

void dfpmin(p,n,ftol,iter,itmax,fret,func,dfunc)
double p[],ftol,*fret,(*func)();
void (*dfunc)();
int n,*iter;
{
  int j,i,its;
  double fp,fae,fad,fac;
  double *xi,*g,*dg,*hdg,*dvector();
  double **hessin,**dmatrix();
  void dlinmin(),nrerror(),free_dmatrix(),free_dvector();
  
  hessin=dmatrix(1,n,1,n);
  xi=dvector(1,n);
  g=dvector(1,n);
  dg=dvector(1,n);
  hdg=dvector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,g);
  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) hessin[i][j]=0.0;
    hessin[i][i]=1.0;
    xi[i] = -g[i];
  }
  for (its=1;its<=itmax;its++) {
    *iter=its;
    /* linmin(p,xi,n,fret,func); */
    dlinmin(p,xi,n,fret,func,dfunc);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      free_dvector(hdg,1,n);
      free_dvector(dg,1,n);
      free_dvector(g,1,n);
      free_dvector(xi,1,n);
      free_dmatrix(hessin,1,n,1,n);
      return;
    }
    fp=(*fret);
    for (i=1;i<=n;i++) dg[i]=g[i];
    *fret=(*func)(p);
    (*dfunc)(p,g);
    for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
    for (i=1;i<=n;i++) {
      hdg[i]=0.0;
      for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=0.0;
    for (i=1;i<=n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
    }
    fac=1.0/fac;
    fad=1.0/fae;
    for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
    for (i=1;i<=n;i++)
      for (j=1;j<=n;j++)
	hessin[i][j] += fac*xi[i]*xi[j]
	  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
    for (i=1;i<=n;i++) {
      xi[i]=0.0;
      for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  free_dvector(hdg,1,n);
  free_dvector(dg,1,n);
  free_dvector(g,1,n);
  free_dvector(xi,1,n);
  free_dmatrix(hessin,1,n,1,n);
  return;
}

#undef EPS
/*----------------------------------------------------------*/
#define ITMAX 100
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent(ax,bx,cx,f,df,tol,xmin)
double ax,bx,cx,tol,*xmin;
double (*f)(),(*df)(); /* ANSI: double (*f)(double),(*df)(double); */
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
	void nrerror();

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	dw=dv=dx=(*df)(x);
	for (iter=1;iter<=ITMAX;iter++) {
	  
/*	  if(DEBUG) printf(" fret in dbrent: %f\n",fx); */

		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx)  d1=(w-x)*dx/(dx-dw);
			if (dv != dx)  d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*f)(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=(*df)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
/*	nrerror("Too many iterations in routine DBRENT"); */
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3
/*----------------------------------------------------------*/
extern int ncom;	/* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)();

double f1dim(x)
double x;
{
	int j;
	double f,*xt,*dvector();
	void free_dvector();

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_dvector(xt,1,ncom);
	return f;
}
/*----------------------------------------------------------*/
extern int ncom;	/* defined in DLINMIN */
extern double *pcom,*xicom,(*nrfunc)();
extern void (*nrdfun)();

double df1dim(x)
double x;
{
	int j;
	double df1=0.0;
	double *xt,*df,*dvector();
	void free_dvector();

	xt=dvector(1,ncom);
	df=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	(*nrdfun)(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	free_dvector(df,1,ncom);
	free_dvector(xt,1,ncom);
	return df1;
}
/*----------------------------------------------------------*/
/*#define TOL 2.0e-4*/
#define TOL 0.02

int ncom=0;	/* defining declarations */
double *pcom=0,*xicom=0,(*nrfunc)();
void (*nrdfun)();

void dlinmin(p,xi,n,fret,func,dfunc)
double p[],xi[],*fret,(*func)();
void (*dfunc)();
int n;
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	double dbrent(),f1dim(),df1dim(),*dvector();
	void mnbrak(),free_dvector();

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	nrdfun=dfunc;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);

	if(DEBUG) printf(" fret from dbrent: %f\n",*fret);

	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}

#undef TOL
/*----------------------------------------------------------*/
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(ax,bx,cx,fa,fb,fc,func)
double *ax,*bx,*cx,*fa,*fb,*fc;
double (*func)();	/* ANSI: double (*func)(double); */
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
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
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT
/*----------------------------------------------------------*/
void nrerror(error_text)
char error_text[];
{
  void exit();

  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
/*--------------------------------------------------------*/
double *dvector(int nl,int nh)
{
  double *v;
  
  v=(double *)malloc((size_t) (nh-nl+1+NR_END)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl+NR_END;
}
/*--------------------------------------------------------*/
double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  double **m;
  
  m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;
  
  m[nrl]=(double *) malloc((size_t) ((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) 
    m[i]=m[i-1]+ncol;
  
  return m;
}

/*--------------------------------------------------------*/
void free_dvector(double *v,int nl,int nh)
{
  free((char*) (v+nl-NR_END));
}
/*--------------------------------------------------------*/
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
  int i;

  free((char*) (m[nrl]+ncl-NR_END));
  free((char*) (m+nrl-NR_END));
}

