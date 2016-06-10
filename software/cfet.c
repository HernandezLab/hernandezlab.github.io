/* cFET: continuous Fisher's Exact Test */
/* Developed by Ryan D. Hernandez, 2014 */
/* Email Ryan with any questions:  ryan.hernandez@ucsf.edu */

/* This program can be compiled at the command-line using: */
/* gcc cfet.c -lm -o cfet -O3 */
/* For help menu, type ./cfet */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <float.h>

#define EPS FLT_EPSILON
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static long double F;  /* density for observed table */
static long double N;  /* sum of cell entries */
static long double r1; /* row 1 sum */
static long double c1; /* column 1 sum */
static long double *table; /* 2x2 table */
static int VERBOSE = 0;

#define SIGN(a,b) ((b) >= 0.0 ? fabsl(a) : -fabsl(a))

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s \n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  abort();
}

long double LnNchooseK(int n, int k)
{
  long double v=0;
  
  v = lgammal(n+1.0) - lgammal(k+1.0) - lgammal(n-k+1.0);
  return v;
}

long double HypergeometricPMF(int a, int r1, int r2, int c1) 
{
  return expl(LnNchooseK(r1, a)+LnNchooseK(r2, c1-a)-LnNchooseK(r1+r2, c1));
}

long double brent(long double ax, long double bx, long double cx, long double (*f)(long double), long double tol, long double *xmin)
/* Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
   (such that bx is between ax and cx, and f(bx) is less than both f(ax) and 
   f(cx)), this routine isolates the minimum to a fractional precision of about
   tol using Brent's method. The abscissa of the minimum is returned as xmin, and
   the minimum function value is returned as brent, the returned function value.*/
{
  int iter;
  long double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  long double e=0.0; /* this will be the distance moved on the step before last*/
  
  a=(ax < cx ? ax : cx);  /*a and b must be in ascending order, */
  b=(ax > cx ? ax : cx);  /*but input abscissas need not be. */
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++){
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabsl(x)+ZEPS);
    if (fabsl(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabsl(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabsl(q);
      etemp=e;
      e=d;
      if (fabsl(p) >= fabsl(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      /*The above conditions determine the acceptability of the parabolic fit.
	Here we take the golden section step into the larger of the two segments*/
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      } 
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabsl(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u); /* This is the one function evaluation per iteration.*/

    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    }
    else { 
      if (u < x) a=u; 
      else b=u; 
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) { 
	v=u;
	fv=fu;
      }
    }
  }
  nrerror("Too many iterations in brent");
  *xmin=x;
  return fx;
}

long double zbrent(long double (*func)(long double), long double x1, long double x2, long double tol)
/* Using Brent's method, find the root of a function func known to
   lie between x1 and x2. The root, returned as zbrent, will be
   refined until its accuracy is tol.*/
{
  int iter;
  long double a=x1,b=x2,c=x2,d=.0,e=.0,min1,min2;
  long double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb >0.0) ||(fa <0.0 && fb < 0.0))
    nrerror("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc> 0.0) || (fb < 0.0 && fc <0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabsl(fc) < fabsl(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabsl(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabsl(xm) <= tol1 || fb== 0.0) return b; 
    if (fabsl(e) >= tol1 && fabsl(fa) >fabsl(fb)) {
      s=fb/fa;
      if (a == c){
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if(p > 0.0) q =-q;
      p=fabsl(p);
      min1=3.0*xm*q-fabsl(tol1*q);
      min2=fabsl(e*q);
      if(2.0*p < (min1 <min2 ?min1 :min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabsl(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=(*func)(b);
  }
  nrerror("Maximum number of iterations exceeded in zbrent"); 
  return 0.0;
}

void polint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy) 
/*Given arrays xa[1..n] and ya[1..n], and given a value x, this routine 
  returns a value y, and an error estimate dy. If P(x) is the polynomial
  of degree N-1 such that P(xai) = yai, i= 1,...,n, then the returned 
  value y=P(x). */
{
  int i,m,ns=1; 
  long double den,dif,dift,ho,hp,w;
  long double *c,*d; 
  
  dif=fabsl(x-xa[1]);
  assert(c = malloc((n+1)*sizeof(*c)));
  assert(d = malloc((n+1)*sizeof(*d)));
  for(i=1; i<=n; i++) { 
    if((dift=fabsl(x-xa[i])) < dif){ 
      ns=i; 
      dif=dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns--];
  for(m=1; m<n; m++){
    for(i=1; i<=n-m; i++){
      ho = xa[i]-x;
      hp = xa[i+m]-x;
      w = c[i+1]-d[i];
      if((den=ho-hp) == 0.0) nrerror("Error in routine polint...");
      den = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(d);
  free(c);
}

long double trapzd(long double (*func)(long double), long double a, long double b, int n)
/* This routine computes the nth stage of refinement of an extended 
   trapezoidal rule. func is input as a pointer to the function to be
   integrated between limits a and b, also input. When called with 
   n=1, the routine returns the crudest estimate of int_a^b{f(x)}dx.
   Subsequent calls with n=2,3,... (in that sequential order) will 
   improve the accuracy by adding 2n-2 additional interior points.*/
{
  long double x, tnm, sum, del;
  static long double s;
  int it, j;

  if(n==1) return(s=0.5*(b-a)*(func(a)+func(b)));
  else{
    for(it=1,j=1; j<n-1; j++) it <<= 1;
    tnm = it;
    del = (b-a)/tnm;
    x = a+0.5*del;
    for(sum=0.0,j=1; j<=it; j++,x+=del) sum += func(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return(s);
  }
}

long double qromb(long double (*func)(long double), long double a, long double b)
/* Returns the integral of the function func from a to b.  Integration is
   performed by Romberg's method of order 2K */
{
  void polint(long double xa[], long double ya[], int n, long double x, long double *y, long double *dy);
  long double trapzd(long double (*func)(long double), long double a, long double b, int n);
  void nrerror(char error_text[]);
  long double ss, dss;
  long double s[JMAXP], h[JMAXP+1]; /* store successive trap. approx & stepsizes */
  int j;

  h[1] = 1.0;
  for(j = 1; j <= JMAX; j++){
    s[j] = trapzd(func, a, b, j);
    if(j >= K){
      polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
      if(fabsl(dss) <= EPS*fabsl(ss)) return ss;
    }
    h[j+1] = 0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0; /* NEVER GET HERE... */
}

long double gammln(long double xx) 
/* Returns the value ln[gamma(xx)] for xx>0. */
{ 
  /* Internal arithmetic will be done in long double precision, a nicety that
     you can omit if five-figure accuracy is good enough. */
  long double x,y,tmp,ser; 
  static long double cof[6]={76.18009172947146,-86.50532032941677, 
			24.01409824083091,-1.231739572450155, 
			0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j; 
  y=x=xx; 
  tmp=x+5.5; 
  tmp -= (x+0.5)*logl(tmp); 
  ser=1.000000000190015; 
  for (j=0;j<=5;j++) ser += cof[j]/++y; 
  return -tmp+logl(2.5066282746310005*ser/x); 
}

long double LNgammaBico(long double z, long double w) 
/* Returns generalized binomial coefficients "z choose w" 
   using gamma function */
{ 
  return gammln(z+1)-gammln(w+1)-gammln(z-w+1); 
}

long double genHyper(long double x)
{
  return(expl(LNgammaBico(r1,x)+LNgammaBico(N-r1, c1-x)-LNgammaBico(N,c1)));
}

long double genHyperDiff(long double x)
{ /* density at x minus density at observed table */
  return(genHyper(x) - F);
}

long double genHyperNeg(long double x)
{
  return(-expl(LNgammaBico(r1,x)+LNgammaBico(N-r1, c1-x)-LNgammaBico(N,c1)));
}

long double FisherExact(int a, int b, int c, int d)
{
  int i, r1, r2, c1, minA=0, maxA=0;
  long double p0, pvalue=0.0, pfinal = 0.0;

  if(VERBOSE){
    printf("running standard Fisher Exact Test!\n");
    fflush(stdout);
  }

  r1 = (a + b);
  r2 = (c + d);
  c1 = (a + c);
  
  minA = (c1-r2 > 0 ? c1-r2 : 0);
  maxA = (r1 > c1 ? c1 : r1);

  pfinal = p0 = HypergeometricPMF(a,r1,r2,c1);
  i = minA;
  if(VERBOSE)  printf("a = %d; minA=%d; maxA=%d; p0 = %Lf\n",a, minA,maxA,p0);
  while((pvalue = HypergeometricPMF(i,r1,r2,c1)) < p0+EPS && i < a){
    if(VERBOSE)  printf("i=%d -> pvalue=%.6Le\n",i,pvalue);
    pfinal += pvalue;
    i++;
  }
  i=maxA;
  while((pvalue = HypergeometricPMF(i,r1,r2,c1)) < p0+EPS && i > a){
    if(VERBOSE)  printf("i=%d -> pvalue=%.6Le\n",i,pvalue);
    pfinal += pvalue;
    i--;
  }
  return (pfinal);
}
 
long double semiContFisher(int x){
  int i, min, max;
  long double PROBx, pval, denom=.0, tmp;

  if(VERBOSE){
    printf("running semi-Continuous Fisher!!\n");
    fflush(stdout);
  }

  min = (int)(0 > c1-(N-r1) ? 0 : c1-(N-r1));
  max = (int)(r1 < c1 ? r1 : c1);
  
  PROBx = genHyper(x+0.0);/* proportional to density for observed data */
  
  /* determine p-value by summing over configurations */
  pval = 0.0;
  i = min;
  while(i <= max){
    tmp = genHyper(i+.0);
    denom += tmp;
    if(tmp <= PROBx){
      pval += tmp;
    }
    i++;
  }
  return(pval/denom);
}

long double ContFisher(long double x){
  long double tmin, min, tmax, max, pval=0.0, denom, s, peak, xopt;
  
  if(VERBOSE){
    printf("running Continuous Fisher!!\n");
    fflush(stdout);
  }
  tmin = min = (0 > c1-(N-r1) ? 0 : c1-(N-r1))+EPS;
  tmax = max = (r1 < c1 ? r1 : c1)-EPS;
  if(VERBOSE){
    printf("min: %.6Le; nmax: %.6Le\n",min, max);
  }
  pval = genHyper(tmin);
  while(!isnan(pval) && pval>4*FLT_EPSILON){
    min = tmin;
    tmin -= (max-min)/1000;
    pval = genHyper(tmin);
  }
  pval = genHyper(tmax);
  while(!isnan(pval) && pval>4*FLT_EPSILON){
    max = tmax;
    tmax += (max-min)/1000;
    pval = genHyper(tmax);
  }  
  if(VERBOSE){
    printf("min -> %.6Le\nmax -> %.6Le\n",min, max);
    for(pval=min; pval<=max; pval+=(max-min)/100){
      printf("%.6Le ",pval);
    }
    printf("\n");
    for(pval=min; pval<=max; pval+=(max-min)/100){
      printf("%.6Le ",genHyper(pval));
    }
    printf("\n");
  }
  if(fabsl(max-min) > EPS){
    F = genHyper(table[0]); /* global variable! */
    peak = brent(min, (max+min)/2, max, *genHyperNeg, EPS, &xopt);
    if(F > (fabsl(peak+EPS))){
      return(1.0); /* optimal table observed! */
    }
    if(VERBOSE)
      printf("min = %Lf;  max = %Lf; r1 = %Lf; c1 = %Lf; N = %Lf;  F=%.10Le(%f:%d); peak(%Lf)=%.10Lf\n", min, max, r1, c1, N, F,EPS,F<EPS,xopt, peak);
    
    denom = qromb(*genHyper, min, max);
    if(VERBOSE)
      printf("denom = %Lf; genHyper(%Le+%Le)=%Le; genHyper(%Le-%Le)=%Le\n",denom,table[0],.2*table[0],genHyper(table[0]+.2*table[0]),table[0],.2*table[0],genHyper(table[0]-.2*table[0]));
    pval = 0.0;
    if(xopt > table[0]){
      /* function increases to the right */
      if(fabsl(table[0] - min) > EPS){
	/* integrate left tail */
	pval = qromb(*genHyper, min, table[0])/denom;
	if(VERBOSE)
	  printf("integrate left tail: (%Le,%Le) pval=%Lf\n",min, table[0], pval);
      }
      if(genHyper(max) < F){
	/*use brent's method to find s s.t. f(s) = f(table[0])*/
	/* then integrate right side */
	s = zbrent(*genHyperDiff, xopt, max, EPS);	  
	pval += qromb(*genHyper, s, max)/denom;
	if(VERBOSE)
	  printf("use brents method left: s=%Le->%Le; pval=%Lf\n",s,max,pval);
      }
    }
    else if(xopt < table[0]){
      /* function increases to the left */
      if(fabsl(table[0] - max) > EPS){
	/* integrate right tail */
	pval = qromb(*genHyper, table[0], max)/denom;
	if(VERBOSE)
	  printf("integrate right tail:(%Le,%Le)=pval=%Lf\n", table[0],max,pval);
      }
      if(genHyper(min) < F){
	/*use brent's method to find s s.t. f(s) = f(table[0])*/
	/* then integrate left side */
	s = zbrent(*genHyperDiff, min, xopt, EPS);	  
	pval += qromb(*genHyper, min, s)/denom;
	if(VERBOSE)
	  printf("use brents method right: s=%Le; pval=%Lf\n",s,pval);
      }
    }
    else{
      printf("not increasing or decreasing!!!\n");
      exit(1);
    }
  }
  else{
    pval = 1;
  }
  return(pval);
}

int main(int argc, char *argv[])
{
  int i, SIMPLE=0, MKPRF=0, arg, TEST=-1;
  int isINT[4];
  char args[100], name[100]="";
  long double tmp, pval=0.0;
  FILE *outfile=stdout, *infile=NULL;
  
  /* initialize global variable table */
  assert(table = malloc(4*sizeof(*table)));
  
  if(argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0){
    printf("program:     CFET (Continuous Fisher's Exact Test\n");
    printf("developer:   Ryan D. Hernandez (ryan.hernandez@ucsf.edu)\n");
    printf("date:        21 February, 2012\n");
    printf("Description: This program implements an analog of Fisher's\n");
    printf("             exact test, adapted to allow for continuous cell\n");
    printf("             entries.  This is done by exchanging factorial\n");
    printf("             functions of the hypergeometric distribution by\n");
    printf("             gamma functions [recall that gamma(n+1) = n! if n\n");
    printf("             is an integer].\n");
    printf("USAGE:       ./cfet [-options]\n\n");
    printf("The most direct way to use this program is to type:\n\n");
    printf("\t./CFET <x11> <x12> <x21> <x22>\n\n");
    printf("where the first index denotes the row, and the second the column\n");
    printf("such that the odds ration = x11*x22/x12/x21.\n");
    
    printf("INPUT DATA FILE (./cfet -i <data.txt> <-s|-m>)\n");
    printf("\tAlternatively, data for one or more tables may be put in a\n");
    printf("\tfile <data.txt>.  There are two file formats accepted:\n\n");
    printf("SIMPLE DATA FILE FORMAT (-s)\n");
    printf("\tFile contains 5 columns, with one table per row:\n");
    printf("\t\tName <x11> <x12> <x21> <x22>\n\tNote: Name can be index\n\n");
    printf("MKPRF DATA FILE INPUT (-m)\n");
    printf("\tInput data as used by mkprf (Bustamante, et. al, Nature 2005)\n");
    printf("\tThis format has one gene per row, and 11 columns for each gene:");
    printf("\n\t1.  Gene class\n\t2.  gene name\n\t3.  FS\n\t4.  PS\n");
    printf("\t5.  FN\n\t6.  PN\n");
    printf("\t7.  n1 (number of sequences sampled in species 1)\n");
    printf("\t8.  n2 (number of sequences sampled in species 2)\n");
    printf("\t9.  TS (total number of silent sites in alignment)\n");
    printf("\t10. TR (total number of replacement sites in alignment\n");
    printf("\t11. Ratio (haploid ratio: 1 for autosomal, .75 for X, .25 for Y");
    printf(")\n\tNOTE:  only columns 2-6 will be used!\n\n");
    printf("OUTPUT (-o <file>)\n");
    printf("\tOutput is printed in the following form: [geneName] p-value\n");
    printf("\tBy default, output is to the screen.  Use -o <file> to print\n");
    printf("\toutput to a file.\n\n");
    printf("VERBOSE OUTPUT  (-v)\n");
    printf("\tTo print out EXCESSIVE details for each table, add -v flag\n");
    printf("SPECIFY TEST (-t <D/S/C>)\n");
    printf("\tTo specify type of test to run (<D>iscrete, <S>emi-continuous,\n\
\tor <C>ontinuous.\n");
    exit(1);
  }
  /* parse argument string */
  arg = 1;
  while(arg < argc){
    if(argv[arg][0] == '-'){
      strcpy(args, &argv[arg][1]);
      if(strcmp(args, "i") == 0){
	arg++;
	assert(infile = fopen(argv[arg],"r"));
	arg++;
      }
      else if(strcmp(args, "m") == 0){
	MKPRF = 1;
	if(SIMPLE == 1){
	  printf("error:  must specify either -s or -m (NOT BOTH)\n");
	  exit(1);
	}
	arg++;
      }
      else if(strcmp(args, "o") == 0){
	arg++;
	assert(outfile = fopen(argv[arg],"w"));
	arg++;
      }
      else if(strcmp(args, "s") == 0){
	SIMPLE = 1;
	if(MKPRF == 1){
	  printf("error:  must specify either -s or -m (NOT BOTH)\n");
	  exit(1);
	}
	arg++;
      }
      else if(strcmp(args, "t") == 0){
	arg++;
	switch(argv[arg][0]){
	case 'D': /* run discrete test */
	  TEST=0;
	  break;
	case 'S': /* run semi-continuous test */
	  TEST=1;
	  break;
	case 'C': /* run continuous test */
	  TEST=2;
	  break;
	default:
	  fprintf(stderr,"CFET error: improper use of option \"-t\" to specify \
test.  please use either <D> for discrete test, <C> for continuous test, or <S>	\
 for semi-continuous test.\n\n");
	  exit(1);
	  break;
	}
	arg++;
      }
      else if(strcmp(args, "v") == 0){
	VERBOSE = 1;
	arg++;
      }
      else{
	printf("Unrecognized option -%s.  For help type ./CFET\n",args);
	exit(1);
      }
    }
    else{
      table[0] = atof(argv[arg++]);
      table[1] = atof(argv[arg++]);
      table[2] = atof(argv[arg++]);
      table[3] = atof(argv[arg++]);
    }
  }
  
  if(infile != NULL && SIMPLE == 0 && MKPRF == 0){
    printf("Must identify file format!  Use -s (simple) or -m (mkprf)\n");
    exit(1);
  }
  
  int numTABLES = 0;
  while(1){ /* loop over possibly unknown number of MK tables */
    name[0] = '\0';
    numTABLES++;
    if(infile != NULL){ /* get MK table from file if necessary */
      if(SIMPLE == 1){
	if((i=fscanf(infile,"%s %Lf %Lf %Lf %Lf\n",&name[0], &table[0], 
		     &table[1], &table[2], &table[3])) == -1) /* EOF */
	  break;
	else if(i != 5){
	  printf("error in input file: only %d elements read\n",i);
	  abort();
	}
      }
      else{ /* read MKPRF file */
	if((i=fscanf(infile,"%s %Lf %Lf %Lf %Lf %*f %*f %*f %*f %*f\n", name,
		     &table[2], &table[0], &table[3], &table[1])) == -1){
	  break;
	}
	else if(i != 5){
	  printf("error reading input file, only read %d elements from %s\n",
		 i,name);
	  abort();
	}
      }
      if(VERBOSE){
	printf("\nGENE:  %s\n",name);
      }
      if(name[0] != '\0')
	strcat(name," "); /* add space after gene name */
    }
    if(table[0] < 0 || table[1] < 0 || table[2] < 0 || table[3] < 0){
      printf("error:  cell entries must be >=0 (read %.1Lf %.1Lf %.1Lf %.1Lf)\n",
	     table[0], table[1], table[2], table[3]);
      exit(1);
    }
    
    r1 = table[0]+table[1]; /* top row */
    c1 = table[0]+table[2]; /* first column */
    N = 0.0;
    for(i=0; i<4; i++)
      N += table[i];
    
    if(VERBOSE){
      printf("data:\n%Lf %Lf | %Lf\n%Lf %Lf | %Lf\n---------\
\n%Lf %Lf\n", table[0], table[1], r1, table[2], table[3], N-r1,
	     c1, N-c1);
    }
    /*
      TO DO:  have 3 implementations of Fisher's exact test:
      (1) all 4 entries integer valued        => discrete Fisher's exact test
      (2) 1 row or column of discrete values  => generalized binomial theorem
      (3) otherwise, run general              => continuous Fisher's exact test
    */
    for(i=0; i<4; i++){
      isINT[i] = 0;
      if(table[i] - floor(table[i]) < EPS)  isINT[i] = 1;
    }
    
    if(isINT[0] + isINT[1] + isINT[2] + isINT[3] == 4 &&
       (TEST == 0 || TEST == -1)){
      /* integer valued data, run standard fisher's exact test */
      pval = FisherExact((int)(floor(table[0])), (int)(floor(table[1])),
			 (int)(floor(table[2])), (int)(floor(table[3])));
    }
    else if(isINT[0] + isINT[1] == 2 && (TEST == 1 || TEST == -1)){ 
      /* row 1 discrete */
      pval = semiContFisher((int)(table[0]));
    }
    else if(isINT[0] + isINT[2] == 2 && (TEST == 1 || TEST == -1)){
      /* column 1 discrete */
      /* transpose matrix */
      tmp = table[1];
      table[1] = table[2];
      table[2] = tmp;
      tmp = r1;
      r1 = c1;
      c1 = tmp;
      pval = semiContFisher((int)(table[0]));
    }
    else if(isINT[2] + isINT[3] == 2 && (TEST == 1 || TEST == -1)){
      /* row 2 discrete */
      /* swap rows */
      tmp = table[0];
      table[0] = table[2];
      table[2] = tmp;
      tmp = table[1];
      table[1] = table[3];
      table[3] = tmp;
      r1 = N-r1;
      pval = semiContFisher((int)(table[0]));
    }
    else if(isINT[1] + isINT[3] == 2 && (TEST == 1 || TEST == -1)){
      /* column 2 discrete */
      /* transpose matrix across anti-diagonal */
      tmp = table[0];
      table[0] = table[3];
      table[3] = tmp;
      r1 = table[0] + table[1];
      c1 = table[0] + table[2];
      pval = semiContFisher((int)(table[0]));
    }
    else if(TEST == 2 || TEST == -1){ /* rows and columns are continuous */
      pval = ContFisher(table[0]);
    }
    else{
      fprintf(stderr, "CFET error:  improper data, cannot run test.  Data read:\
\t%Lf  %Lf\n\t%Lf %Lf\nIf specific test specified (-t flag), check that the data \
has the right format (i.e. cannot run discrete test on continuous data).\n",
	      table[0], table[1], table[2], table[3]);
      exit(1);
    }
    if(name[0] != '\0')
      fprintf(outfile,"%s ",name);
    else
      fprintf(outfile,"%d ", numTABLES);
    for(i=0;i<4;i++){
      if(isINT[i] == 1)
	fprintf(outfile, "%d ", (int)(table[i]+EPS));
      else
	fprintf(outfile, "%Lf ", table[i]);
    }
    fprintf(outfile, "%.6Le\n", pval);
    if(infile == NULL)  break; /* only one MK table entered at command line */
  }
  free(table);
  if(infile != NULL)
    fclose(infile);
  fflush(stdout);
  return(0);
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K
#undef ITMAX
