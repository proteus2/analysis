/* sidi.c						*/
/* Testprogramm fuer sidi()		*/
/* (C) 1992/2011 H. Weber       */
/* ANSI/ISO C                   */
/* Komplexe Version             */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

static complex double F(complex double s);
static double f(double t);
double acosh(double x);
double U(double x);
double *davec(long nl, long nh);
void dfvec(double *v, long nl, long nh);
double j0(double x);
static void wahl();
static double pi = M_PI;

int linvsidi(complex double (*F)(complex double s),
			 int p,
		     double c,
			 int mm,
			 double *t,
			 double *y);

int cinitausd(int pn, char *instring);
int causwert(int pn, complex double parval[], complex double *result);
void cparinfo(int pn, int wahl, int *anzahl, int i, char *name);
int cinitexpr();
void cfreeexpr();
void cinstrukt();
int cnewexpr(int p, char *str, char *para);
int cnewexprx(int p, char *str, int nparas, char *para[]);
int cnewexprs(int p, int nparas, char *para[], char *arg);
int cpermut(int expid, int parid);

void mkctime();

int beispiel;

int main()
{
	double *t, *y;
	double a0, b0, h, c, erg, wahr, fehl, maxfehl;
	int i, m, p = 1, rv;
	char s[2];
	clock_t time1, time2;

	printf("Inversion of the Laplace Transform for given values of the transformed function\n");
	printf("for Complex Arguments (Bromwich Integral, Algorithm of Sidi)\n");
	printf("Version with Formula Parser\n\n");

	do {
		wahl("Sidi");
		if (beispiel == 0) break;

		fflush(stdin);
		if (beispiel == -1)
		{	  
			cinitexpr();
			rv = cnewexpr(1, "Expression for f(s): " , "S");
			if (rv > 0) {
				printf("Error in Expression for f\n");
				continue;
			}         
		}

		printf("\nEnter p, c (> conv. absc.), a0, b0, m: ");
		scanf("%d %lf %lf %lf %d", &p, &c, &a0, &b0, &m);

		t = davec(1, m);
		y = davec(1, m);

		h = (b0 - a0)/(m - 1);
		for (i = 1; i <= m; i++)
		{
			double tx;
			tx = a0 + (i-1)*h;
			t[i] = tx;
		}

		time1 = clock();
		rv = linvsidi(F, p, c, m, t, y);
		time2 = clock();

		printf("\nError %d : ", rv);
		if (rv == 1)
			printf("m <= 0\n");
		else if (rv == 2)
			printf("there is a t[i] <= 0\n");
		else if (rv == 3)
			printf("p <= 0\n");
		if (rv != 0) goto l99;

		if (beispiel == -1) {
			printf("\n     t         Result         \n");
			for (i = 1; i <= m; i++) 
			{
				erg = y[i];
				printf(" %11.4e %23.15e \n", t[i], erg);
			}
			printf("\n");
			printf("Computing time : %lf sec\n", (double)(time2 - time1)/((double)CLOCKS_PER_SEC));
		} else {
			maxfehl = 0;
			printf("\n     t         Result                 True Solution           Error  \n");
			for (i = 1; i <= m; i++) 
			{
				wahr = f(t[i]);
				erg = y[i];
				fehl = fabs(erg - wahr);
				maxfehl = (maxfehl > fehl) ? maxfehl : fehl;
				printf(" %11.4e %23.15e %23.15e %11.4e \n", t[i], erg, wahr, fehl);
			}
			printf("\nMax. abs. Error: %11.4e\n", maxfehl);
			printf("Computing time : %lf sec\n", (double)(time2 - time1)/((double)CLOCKS_PER_SEC));
		}
		mkctime();

l99:
		fflush(stdin);
		gets(s);
        if (beispiel == -1)
			cfreeexpr();
		dfvec(t, 1, m);
		dfvec(y, 1, m);

	} while (1);

	printf("\nEnd ...\n");
	fflush(stdin);
	gets(s);
	gets(s);

	return 0;
}

void mkctime()
{
    time_t t;
    char *ss;

    time(&t);
	ss = ctime(&t);
	printf("%s\n", ss);
   	return;
}

/* --------------------------------------------------------------------

            DOKUMENTATION

Funktion:   linvsidi

Zweck:      linvsidi berechnet die Inverse f der Laplace-Transformierten

					  /00
			F(s) =   | exp(-s*t)f(t)dt, Re(s) > c
                    0/
            
			mit dem Algorithmus von Dubner&Abate und Sidi

Parameter
Eingabe:    F		  - Pointer auf Funktion, die die Transformierte 
						F(s) berechneten, fuer komplexes Argument s 
			p         - Genauigkeit der iterierten Gauss-Integration:
			            p mal Formel vom Grad 64 (p >= 1)
			c		  - Parameter c des Algorithmus, muss groesser als
			            der Realteil der groessten Singularitaet von
						F sein
			mm        - Anzahl der Punkte, fuer die f(t) berechnet 
				        werden soll (mm >= 1)
            t         - eindim. Array (1 ... m) mit den t-Werten, fuer
					    die f(t) zu berechnen ist. 

Ausgabe		y 		  -	eindim. Array (1 ... m) mit den Resultaten
					    der Berechnungen fuer die Werte t[i]

Resultat	- Fehlerindikator
					  1 : m <= 0
					  2 : fuer i ist ein t[i] <= 0

Literatur:  H. Dubner, J. Abate, Numerical Inversion of Laplace 
			Transforms and the Finte Fourier Transform, JACM 15(1968), 
			115-123
			A.M. Cohen, Numerical Methods for Laplace Transform 
			Inversion, Springer-Verlag, New York 2007

Bemerkung:	Version mit komplexer Arithmetik (C99)  
--------------------------------------------------------------------- */

int linvsidi(complex double (*F)(complex double sr),
			 int p,
			 double c,
			 int mm,
			 double *t,
			 double *y)
{
	void gauss(int n, double *x, double *w);
	void itgauss(int n, int p, double *x, double *w);
	void extrap(int m, double *psi, double *f, double *z, double *a);
    int ngauss = 64;
	int n;
	const int m = 40;
	double *x, *w;
	double *psia, *fa, *za;
	double *psib, *fb, *zb;
	double a, b, tt, xa, xb;
	int i, k, l, rv;

	if (p <= 0) return 3;
	if (mm <= 0) return 1;

	/* Modifikation, p-fach iterierte Quadraturformel */
	n = ngauss*p;  

	x = davec(0, n);
	w = davec(0, n);
	psia = davec(0, m);
	fa = davec(0, m);
	za = davec(0, m+1);
	psib = davec(0, m);
	fb = davec(0, m);
	zb = davec(0, m+1);

	rv = 0;

    gauss(ngauss, x, w);
	/* Modifikation, p-fach iterierte Quadraturformel */
	if (p > 0) itgauss(ngauss, p, x, w);

    for (l = 1; l <= mm; l++)
	{

		 tt = t[l];

		 if (tt <= 0) { 
			 rv = 2;
			 goto l99;
		 }

		 for (i=1; i<=m+1; i++) {
			 za[i-1] = (i*pi)/tt;
             zb[i-1] = za[i-1];
         }

         za[0] = za[0]/2;

         fa[0] = 0;
         k = 0;
         psia[k] = 0;

         for (i=1; i<=n; i++) 
		 {
			 double valr;
			 complex double val;

             xa = k*za[0] + za[0]*x[i];
			 val = F(c + I*xa);
			 valr = creal(val);
			 valr = valr*cos(tt*xa);
             psia[k] = psia[k] + w[i]*valr;
         }

         psia[k] = (pi*psia[k])/2;
         fb[0] = 0; 
         psib[k] = 0;

         for (i=1; i<=n; i++) 
		 {
			 double vali;
			 complex double val;
         
			 xb = k*zb[0] + zb[0]*x[i];
             val = F(c + I*xb);
			 vali = cimag(val);
			 vali = -vali*sin(tt*xb);
			 psib[k] = psib[k] + w[i]*vali;
         }

         psib[k] = pi*psib[k];
 
         for (k=1; k<=m; k++) 
		 {
			 psia[k] = 0;
             psib[k] = 0;

             for (i=1; i<=n; i++) 
			 {
				 double valr, vali;
				 complex double val;

                 xa = (2*k-1)*za[0] + 2*za[0]*x[i];
                 xb = k*zb[0] + zb[0]*x[i];

                 val = F(c + I*xa);
				 valr = creal(val);
			     valr = valr*cos(tt*xa);
                 psia[k] = psia[k] + w[i]*valr;

                 val = F(c + I*xb);
				 vali = cimag(val);
				 vali = -vali*sin(tt*xb);
                 psib[k] = psib[k] + w[i]*vali;
			 }

             psia[k] = psia[k]*pi;
             fa[k] = fa[k-1] + psia[k-1];
             psib[k] = psib[k]*pi;
             fb[k] = fb[k-1] + psib[k-1];
         }

         extrap(m, psia, fa, za, &a);
         a = 2*a/pi;
         a = (a*exp(c*tt))/tt;
         
		 extrap(m, psib, fb, zb, &b);
         b = 2*b/pi;
         b = (b*exp(c*tt))/tt;
         y[l] = (a+b)/2.0;

	}

l99:
	dfvec(x, 0, n);
	dfvec(w, 0, n);
	dfvec(psia, 0, m);
	dfvec(fa, 0, m);
	dfvec(za, 0, m+1);
	dfvec(psib, 0, m);
	dfvec(fb, 0, m);
	dfvec(zb, 0, m+1);
	
	return rv;
}

void extrap(int m, 
			double *psi, 
			double *f, 
			double *z, 
			double *a)
{
	double *p, *q;
	double b;
	int j, k;
	
	p = davec(0, m);
	q = davec(0, m);

	for (j=0; j<=m; j++) 
	{
		p[j] = f[j]/psi[j];
        q[j] = 1.0/psi[j];
    }
    
	for (k=1; k<=m; k++) 
	{
        for (j=0; j<=m-k; j++) 
		{
            p[j] = p[j+1] - p[j];
            q[j] = q[j+1] - q[j];         
			b = (1.0/z[j+k]) - (1.0/z[j]);
            p[j] = p[j]/b;
            q[j] = q[j]/b;
        }
          
		*a = p[0]/q[0];
	}

    return;
}

void gauss(int n, 
		   double *x, 
		   double *w)
{
	int i, j;

    x[1] = 0.999652520867886069728452812172818; 
    w[1] = 0.0008916403608482164736480395724898344;
    x[2] = 0.998170058385977639673462250338200;
    w[2] = 0.002073516630281233817643767864276530;
    x[3] = 0.995506685738372160369691191721652;
    w[3] = 0.003252228984489181428058680199989531;
    x[4] = 0.991668126942312958465649651078416;
    w[4] = 0.004423379913181973861515457329864311;
    x[5] = 0.986663413894955481870926753676136;
    w[5] = 0.005584069730065564409295246509604289;
    x[6] = 0.980504399826026859459307060948579;
    w[6] = 0.006731523948359321299030383342978218;
    x[7] = 0.973205687429201408031240745673632;
    w[7] = 0.007863015238012359660982997648769673;
    x[8] = 0.964784586065969787910745077279613;
    w[8] = 0.008975857887848671542522651000559224;
    x[9] = 0.955261068539251402878190334004165;
    w[9] = 0.01006741157676510468617015836427180;
    x[10] = 0.944657722997557052926702019136426;
    w[10] = 0.01113508690419162707964916519207727;
    x[11] = 0.932999699077046409880391692535079;
    w[11] = 0.01217635128435543666908877520453430;
    x[12] = 0.920314648126290181375845772347937;
    w[12] = 0.01318873485752732933584589631261280;
    x[13] = 0.906632657561398779870961669043152;
    w[13] = 0.01416983630712974161375565260011866;
    x[14] = 0.891986179471670703805110262606884;
    w[14] = 0.01511732853620123943398702990977421;
    x[15] = 0.876409953630265948305931887442847;
    w[15] = 0.01602896417742577679273375217394922;
    x[16] = 0.859940925085805413424470108915974;
    w[16] = 0.01690258091857080469578274105536263;
    x[17] = 0.842618156527116621281779185515688;
    w[17] = 0.01773610662844119190534657335762311;
    x[18] = 0.824482735627328669928880615996702;
    w[18] = 0.01852756427012002302020755090479165;
    x[19] = 0.805577677586196625124426485509274;
    w[19] = 0.01927507658930781456448124847340455;
    x[20] = 0.785947823101317017141939058329594;
    w[20] = 0.01997687056636017069332846306416800;
    x[21] = 0.765639732009947272829006951772228;
    w[21] = 0.02063128162131176430507814873681898;
    x[22] = 0.744701572853526478739263153510961;
    w[22] = 0.02123675756182679450366988395440869;
    x[23] = 0.723183008626732043992473857379458;
    w[23] = 0.02179186226466172668841393048686875;
    x[24] = 0.701135078981995801847883385630079;
    w[24] = 0.02229527908187828153006735501547241;
    x[25] = 0.678610079168834057975221307523101;
    w[25] = 0.02274581396370907223988549848563456;
    x[26] = 0.655661435995105478078756349280078;
    w[26] = 0.02314239829065720864797662461613062;
    x[27] = 0.632343581104383708186982086255010;
    w[27] = 0.02348409140810500866266314287729056;
    x[28] = 0.608711821870003542074824374494411;
    w[28] = 0.02377008285741515433114110347211160;
    x[29] = 0.584822210211996409018656814874135;
    w[29] = 0.02399969429822915386406308993567301;
    x[30] = 0.560731409648060277235188231746124;
    w[30] = 0.02417238111740147858488476357900890;
    x[31] = 0.536496560893899519724771470970169;
    w[31] = 0.02428773372075171346739953339198905;
    x[32] = 0.512175146331712216254477921426858;
    w[32] = 0.02434547850456986019168269536737497;
      
    j = n/2;
	for (i=1; i<=j; i++) {
        x[j+i] = 1.0 - x[j-i+1];
        w[j+i] = w[j-i+1];
    }
      
	return;
}

void itgauss(int n, 
			 int p,
			 double *x,
			 double *w)
{
	int k, i;
	double h;
	double dp = (double) p;

	/* Urspruengliche Abszissen u. Gewichte stauchen */
	for (i = 1; i <= n; i++) 
	{
		x[i] = x[i]/dp;
		w[i] = w[i]/dp;
	}

	h = 1.0/dp;

	for (k = 2; k <= p; k++)
	{
		int kk = k - 1;
		double hh = kk*h;	

		for (i = 1; i <= n; i++) 
		{
			x[kk*n + i] = hh + x[i];
			w[kk*n + i] = w[i];	
		}
	}

	return;
}


static void wahl(char *s)
{
	printf("Choice of Example (%s)\n", s);
   	printf("(-1) Formula parser\n");
	printf("(0)  End \n");
    printf("(1)  F(s) = 1/(s + 0.5)                f(t) = exp(-t/2) \n");
    printf("(2)  F(s) = (s+0.2)/((s+0.2)^2 + 1)    f(t) = exp(-0.2*t)*cos(t) \n");
    printf("(3)  F(s) = 1/s                        f(t) = 1 \n");
    printf("(4)  F(s) = 1/(s*s)                    f(t) = t \n");
    printf("(5)  F(s) = 1/(s*s*s)                  f(t) = t*t/2 \n");
    printf("(6)  F(s) = 1/sqrt(s)                  f(t) = 1.0/sqrt(pi*t) \n");
	printf("(7)  F(s) = 1/(s*s + 1)                f(t) = sin(t) \n");
	printf("(8)  F(s) = 1/sqrt(s*s + 1)            f(t) = J0(t) \n");
	printf("(9)  F(s) = exp(-1/s)/s                f(t) = J0(2*srqt(t)) \n");
	printf("(10) F(s) = 1/(s*sqrt(s))              f(t) = 2*sqrt(t/pi) \n");
	printf("(11) F(s) = log((s*s + 1)/(s*s + 4))   f(t) = 2*(cos(2*t) - cos(t))/t \n");
	printf("(12) F(s) = atan(1/s)                  f(t) = sin(t)/t \n");
	printf("(13) F(s) = exp(-s)/s                  f(t) = u(t-1) \n");
	printf("(14) F(s) = exp(-s)/(s*(s+1))          f(t) = 0 (t<=1), 1-exp(-(t-1)), (1<=t)\n");
	printf("(15) F(s) = 1/(s^2 + s + 1)            f(t) = 2/sqrt(3)*exp(-t/2)*sin(t*sqrt(3)/2)\n");
	printf("(16) F(s) = log(s)/s                   f(t) = -gamma - log(t)\n");
   	printf("(17) F(s) = 1/(s - 1)                  f(t) = exp(t)\n");
	printf("(18) F(s) = exp(-4*sqrt(s))            f(t) = 2*exp(-4/t)/sqrt(pi*t^3)\n");
	printf("(19) F(s) = 1/(s^3 - 8)                f(t) = 1/12*exp(-t)(exp(3*t) - cos(sqrt(3)*t)  \n");
	printf("                                              - sqrt(3)*sin(sqrt(3)*t)) \n");
	printf("(20) F(s) = 1/(s*(1+exp(s)))           f(t) = 0 (2k < t < 2k+1), 1 (2k+1 < t < 2k+2)\n");
	printf("(21) F(s) = 1/(s*(s+1))*(1+exp(-pi*s)) f(t) = |sin(t)|   \n");
	printf("            /(1-exp(-pi*s))                              \n");  
	printf("(22) F(s) = 1/((s-1)*sqrt(s))          f(t) = exp(t)*erf(sqrt(t))\n");
	printf("(23) F(s) = exp(-5*sqrt(s))/s          f(t) = 1-erf(5/(2*sqrt(t)))\n");

	do {
		char ss[10];
        printf("Choose : "); 
		gets(ss);
		beispiel = atoi(ss);
    } while (beispiel < -1 || beispiel > 23);
}

/*
double j0(double x)
{
	double ax, z, xx, y, ans, ans1, ans2;

	if ((ax=fabs(x)) < 8.0) {
	   y=x*x;
	   ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7
			+ y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))));
	   ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718
			+ y*(59272.64853 + y*(267.8532712 + y*1.0))));
	   ans = ans1/ans2;
	} else {
		z = 8.0/ax;
		y = z*z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4
			+ y*(-0.2073370639e-5 + y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y*(0.1430488765e-3
			+ y*(-0.6911147651e-5 + y*(0.7621095161e-6
			- y*0.934945152e-7)));
		ans = sqrt(0.636619772/ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
	}
	return ans;
}
*/


static complex double FF(complex double s)
{
	complex double parval[5];
	complex double result;
   
	parval[1] = s;
	causwert(1, parval, &result);
	return result;
}

static complex double F(complex double s)
{
    switch (beispiel) {
	   
	   case -1: return FF(s); 

       case 1 : return 1.0/(s + 0.5);

       case 2 : return (s + 0.2)/((s + 0.2)*(s + 0.2) + 1.0);

       case 3 : return 1.0/s;

       case 4 : return 1.0/(s*s);

       case 5 : return 1.0/(s*s*s);

       case 6 : return 1.0/csqrt(s);

       case 7 : return 1.0/(s*s + 1.0);

       case 8 : return 1.0/csqrt(s*s + 1.0);

       case 9 : return cexp(-1.0/s)/s;

       case 10: return 1.0/(s*csqrt(s));

       case 11: return clog((s*s + 1.0)/(s*s + 4.0));
       
	   case 12: return catan(1.0/s);

       case 13: return cexp(-s)/s;

	   case 14: return cexp(-s)/(s*(s + 1.0));

       case 15: return 1.0/(s*s + s + 1.0);

       case 16: return clog(s)/s;

       case 17: return 1.0/(s - 1.0);

	   case 18: return cexp(-4.0*csqrt(s));	

	   case 19: return 1.0/(s*s*s - 8.0);

	   case 20: return 1.0/(s*(1.0 + cexp(s)));

	   case 21: return 1.0/(s*s + 1.0)*(1.0 + cexp(-pi*s))/(1.0 - cexp(-pi*s));

	   case 22: return 1/((s-1.0)*csqrt(s));
	   
	   case 23: return cexp(-5.0*csqrt(s))/s;
	   
	   default: return 0 + I*0;
    }
}


static double f(double t)
{
	double gam = 0.577215664901532;
	double sqwave(double t);

    switch (beispiel) 
	{
        case 1 : return exp(-t*0.5);
        case 2 : return exp(-0.2*t)*cos(t);
        case 3 : return 1;
        case 4 : return t;
        case 5 : return t*t/2;
        case 6 : return 1.0/sqrt(pi*t);
	    case 7 : return sin(t);
	    case 8 : return j0(t);
	    case 9 : return j0(2*sqrt(t));
	    case 10: return 2*sqrt(t/pi);
	    case 11: return (t == 0) ? 0 : (2*(cos(2*t) - cos(t))/t);
	    case 12: return (t == 0) ? 1 : (sin(t)/t);
	    case 13: return U(t-1);
	    case 14: return (t <= 1) ? 0 : 1 - exp(-(t-1));
		case 15: return 2.0/sqrt(3)*exp(-t/2)*sin(t*sqrt(3)/2); 
		case 16: return -gam-log(t);
		case 17: return exp(t);
 		case 18: if (fabs(t) < 0.001) return 0.0; else return 2*exp(-4.0/t)/sqrt(pi*t*t*t);
 		case 19: return 1.0/12.0*exp(-t)*(exp(3.0*t) - cos( sqrt(3.0)*t) - sqrt(3.0)*sin(sqrt(3.0)*t) );
		case 20: return sqwave(t);
		case 21: return fabs(sin(t));
		case 22: return exp(t)*erf(sqrt(t));
		case 23: return 1.0 - erf(5.0/(2.0*sqrt(t)));

		default: return 0;
    }
}

double sqwave(double t)
{
	int r, tf;
	double tfl = floor(t);
	double eps = 1.0e-7;
	tf = (int) tfl;
	r = tf % 2;

	if (fabs(t - tfl) <= eps) 
		return 0.5;
	else if (r == 0) 
		return 0;	
	else
		return 1;
}


#define NR_END 1
#define FREE_ARG char*

static void aerror(char error_text[])
{
	fprintf(stderr,"NUMLIB Runtime error...\n");
	fprintf(stderr,"%s\n", error_text);
	exit(1);
}

/* --------------------------------------------------------------------*/
/* allocate a double vector with subscript range v[nl..nh] */

double *davec(long nl, long nh)
{
	double *v;

	v= (double *) malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) aerror("allocation failure in davec()");
	return v-nl+NR_END;
}

/* --------------------------------------------------------------------*/
/* free a double vector allocated with davec() */

void dfvec(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}


/*--------------------------------------------------------------------

            DOKUMENTATION

Funktion:   acosh
Zweck:      Berechnung des Area Cosinus Hyperbolicus

Parameter
Eingabe:    x   - Argument
---------------------------------------------------------------------*/
double acosh(double x)
{
    double sign;

    sign = 1.0;
    if (x < 0) {
      sign = -1.0;
      x = -x;
    }
    return log(x + sqrt(x*x - 1))*sign;
}


/*--------------------------------------------------------------------

            DOKUMENTATION

Funktion:   U
Zweck:      Berechnung der Einheits-Schrittfunktiom
Funktionswert: s.o.

---------------------------------------------------------------------*/
#define weps 1.0e-10

double U(double x)
{
	if (x < -weps) return(0);
	else if (x > weps) return(1);
	else return(0.5);
}
