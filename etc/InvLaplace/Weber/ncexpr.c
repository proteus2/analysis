/************************************************}
{                                                }
{   Numerical Work Bench                         }
{   ====================                         }
{                                                }
{   Helmut Weber Software, Hochheim              }
{   ANSI/ISO C                                   }
{   Modul ncexpr.c                               }
{												 }
{   Numerik-Programmbibliothek                   }
{   Version 5.01                                 }
{												 }
{   Komplexer Formelparser					     }
{	----------------------                       }
{												 }
{   Jan./Febr. 2012								 }
{   Copyright (C) 1988/2012                      }
{                                                }
{   Autor   : H. Weber                           }
{   						                     }
{================================================}

Das Modul NCExpr enthaelt die Funktionen

   cinitexpr
   cinitausd
   cfreeexpr
   causwert
   cinstrukt
   cparinfo
   cnewexpr
   cnewexprx
   cnewexprs
   cpermut

-------------------------------------------------*/

/*-----------------------------------------------------------------

            DOKUMENTATION

Funktion:   cinitexpr

Zweck:      Es werden die Datenstrukturen für die Auswertung von
            Ausdrücken angelegt.

Literatur:  R.L. Kruse, Data Structures & Program Design (2nd ed.)
            Prentice-Hall, Englewood Cliffs

Resultat:   0 : ok
            1 : zuwenig Speicherplatz

            
-------------------------------------------------------------------

            DOKUMENTATION

Funktion:   cinitausd

Zweck:      cinitausd initialisiert die Datenstrukturen fuer
            die Auswertung eines Ausdrucks. Es folgt das Parsing
            des Ausdrucks, gegebenfalls mit Fehlerausgabe, und
            die Erzeugung von Postfixcode zur Auswertung des
            Ausdrucks. Die eigentliche Auswertung geschieht mit
            Hilfe der Prozedur causwert. Zur Bildung der Ausdruecke
            koennen folgende Operatoren benutzt werden:

                         Prioritaet
               +             4
               -             4
               *             5
               /             5
               ^             6
               - (unaer)     6

            Ferner stehen folgende Standardfunktionen zur
            Verfuegung:

               abs           6
               sqr           6
               sqrt          6
               exp           6
               log           6
               lg (Basis 2)  6
               sin           6
               cos           6
               tan           6
               cot           6
               asin			 6
               acos			 6
               atan          6
               acot          6
               sinh          6
               cosh          6
               tanh          6
               coth          6

             Bekannte Konstanten sind

               e
               pi            
               pio2          pi/2
               pisq          pi*pi
               log10         log 10

Parameter
Eingabe:     pn       - Nummer der fuer den neuen Ausdruck zu
                        erzeugende Datenstruktur
             instring - String mit einzugebendem Ausdruck

Resultat:    fehler   - Fehlerindikator:
                        0 : o.k.
                        1 : Ausdruck falsch
						2 : pn falsch

Literatur:   R.L. Kruse, Data Structures & Program Design (2nd ed.)
             Prentice-Hall, Englewood Cliffs

---------------------------------------------------------------------

            DOKUMENTATION

Funktion:   cfreeexpr

Zweck:      cfreeexpr loescht die mit initexpr angelegten
            Datenstrukturen wieder

Literatur:  R.L. Kruse, Data Structures & Program Design (2nd ed.)
            Prentice-Hall, Englewood Cliffs

---------------------------------------------------------------------

            DOKUMENTATION

Funktion:   causwert

Zweck:      causwert wertet den mit cinitausd erzeugten Postfixcode
            fuer den eingegebenen Ausdruck aus. In der Regel
            werden dabei bestimmte Parameter benoetigt, z.B.
            fuer den Ausdruck x*x - 3*sin(x) wird der Wert des
            Parameters x benoetigt.

Parameter
Eingabe:    pn       - Nummer der vorher mit cinitausd
                       erzeugte Datenstruktur

            parval   - Parametervektor fuer die Auswertung.
                       Indexposition 1 enthaelt den ersten
                       Parameter, Position 2 den zweiten,
                       etc. (komplexe Groessen)

Ausgabe:    result   - Resultat der Auswertung
            fehler   - Fehlerindikator:
                       0 : o.k.
                       1 : Fehler bei der Auswertung

---------------------------------------------------------------------

            DOKUMENTATION

Prozedur:   cparinfo

Zweck:      cparinfo gibt Aufschluss ueber die Parameter des
            vorher eingegebenen Ausdrucks.

Parameter
Eingabe:    pn       - Nummer der vorher mit cinitausd
                       erzeugte Datenstruktur
            wahl     - 1 : es wird ueber die Anzahl der
                           Parameter informiert (anzahl)
                       2 : es wird der Name des Parameters
                           Nr. i zurueckgegeben (name)
            i        - wird nur im Fall wahl=2 gebraucht,
                       spezifiziert die Nummer des Parameters
                       (i >= 1, i <= anzahl)

Ausgabe:    anzahl   - Anzahl der verwendeten Parameter
            name     - Name des durch i spezifizierten Parameters

--------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <setjmp.h>
#include <complex.h>

#define M_PI    3.14159265358979323846
#define M_PI_2  1.57079632679489661923
#define M_LN10  2.30258509299404568402
                
#define True 1
#define False 0

#define firstunary 4
#define lastunary 25
#define firstbinary 26
#define lastbinary 32
#define firstoperand 33
#define lastoperand 38
#define maxstring 60
#define maxstack 100
#define maxtoken 50
#define maxpriority 6
#define maxparameter 10
#define namelength 10

#define hashsize 101
#define NEXPRESSIONS 10

typedef enum {
	operand, unaryop, binaryop, endexpression,
	leftparen, rightparen
} 
tokenkind;

typedef struct {
	char nm[namelength+1];
	tokenkind k;
	complex double Val;
	int pri;
} deftoken;

typedef struct {
	int postfix[maxstring+1];
	deftoken lexicon[maxtoken+1];
	int parameter[maxparameter+1];
	int nparameter;
} expinfo;

char errormsg[100];
extern int calc_error;

static expinfo *p, *parr[NEXPRESSIONS];
static jmp_buf env;

/*
double cot(double x);
double acot(double x);
double coth(double x);
double U(double x);
double round(double x);
*/

int cinitexpr()
{
	int i;

	for (i = 0; i < NEXPRESSIONS; i++) {
		parr[i] = (expinfo *) malloc(sizeof(expinfo));
		if (parr[i] == NULL) return 1;
	}

    return 0;
}


void cfreeexpr()
{
	int i;

	for (i = 0; i < NEXPRESSIONS; i++) 
       free(parr[i]);
}


static tokenkind kind(int t)
{
	return(p->lexicon[t].k);
}   /* kind */


static void error(char *text)
{
	strcpy(errormsg, text);
    strcat(errormsg, "\n\n");
    fprintf(stderr, errormsg);
	longjmp(env, 1);
}


static void definetokens()
{
	char *test;

	test = p->lexicon[1].nm;
	*test = '\0';
	p->lexicon[1].k = endexpression;

	strcpy(p->lexicon[2].nm, "(");
	p->lexicon[2].k = leftparen;

	strcpy(p->lexicon[3].nm, ")");
	p->lexicon[3].k = rightparen;

	strcpy(p->lexicon[4].nm, "~");      /* unary - */
	p->lexicon[4].k = unaryop;
	p->lexicon[4].pri = 6;

	strcpy(p->lexicon[5].nm, "ABS");
	p->lexicon[5].k = unaryop;
	p->lexicon[5].pri = 6;

	strcpy(p->lexicon[6].nm, "SQR");
	p->lexicon[6].k = unaryop;
	p->lexicon[6].pri = 6;

	strcpy(p->lexicon[7].nm, "SQRT");
	p->lexicon[7].k = unaryop;
	p->lexicon[7].pri = 6;

	strcpy(p->lexicon[8].nm, "EXP");
	p->lexicon[8].k = unaryop;
	p->lexicon[8].pri = 6;

	strcpy(p->lexicon[9].nm, "LOG");
	p->lexicon[9].k = unaryop;
	p->lexicon[9].pri = 6;

	strcpy(p->lexicon[10].nm, "LG");
	p->lexicon[10].k = unaryop;
	p->lexicon[10].pri = 6;

	strcpy(p->lexicon[11].nm, "SIN");
	p->lexicon[11].k = unaryop;
	p->lexicon[11].pri = 6;

	strcpy(p->lexicon[12].nm, "COS");
	p->lexicon[12].k = unaryop;
	p->lexicon[12].pri = 6;

	strcpy(p->lexicon[13].nm, "ATAN");
	p->lexicon[13].k = unaryop;
	p->lexicon[13].pri = 6;

#ifdef DEBUG
   fprintf(stderr, "definetokens(): nach atan\n");
#endif

	strcpy(p->lexicon[14].nm, "ROUND");
	p->lexicon[14].k = unaryop;
	p->lexicon[14].pri = 6;

	strcpy(p->lexicon[15].nm, "TRUNC");
	p->lexicon[15].k = unaryop;
	p->lexicon[15].pri = 6;

	strcpy(p->lexicon[16].nm, "TAN");
	p->lexicon[16].k = unaryop;
	p->lexicon[16].pri = 6;

	strcpy(p->lexicon[17].nm, "COT");
	p->lexicon[17].k = unaryop;
	p->lexicon[17].pri = 6;

	strcpy(p->lexicon[18].nm, "SINH");
	p->lexicon[18].k = unaryop;
	p->lexicon[18].pri = 6;

	strcpy(p->lexicon[19].nm, "COSH");
	p->lexicon[19].k = unaryop;
	p->lexicon[19].pri = 6;

	strcpy(p->lexicon[20].nm, "ACOT");
	p->lexicon[20].k = unaryop;
	p->lexicon[20].pri = 6;

	strcpy(p->lexicon[21].nm, "ASIN");
	p->lexicon[21].k = unaryop;
	p->lexicon[21].pri = 6;

	strcpy(p->lexicon[22].nm, "ACOS");
	p->lexicon[22].k = unaryop;
	p->lexicon[22].pri = 6;

	strcpy(p->lexicon[23].nm, "TANH");
	p->lexicon[23].k = unaryop;
	p->lexicon[23].pri = 6;

	strcpy(p->lexicon[24].nm, "COTH");
	p->lexicon[24].k = unaryop;
	p->lexicon[24].pri = 6;

	strcpy(p->lexicon[25].nm, "U");
	p->lexicon[25].k = unaryop;
	p->lexicon[25].pri = 6;

	/* binary operations */

	strcpy(p->lexicon[26].nm, "+");
	p->lexicon[26].k = binaryop;
	p->lexicon[26].pri = 4;

	strcpy(p->lexicon[27].nm, "-");
	p->lexicon[27].k = binaryop;
	p->lexicon[27].pri = 4;

	strcpy(p->lexicon[28].nm, "*");
	p->lexicon[28].k = binaryop;
	p->lexicon[28].pri = 5;

	strcpy(p->lexicon[29].nm, "/");
	p->lexicon[29].k = binaryop;
	p->lexicon[29].pri = 5;

	strcpy(p->lexicon[30].nm, "DIV");
	p->lexicon[30].k = binaryop;
	p->lexicon[30].pri = 5;

	strcpy(p->lexicon[31].nm, "MOD");
	p->lexicon[31].k = binaryop;
	p->lexicon[31].pri = 5;

	strcpy(p->lexicon[32].nm, "^");
	p->lexicon[32].k = binaryop;
	p->lexicon[32].pri = 6;

	/* operands */
	
	strcpy(p->lexicon[33].nm, "PI");
	p->lexicon[33].k = operand;
	p->lexicon[33].Val = M_PI;

	strcpy(p->lexicon[34].nm, "E");
	p->lexicon[34].k = operand;
	p->lexicon[34].Val = exp(1);

	strcpy(p->lexicon[35].nm, "PIO2");
	p->lexicon[35].k = operand;
	p->lexicon[35].Val = M_PI_2;

	strcpy(p->lexicon[36].nm, "PISQ");
	p->lexicon[36].k = operand;
	p->lexicon[36].Val = M_PI*M_PI;

	strcpy(p->lexicon[37].nm, "LOG10");
	p->lexicon[37].k = operand;
	p->lexicon[37].Val = M_LN10;

	strcpy(p->lexicon[38].nm, "I");
	p->lexicon[38].k = operand;
	p->lexicon[38].Val = 0 + I*1;
}                            /* definetokens */


/* Funktionen zu cinitausd */

static int hash(char *x, int h[])
{
	int a;
	char ch, *posit;
	int found;
	char s[255], *test, test2[255];

	test = (char *)strdup(x);
	if (test == NULL) {
		error("strdup failed -> not enough memory");
	}
	posit = (char *)strstr(x, " ");

	if (posit !=NULL) {
		strncpy(test, x, (posit - x));
		*(test + (posit - x)) = '\0';
	}

	ch = *test;
	a = ch % hashsize;

	do {
		if (h[a] == 0)
			found = True;
		else {
			strcpy(s, p->lexicon[h[a]].nm);
			posit = (char *)strstr(s, " ");
			if (posit == NULL)
				strcpy(test2, s);
			else {
				strncpy(test2, s, (posit - s));
				test2[posit - s] = '\0';
			}
			if (!strcmp(test2, test))
				found = True;
			else {
				if (strlen(test) == 1)
					ch = ' ';
				else
					ch = *(test+1);
				a += ch;
				if (a > hashsize) a = a % hashsize;
				found = False;
			}
		}
	} while (!found);

	free(test);
	return(a);
}                           /* hash */


static void makehashtable(int h[])
{
	int a;
	int t;

	for (a = 0; a <= hashsize; a++) h[a] = 0;

	for (t = 2; t <= lastoperand; t++) {
		h[hash(p->lexicon[t].nm, h)] = t;
	}
}                           /* makehashtable */


/* Unterfunktionen zu readexpression */

static void re_puttoken(int t, int *exprlength, int infix[])
{
	(*exprlength)++;
	if (*exprlength > maxstring) {
		error("Too many variables, operators and constants.");
	}
	infix[*exprlength] = t;
}                         /* re_puttoken */


static int leading(int *exprlength, int infix[])
{
	tokenkind k;
	int erg;

	if (*exprlength == 0)
		erg = True;
	else {
		k = kind(infix[*exprlength]);
		if ((k == leftparen) || (k == unaryop) || (k == binaryop))
			erg = True;
		else erg = False;
	}

	return(erg);
}                         /* leading */


static void extractword(char *instring, int *position, char *Word)
{
	int i;
	int start, stop;
	char ch[2];

	start = *position;
	while(isalnum(*(instring+(*position)))) (*position)++;

	stop = (*position)-1;
	if (stop >= start+namelength) {
/*
	  printf("Warning: The name ");
	  for (i = start; i <= stop; i++)
	    printf("%c", *(instring+i));
	  printf(" was shortened to %2d characters.\n", namelength);
*/
		stop = start+namelength-1;
	}
	*Word = '\0';  /* evtl. Fehler !!! */
	ch[1] = '\0';
	for (i = 1; i <= namelength; i++) {
		if (start+i-1 <= stop) {
			ch[0] = *(instring+start+i-1);
			ch[0] = toupper(ch[0]);
			strcat(Word, ch);
		}
		else strcat(Word, " ");
	}
}                         /* extractword */


static void findword(char *instring, int *position, int h[], int *exprlength, 
					 int infix[], int *tokencount)
{
	char Word[namelength+1];
	int t;

	extractword(instring, position, Word);
	t = h[hash(Word, h)];

	if (t != 0)
		if (leading(exprlength, infix))
			if (kind(t) == binaryop) {
				error("Binary operator at wrong position");
			}
			else {
				re_puttoken(t, exprlength, infix);
			}
		else
			if (kind(t) != binaryop) {
				error("Binary operator expected.");
			}
			else {
				re_puttoken(t, exprlength, infix);
			}
	else
		if (*tokencount >= maxtoken) {
			error("Too many different variables and constants.");
		}
		else if (!leading(exprlength, infix)) {
			error("Operand follows ) or other operand.");
		}
		else {
			(*tokencount)++;
			h[hash(Word, h)] = *tokencount;
			strcpy(p->lexicon[*tokencount].nm,  Word);
			p->lexicon[*tokencount].k = operand;
			if (p->nparameter >= maxparameter) {
				error("Expression contains too much parameters.");
			}
			else {
				(p->nparameter)++;
				p->parameter[p->nparameter] = *tokencount;
				re_puttoken(*tokencount, exprlength, infix);
			}
		}
}                         /* findword */


static double convertreal(char *x, int start, int stop)
{
	char ch;
	int i, j = 0;
	char zahl[60];

	for (i = start; i <= stop; i++) {
		ch = *(x+i);
		if (!isdigit(ch)) {
			error("Wrong character in number.");
		}
		zahl[j++] = ch;
	}
	zahl[j] = '\0';
	return(strtod(zahl, (char **)NULL));
}                         /* convertreal */


static void findnumber(int *exprlength, int infix[], int *tokencount, 
					   int *position, char *instring)
{
	int decpoint,
	newposition;
	double fraction, r;
	int i;

	if (!leading(exprlength, infix)) {
		error("constant in wrong position.");
	}
	else if (*tokencount >= maxtoken) {
		error("Too much constants and variables.");
	}
	else {
		newposition = *position;
		while (isdigit(*(instring+newposition)))
			newposition++;
		r = convertreal(instring, *position, newposition-1);

		if (*(instring+newposition) == '.') {
			decpoint = newposition;
			do {
				newposition++;
			} while (isdigit(*(instring+newposition)));
			fraction = convertreal(instring, decpoint+1, newposition-
			    1);
			for (i = 1; i < (newposition-decpoint); i++)
				fraction = fraction/10;
			r = r+fraction;
		}

		if (*(instring+newposition) == 'E'||
		    *(instring+newposition) == 'e') {
			error("Scientific notation not allowed here.");
		}
		else {
			(*tokencount)++;
			strcpy(p->lexicon[*tokencount].nm, "number");
			p->lexicon[*tokencount].k = operand;
			p->lexicon[*tokencount].Val = r;
			re_puttoken(*tokencount, exprlength, infix);
			*position = newposition;
		}
	}
}                         /* findnumber */

static void findsymbol(char *instring, int h[], int *exprlength, int infix[], 
					   int *position, int *parencount)
{
	char x[namelength+1];
	int t;

	x[0] = *(instring+(*position));
	x[1] = '\0';
	t = h[hash(x, h)];

	if (!t) {                      /* evtl. Fehler !!! */
		error("Unknown symbol in expression.");
	}
	else if (leading(exprlength, infix))
		if (kind(t) == rightparen) {
			error("Illegal place for closing parenthesis.");
		}
		else if (kind(t) == binaryop)
			if (!strcmp(x, "+")) t = t; /* dummy */
			else if (!strcmp(x, "-")) {
				strcpy(x, "~");           /* unary - */
				t = h[hash(x,h)];
				re_puttoken(t, exprlength, infix);
			}
			else {
				error("Binary operator in wrong position.");
			}
		else {
			re_puttoken(t, exprlength, infix);
		}
	else
		if (kind(t) == rightparen ||
		    kind(t) == binaryop) {
			re_puttoken(t, exprlength, infix);
		}
		else {
			error("Binary operator or ) expected.");
		}
	if (kind(t) == leftparen)
		(*parencount)++;
	else if (kind(t) == rightparen) {
		(*parencount)--;
		if (*parencount < 0) {
			error("More right than left parentheses.");
		}
	}
	(*position)++;
}                         /* findsymbol */


static void readexpression(char *instring, int *tokencount, int h[], 
						   int infix[])
{
	int lengthinstring;
	int exprlength;
	int position;
	int parencount;

	*tokencount = lastoperand;
	strcat(instring, " ");
	lengthinstring = strlen(instring);
	exprlength = 0;
	p->nparameter = 0;
	parencount = 0;
	position = 0;   /* war vorher 1 !!! */

	while (position < lengthinstring) { /* war vorher <= !!! */
		if (*(instring+position) == ' ')
			position++;
		else if (isalpha(*(instring+position)))
			findword(instring, &position, h, &exprlength, infix,
			    tokencount);
		else if (isdigit(*(instring+position)) ||
		    *(instring+position) == '.')
			findnumber(&exprlength, infix, tokencount, &position,
			    instring);
		else
			findsymbol(instring, h, &exprlength, infix, &position,
			    &parencount);
	}

	if (parencount) {
		error("Number of left parentheses not equal to number of right parentheses.");
	}
	if (leading(&exprlength, infix)) {
		error("Incomplete expression.");
	}
	re_puttoken(1, &exprlength, infix);
}                           /* readexpression */


static void tr_push(int t, int stack[], int *nstack)
{
	if (*nstack == maxstack) {
		error("Stack overflow.");
	}
	(*nstack)++;
	stack[*nstack] = t;
}                         /* tr_push */


static void tr_pop(int *t, int stack[], int *nstack)
{
	if (*nstack <= 0) {
		error("Stack empty.");
	}
	*t = stack[*nstack];
	(*nstack)--;
}                         /* tr_pop */


static void tr_gettoken(int *t, int *nget, int infix[])
{
	(*nget)++;
	if (*nget > maxstring) {
		error("Too many variables, constants and operators.");
	}
	*t = infix[*nget];
}                         /* tr_gettoken */


static void tr_puttoken(int t, int *nput)
{
	(*nput)++;
	if (*nput > maxstring) {
		error("Too many variables, constants and operators.");
	}
	p->postfix[*nput] = t;
}                         /* tr_puttoken */


static int tr_priority(int t)
{
	return(p->lexicon[t].pri);
}                         /* tr_priority */


static void translate(int infix[])
{
	int stack[maxstack+1];
	int nstack;
	int t, x;
	int endright;
	int nget, nput;

	nget = 0;
	nput = 0;
	nstack = 0;
	do {
		tr_gettoken(&t, &nget, infix);

		switch (kind(t)) {
		case  operand : 
			tr_puttoken(t, &nput); 
			break;
		case  leftparen : 
			tr_push(t, stack, &nstack); 
			break;
		case  rightparen : 
			tr_pop(&t, stack, &nstack);
			while (kind(t) != leftparen) {
				tr_puttoken(t, &nput);
				tr_pop(&t, stack, &nstack);
			}
			break;
		case unaryop:
		case binaryop : 
			do {
				if (!nstack) endright = True;
				else if (kind(stack[nstack]) == leftparen)
					endright = True;
				else if (tr_priority(stack[nstack]) <
					 tr_priority(t))
					endright = True;
				else if (tr_priority(stack[nstack])==
					 tr_priority(t) && tr_priority(t) == 6)
					endright = True;
				else {
					endright = False;
					tr_pop(&x, stack, &nstack);
					tr_puttoken(x, &nput);
				}
			} while (!endright);
			tr_push(t, stack, &nstack);
			break;
		case endexpression : 
			while (nstack > 0) {
				tr_pop(&x, stack, &nstack);
				tr_puttoken(x, &nput);
			}
			break;
		}
	} while (kind(t) != endexpression);
	tr_puttoken(t, &nput);
}                           /* translate */


int cinitausd(int pn, char *instring)
{
	int infix[maxstring+1];
	int h[hashsize+1];
	int tokencount;
	int xjump;

    if (pn < 0 || pn > NEXPRESSIONS-1) 
        return 2;

	p = parr[pn];
	definetokens();

	xjump = setjmp(env);
	/* Bei einem Fehler in einer der folgenden Subroutinen wird */
	/* an diese Stelle zurueckgesprungen und die Funktion beendet */
	if (xjump == 1) return(True);

	makehashtable(h);
	readexpression(instring, &tokencount, h, infix);
	translate(infix);

	return False;
}                            /* cinitausd */


/* Funktionen zu causwert */

static void aw_error(char *text)
{
	char buffer[150];

	sprintf(errormsg, "%s\n", text);
	strcpy(buffer, errormsg);
	fprintf(stderr, buffer);

	longjmp(env, 1);
    return;
}


static complex double getvalue(int t)
{
	if (kind(t) != operand) {
		aw_error("Non operand.");
	}
	else
		return p->lexicon[t].Val;
	return 0;
}                           /* getvalue */


static complex double dounary(int t, complex double x)
{
	complex double ergebnis = 0;

	if (t < firstunary || t > lastunary) {
		aw_error("Code for unary operator wrong.");
	}
	else
		switch (t) {
		case   4: 
			ergebnis = -x; 
			break;
		case   5: 
			ergebnis = cabs(x); 
			break;
		case   6: 
			ergebnis = x*x; 
			break;
		case   7: 
			ergebnis = csqrt(x); 
			break;
		case   8: 
			ergebnis = cexp(x); 
			break;
		case   9: 
			ergebnis = clog(x); 
			break;
		case  10: 
			ergebnis = clog(x)/log(2); 
			break;
		case  11: 
			ergebnis = csin(x); 
			break;
		case  12: 
			ergebnis = ccos(x); 
			break;
		case  13: 
			ergebnis = catan(x); 
			break;
		case  14: 
			ergebnis = round(x); 
			break;
		case  15: 
			ergebnis = floor(x); 
			break;
		case  16: 
			ergebnis = ctan(x); 
			break;
		case  17: 
			ergebnis = ccos(x)/csin(x); 
			break;
		case  18: 
			ergebnis = csinh(x); 
			break;
		case  19: 
			ergebnis = ccosh(x); 
			break;
		case  20: 
			ergebnis = M_PI_2 - catan(x); 
			break;
		case  21: 
			ergebnis = casin(x); 
			break;
		case  22: 
			ergebnis = cacos(x); 
			break;
		case  23: 
			ergebnis = ctanh(x); 
			break;
/*
		case  24: 
			ergebnis = ccoth(x); 
			break;
		case  25: 
			ergebnis = U(x); 
			break;
*/
		}

	return ergebnis;
}                           /* dounary */


/* Unterfunktion zu dobinary */
/*
static double exponent(complex double x, complex double y)
{
	double epsilon = 0.000001;
	int i;
	double p;
	double erg = 0;

	if (fabs(y-round(y)) < epsilon) {
		p = 1.0;
		if (y >= 0.0)
			for (i = 1; i <= round(y); i++) p = p*x;
		else if (x == 0.0) {
			aw_error("Negative power of 0.0");
		}
		else for (i = -1; i >= round(y); i--) p = p/x;
		erg = p;
	}
	else if (x > 0.0)
		erg = exp(y*log(x));
	else if (fabs(x) < epsilon)
		erg = 0.0;
	else
		aw_error("Attempt to compute non-integral power of negative number.");
	return(erg);
}                         // exponent
*/

static complex double dobinary(int t, complex double x, complex double y)
{
	complex double ergebnis = 0;

	if (t < firstbinary || t > lastbinary) {
		aw_error("Binary operator code wrong.");
	}
	else
		switch (t) {
		case 26: 
			ergebnis = x + y; 
			break;
		case 27: 
			ergebnis = x - y; 
			break;
		case 28: 
			ergebnis = x * y; 
			break;
		case 29: 
			if (y !=0) ergebnis = x / y;
			else aw_error("Division by zero");
			break;
/*
		case 30: 
			if (round(y) !=0)
				ergebnis = floor(round(x) / round(y));
			else aw_error("Division by zero");
			break;
		case 31: 
			if (round(y) != 0)
				ergebnis = fmod(round(x), round(y));
			else aw_error("Division by zero");
			break;
*/
		case 32: 
			ergebnis = cpow(x, y); 
			break;
		}
		
	return ergebnis;
}                          /* dobinary */


static void aw_push(complex double v, complex double stack[], int *nstack)
{
	if (*nstack == maxstack) {
		aw_error("Stack overflow.");
	}
	(*nstack)++;
	stack[*nstack] = v;
}                          /* aw_push */


static void aw_pop(complex double *v, complex double stack[], int *nstack)
{
	if (*nstack <= 0) {
		aw_error("Stack empty.");
	}
	*v = stack[*nstack];
	(*nstack)--;
}                             /* aw_pop */


static void aw_gettoken(int *t, int *nget)
{
	(*nget)++;
	if (*nget > maxstring) {
		aw_error("Too many variables, constants and operators.");
	}
	*t = p->postfix[*nget];
}                             /* aw_gettoken */


static void aw_setparameter(int i, complex double v)
{
	p->lexicon[p->parameter[i]].Val = v;
}                            /* aw_setparameter */


int causwert(int pn, complex double parval[], complex double *result)
{
//	double stack[maxstack+1];
	complex double stack[maxstack+1];
	int nstack;
	int t;
//	double x, y;
	complex double x, y;
	int i, nget;
	int xjump;

	p = parr[pn];

	for (i = 1; i <= p->nparameter; i++)
		aw_setparameter(i, parval[i]);

	calc_error = 0;

	xjump = setjmp(env);
	/* Bei einem Fehler in einer der folgenden Subroutinen wird */
	/* an diese Stelle zurueckgesprungen und die Funktion beendet */
	if (xjump == 1) return(True);

	nget = 0;
	nstack = 0;
	do {
		aw_gettoken(&t, &nget);
		switch (kind(t)) {
		case operand :
			aw_push(getvalue(t), stack, &nstack); 
			break;
		case unaryop : 
			aw_pop(&x, stack, &nstack);
			aw_push(dounary(t, x), stack, &nstack);
			break;
		case binaryop : 
			aw_pop(&y, stack, &nstack);
			aw_pop(&x, stack, &nstack);
			aw_push(dobinary(t, x, y), stack, &nstack);
			break;
		case endexpression :
			if (nstack == 1) 
                aw_pop(result, stack, &nstack);
			else
				aw_error("Error in execution.");
		case leftparen : ;
		case rightparen : ;
        }
		if (calc_error) break;
	} while (kind(t) != endexpression);

	return(calc_error);
}                             /* causwert */


void cparinfo(int pn, int wahl, int *anzahl, int i, char *name)
{
	char s[namelength], *posit;

	p = parr[pn];

	if (wahl == 1)
		*anzahl = p->nparameter;
	else {
		if (!p->nparameter)
			*name = '\0';
		else {
			strcpy(s, p->lexicon[p->parameter[i]].nm);

			posit = (char *)strstr(s, " ");
			if (posit !=NULL) {
				strncpy(name, s, (posit - s));
				*(name + (posit - s)) = '\0';
			}
			else strcpy(name, s);
		}
	}
}                             /* cparinfo */


static int perm[NEXPRESSIONS][10];

int cpermut(int expid, int parid)
{
    return perm[expid][parid];    
}                             /* cpermut */


int cnewexprx(int p, char *str, int nparas, char *para[])
{
   char ausdr[50];
   char parname[10][10];
   int i, j, r;
   int paranzahl;
   
   if (p < 0 || p > 9) {   
      printf("Expression-Identifier falsch\n");
	  return 1;
   }

   printf("%s", str);
   gets(ausdr);
   r = cinitausd(p, ausdr);
   if (r == 1) 
   {
      printf("Ausdruck fehlerhaft\n");
	  return 1;
   }

   cparinfo(p, 1, &paranzahl, 1, parname[0]);
/*
   printf("Zahl der Parameter = %d\n", paranzahl);
*/
 
   for (i = 1; i <= paranzahl; i++) {
      cparinfo(p, 2, &paranzahl, i, parname[i]);
/*    printf("Parameter Nr. %d = %s\n", i, parname[i]); */
      j = 1; 
      do { 
         if (strcmp(parname[i], para[j]) == 0) break;
         j++;
      } while (j <= nparas);  
      if (j > nparas) { 
         printf("Falscher Paramater in Ausdruck \n");
         return 2; 
      }
   }   
   
   for (i = 1; i <= nparas; i++) perm[p][i] = 0;   
   for (i = 1; i <= nparas; i++) {
      for (j = 1; j <= paranzahl; j++) 
         if (strcmp(para[i], parname[j]) == 0)
            perm[p][i] = j;    
   }

/*   
   for (i = 1; i <= nparas; i++) {
      printf("perm[%d][%d] = %d\n", p, i, perm[p][i]);
   }
*/   
   return 0;     
}                          /* cnewexprx */


int cnewexpr(int p, char *str, char *para)
{
   char ausdr[100];
   char parname[1][10];
   int r;
   int paranzahl;
   
   if (p < 0 || p > 9) 
   {
      printf("Expression-Identifier falsch\n");
	  return 1;
   }

   printf("%s", str);
   gets(ausdr);
   r = cinitausd(p, ausdr);
   if (r == 1) 
   {
      printf("Ausdruck fehlerhaft\n");
	  return 1;
   }

   cparinfo(p, 1, &paranzahl, 1, parname[0]);
   if (paranzahl != 1) 
   {  
      printf("Zahl der Parameter ungleich 1\n");
	  return 1;
   }
   
   cparinfo(p, 2, &paranzahl, 1, parname[0]);
   if (strcmp(parname[0], para) != 0) 
   {  
      printf("Parametername falsch\n");
	  return 1;
   }
   
   return 0;     
}                           /* cnewexpr */


int cnewexprs(int p, int nparas, char *para[], char *arg)
{
   char ausdr[100];
   char parname[10][10];
   int i, j, r;
   int paranzahl;
   
   if (p < 0 || p > 9) {   
      printf("Expression-Identifier falsch\n");
	  return 1;
   }

   strcpy(ausdr, arg);
   
   r = cinitausd(p, ausdr);
   if (r == 1) 
   {
      printf("Ausdruck fehlerhaft\n");
	  return 1;
   }

   cparinfo(p, 1, &paranzahl, 1, parname[0]);
/*
   printf("Zahl der Parameter = %d\n", paranzahl);
*/
 
   for (i = 1; i <= paranzahl; i++) {
      cparinfo(p, 2, &paranzahl, i, parname[i]);
/*    printf("Parameter Nr. %d = %s\n", i, parname[i]); */
      j = 1; 
      do { 
         if (strcmp(parname[i], para[j]) == 0) break;
         j++;
      } while (j <= nparas);  
      if (j > nparas) { 
         printf("Falscher Paramater in Ausdruck \n");
         return 2; 
      }
   }   
   
   for (i = 1; i <= nparas; i++) perm[p][i] = 0;   
   for (i = 1; i <= nparas; i++) {
      for (j = 1; j <= paranzahl; j++) 
         if (strcmp(para[i], parname[j]) == 0)
            perm[p][i] = j;    
   }

/*   
   for (i = 1; i <= nparas; i++) {
      printf("perm[%d][%d] = %d\n", p, i, perm[p][i]);
   }
*/   
   return 0;     
}                          /* cnewexprs */


void cinstrukt()
{
    printf("\nMoegliche Operatoren des Formelparsers sind:\n");
    printf("  Operator     Prioritaet\n");
    printf("   +, -           4 \n");
    printf("   *, /           5 \n");
//  printf("   div, mod       5 \n");
//  printf("   ^              6 \n");
    printf("   - (unaer)      6 \n");
    printf("   abs            6 \n");
    printf("   sqr, sqrt      6 \n");
    printf("   exp, log       6 \n");
    printf("   lg (Basis 2)   6 \n");
    printf("   sin, cos       6 \n");
	printf("   tan, cot       6 \n");
	printf("   asin, acos     6 \n");
	printf("   atan, acot     6 \n");
    printf("   sinh, cosh     6 \n");
	printf("   tanh, coth     6 \n");
//  printf("   round, trunc   6 \n");
//  printf("   u (unit step)  6 \n");
    printf("Bekannte Konstanten sind:\n");
    printf("   e, pi, pio2      \n");
    printf("   pisq, log10, i \n\n");
	return;
}                             /* cinstrukt */

int calc_error;
