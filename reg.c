/* This is a plain and simple linear regression solver.

To build and run:

     gcc reg.c -lm
     ./a.out infile

Infile has the following format:

1. The first line has a model statement such as Y = A B A*B
   An asterisk indicates an interaction term or multiplication.
   Powers of numerical variables are written as A*A, etc.
   The base-ten log of numerical variable A is written log(A).

2. The second line has a list of variable names in the same order as the data.
   A dollar sign is needed after each categorical variable.

3. All of the remaining lines are data lines.

4. Blank lines and lines that start with an asterisk are skipped.

For example, the infile

log(Volume) = log(Girth) log(Height)
Girth Height Volume
  8.3     70   10.3
  8.6     65   10.3
  8.8     63   10.2
 10.5     72   16.4
 10.7     81   18.8
 10.8     83   19.7
 11.0     66   15.6
 11.0     75   18.2
 11.1     80   22.6
 11.2     75   19.9
 11.3     79   24.2
 11.4     76   21.0
 11.4     76   21.4
 11.7     69   21.3
 12.0     75   19.1
 12.9     74   22.2
 12.9     85   33.8
 13.3     86   27.4
 13.7     71   25.7
 13.8     64   24.9
 14.0     78   34.5
 14.2     80   31.7
 14.5     74   36.3
 16.0     72   38.3
 16.3     77   42.6
 17.3     81   55.4
 17.5     82   55.7
 17.9     80   58.3
 18.0     80   51.5
 18.0     80   51.0
 20.6     87   77.0

produces the following output:

                              Analysis of Variance

    Source     DF     Sum of Squares     Mean Square     F Value     Pr > F
    Model       2         1.53213547      0.76606773      613.19     0.0000
    Error      28         0.03498056      0.00124931                       
    Total      30         1.56711603                                       

               Root MSE           0.03535     R-Square     0.9777
               Dependent Mean     1.42133     Adj R-Sq     0.9761
               Coeff Var          2.48679                        

                              Parameter Estimates

         Parameter       Estimate     Std Err     t Value     Pr > |t|
         (Intercept)     -2.88007     0.34734       -8.29       0.0000
         log(Girth)       1.98265     0.07501       26.43       0.0000
         log(Height)      1.11712     0.20444        5.46       0.0000
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#define __USE_ISOC99
#include <math.h>

#ifndef NAN
#define NAN nan("0")
#endif

#define MAXOBS 10000
#define MAXVAR 20
#define MAXLVL 20
#define MAXPAR 100
#define MAXINT 4

//	data	Input data array
//
//	nobs	Number of observations (rows in data array)
//
//	nvar	Number of variables (columns in data array)
//
//	nexp	Number of explanatory variables in the model
//
//	npar	Number of columns in the design matrix X

double data[MAXOBS][MAXVAR];
int nobs;
int nvar;
int nexp;
int npar;

//	B	Regression coefficients
//
//	G	Inverse of X'X
//
//	T	Temporary matrix
//
//	X	Design matrix
//
//	Y	Response vector
//
//	Yhat	Predicted response X * B

double B[MAXPAR];
double G[MAXPAR][MAXPAR];
double T[MAXPAR][MAXPAR];
double X[MAXOBS][MAXPAR];
double Y[MAXOBS];
double Yhat[MAXOBS];

//	alpha	Level of significance
//
//	css	Corrected sum of squares
//
//	df	Degrees of freedom
//
//	dfe	Degrees of freedom error
//
//	fval	Summary F-statistic
//
//	mse	Mean square error (estimate of model variance)
//
//	msr	Mean square regression
//
//	pval	p-value for F-statistic
//
//	ss	Sequential sum of squares (Type I)
//
//	ssr	Sum of squares regression
//
//	sse	Sum of squares error
//
//	ybar	Mean of response variable Y

double alpha = 0.05;
double css;
int df[MAXVAR];
double dfe;
double fval;
double mse;
double msr;
double pval;
double ss[MAXVAR];
double ssr;
double sse;
double ybar;

//	adjrsq	Adjusted R-squared
//
//	cv	Coefficient of variation
//
//	rsq	R-squared
//
//	sqrtmse	Square root of MSE

double adjrsq;
double cv;
double rsq;
double sqrtmse;

// parameter table

struct {
	double est;
	double se;
	double tval;
	double pval;
} param[MAXPAR];

int intercept = 1;
char *var_name[MAXVAR];
int var_type[MAXVAR];
char *level_name[MAXVAR][MAXLVL];
int level_count[MAXVAR];
#define MAXTOK 100
char *ttab[MAXTOK];
char **token;
char *model;
int response_index;
int explanatory_index[MAXVAR];
char buf[10000];
int nrow;
int ncol;
#define MAXROW 1000
#define MAXCOL 10
int fmt[MAXCOL];
char *atab[MAXROW][MAXCOL];
#define A(i, j) atab[i][j]

void read_infile(char *filename);
int get_next_line(FILE *f);
void scan_var_names();
void tokenize(char *s);
int isname(int c);
void scan_obs();
void scan_number(char *s, int v);
void scan_level_name(char *s, int v);
void parse_model();
int parse_model_term();
int parse_log();
int parse_interaction();
void numerical_interaction(int vtab[], int n);
void treatment_interaction(int vtab[], int n);
int get_var_index(char *s);
void stop(char *s);
void stop_and_print_buf(char *s);
void compute();
void prelim();
void compute_model_parameters();
void fit_variable(int v);
void numerical_column(int v);
void indicator_column(int v, int k);
int compute_G();
void compute_B();
void compute_Yhat();
void compute_ss();
void compute_param_table();
void print_results();
void print_anova_table_part1();
void print_anova_table_part2();
void print_param_table();
void print_table_and_free();
void print_table();
void emit_line(char *s);
void emit_line_center(char *s);
double qt(double p, double df);
double tdist(double t, double df);
double qf(double p, double df1, double df2);
double fdist(double t, double df1, double df2);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);

int
main(int argc, char *argv[])
{
	if (argc < 2)
		stop("Input file?");
	read_infile(argv[1]);
	compute();
	print_results();
}

void
read_infile(char *filename)
{
	FILE *f;
	f = fopen(filename, "r");
	if (f == NULL)
		stop("Cannot open file");
	if (get_next_line(f) == 0)
		stop("File format?");
	model = strdup(buf);
	if (get_next_line(f) == 0)
		stop("File format?");
	scan_var_names();
	while (get_next_line(f))
		scan_obs();
	fclose(f);
	parse_model();
}

// Skip blank lines and comment lines

int
get_next_line(FILE *f)
{
	char *s;
	while (fgets(buf, sizeof buf, f)) {
		s = buf;
		while (*s && *s <= ' ')
			s++;
		if (*s && *s != '*')
			return 1;
	}
	return 0;
}

void
scan_var_names()
{
	tokenize(buf);
	while (**token) {
		if (nvar == MAXVAR)
			stop("Try increasing MAXVAR");
		var_name[nvar] = strdup(*token++);
		if (strcmp(*token, "$") == 0) {
			var_type[nvar] = 1;
			token++;
		}
		nvar++;
	}
}

void
tokenize(char *s)
{
	int c, i = 0;
	char *t;
	for (;;) {
		while (*s && *s <= ' ')
			s++;
		if (*s == 0)
			break;
		if (i == MAXTOK - 1)
			stop("Try increasing MAXTOK");
		t = s;
		if (isname(*t))
			while (isname(*++t));
		else
			t++;
		c = *t;
		*t = 0;
		if (ttab[i])
			free(ttab[i]);
		ttab[i++] = strdup(s);
		*t = c;
		s = t;
	}
	if (ttab[i])
		free(ttab[i]);
	ttab[i] = strdup("");
	token = ttab;
}

int
isname(int c)
{
	return (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
}

// Scan one line of observational data

void
scan_obs()
{
	int i;
	char *s;

	if (nobs == MAXOBS)
		stop("Try increasing MAXOBS");

	s = buf;

	for (i = 0; i < nvar; i++) {

		// skip leading space

		while (*s && *s <= ' ')
			s++;

		if (*s == 0)
			stop_and_print_buf("Missing data?");

		switch (var_type[i]) {
		case 0:
			scan_number(s, i);
			break;
		case 1:
			scan_level_name(s, i);
			break;
		}

		// skip over datum

		while (*s > ' ' && *s != ',' && *s != ';')
			s++;

		// skip trailing space

		while (*s && *s <= ' ')
			s++;

		// skip delimiter

		if (*s == ',' || *s == ';')
			s++;
	}

	nobs++;
}

void
scan_number(char *s, int v)
{
	if (sscanf(s, "%lf", &data[nobs][v]) < 1)
		stop_and_print_buf("Number format?");
}

void
scan_level_name(char *s, int v)
{
	int c, i, k;
	char *t;

	t = s;

	while (*t > ' ' && *t != ',' && *t != ';')
		t++;

	if (t == s)
		stop_and_print_buf("Missing data?");

	c = *t;
	*t = 0;

	k = level_count[v];

	for (i = 0; i < k; i++)
		if (strcmp(level_name[v][i], s) == 0)
			break;

	if (i == k) {
		if (k == MAXLVL)
			stop("Try increasing MAXLVL");
		level_name[v][k] = strdup(s);
		level_count[v]++;
	}

	*t = c;

	data[nobs][v] = i;
}

void
parse_model()
{
	tokenize(model);
	response_index = parse_model_term();
	if (var_type[response_index] != 0)
		stop("Response variable is not numerical");
	if (strcmp(*token++, "=") != 0)
		stop("Missing or misplaced '=' in model statement");
	while (**token) {
		if (nexp == MAXVAR)
			stop("Try increasing MAXVAR");
		explanatory_index[nexp++] = parse_model_term();
	}
}

int
parse_model_term()
{
	if (strcmp(token[0], "log") == 0 && strcmp(token[1], "(") == 0)
		return parse_log();
	else if (strcmp(token[1], "*") == 0)
		return parse_interaction();
	else
		return get_var_index(*token++);
}

int
parse_log()
{
	int i, v;

	if (nvar == MAXVAR)
		stop("Try increasing MAXVAR");

	v = get_var_index(token[2]);

	if (var_type[v] != 0)
		stop("Log of categorical variable?");

	if (strcmp(token[3], ")") != 0)
		stop("Missing or misplaced ')' in model statement");

	// name of new variable

	strcpy(buf, "log(");
	strcat(buf, token[2]);
	strcat(buf, ")");

	token += 4;

	var_name[nvar] = strdup(buf);
	var_type[nvar] = 0;
	
	for (i = 0; i < nobs; i++)
		data[i][nvar] = log10(data[i][v]);

	return nvar++;
}

int
parse_interaction()
{
	int i, n = 0, vtab[MAXINT];

	if (nvar == MAXVAR)
		stop("Try increasing MAXVAR");

	// stitch together the name of the interaction variable

	strcpy(buf, *token);
	vtab[n++] = get_var_index(*token++);

	while (strcmp(*token, "*") == 0) {
		if (n == MAXINT)
			stop("Try increasing MAXINT");
		strcat(buf, *token++);
		strcat(buf, *token);
		vtab[n++] = get_var_index(*token++);
	}

	var_name[nvar] = strdup(buf);
	var_type[nvar] = var_type[vtab[0]];

	// check that var types match

	for (i = 1; i < n; i++)
		if (var_type[vtab[0]] != var_type[vtab[i]])
			stop("Mixed var types in interaction term");

	if (var_type[vtab[0]] == 0)
		numerical_interaction(vtab, n);
	else
		treatment_interaction(vtab, n);

	return nvar++;
}

void
numerical_interaction(int vtab[], int n)
{
	int i, j, v;
	double x;
	for (i = 0; i < nobs; i++) {
		v = vtab[0];
		x = data[i][v];
		for (j = 1; j < n; j++) {
			v = vtab[j];
			x *= data[i][v];
		}
		data[i][nvar] = x;
	}
}

void
treatment_interaction(int vtab[], int n)
{
	int c, count[MAXINT], i, j, k, v;
	char *s;

	// calculate the overall number of levels

	k = 1;

	for (i = 0; i < n; i++) {
		v = vtab[i];
		if (var_type[v] == 0)
			stop("Type of variable?");
		k *= level_count[v];
	}

	if (k >= MAXLVL)
		stop("Try increasing MAXLVL");

	level_count[nvar] = k;

	for (i = 0; i < n; i++)
		count[i] = 0;

	for (i = 0; i < k; i++) {

		*buf = 0;

		// stitch together the level name

		for (j = 0; j < n; j++) {
			v = vtab[j];
			c = count[j];
			s = level_name[v][c];
			strcat(buf, s);
			if (j < n - 1)
				strcat(buf, "*");
		}

		level_name[nvar][i] = strdup(buf);

		// increment level indices

		for (j = n - 1; j >= 0; j--) {
			v = vtab[j];
			if (++count[j] < level_count[v])
				break;
			count[j] = 0;
		}
	}

	// add data column

	for (i = 0; i < nobs; i++) {
		k = 0;
		for (j = 0; j < n; j++) {
			v = vtab[j];
			k += data[i][v];
			if (j < n - 1) {
				v = vtab[j + 1];
				k *= level_count[v];
			}
		}
		data[i][nvar] = k;
	}
}

int
get_var_index(char *s)
{
	int i;
	for (i = 0; i < nvar; i++)
		if (strcmp(s, var_name[i]) == 0)
			return i;
	stop("Misplaced or incorrect model variable");
	return 0;
}

void
stop(char *s)
{
	printf("%s\n", s);
	exit(1);
}

void
stop_and_print_buf(char *s)
{
	printf("%s\n", s);
	printf("%s\n", buf);
	exit(1);
}

void
compute()
{
	double df1, df2;
	prelim();
	compute_model_parameters();
	dfe = nobs - npar;
	mse = sse / dfe;
	sqrtmse = sqrt(mse);
	msr = ssr / (npar - 1);
	fval = msr / mse;
	df1 = npar - 1;
	df2 = nobs - npar;
	pval = 1.0 - fdist(fval, df1, df2);
	rsq = 1.0 - sse / css;
	adjrsq = 1.0 - (nobs - 1.0) / (nobs - npar) * sse / css;
	cv = 100.0 * sqrtmse / ybar;
	compute_param_table();
}

void
prelim()
{
	int i;

	// mean of response variable

	ybar = 0;

	for (i = 0; i < nobs; i++) {
		Y[i] = data[i][response_index];
		ybar += (Y[i] - ybar) / (i + 1);
	}

	// corrected sum of squares

	css = 0;

	for (i = 0; i < nobs; i++)
		css += (Y[i] - ybar) * (Y[i] - ybar);
}

void
compute_model_parameters()
{
	int i, v;

	if (intercept) {
		npar = 1;
		for (i = 0; i < nobs; i++)
			X[i][0] = 1;
	}

	// fit explanatory variables

	for (i = 0; i < nexp; i++) {
		v = explanatory_index[i];
		fit_variable(v);
		df[i] = npar;
		compute_B();
		compute_Yhat();
		compute_ss();
		ss[i] = ssr;
	}

	// convert cumulative df[] and ss[] to differential

	for (i = nexp - 1; i > 0; i--) {
		df[i] -= df[i - 1];
		ss[i] -= ss[i - 1];
	}

	if (intercept)
		df[0]--; // subtract 1 for intercept
}

void
fit_variable(int v)
{
	int k, n;
	switch (var_type[v]) {
	case 0:
		numerical_column(v);
		break;
	case 1:
		n = level_count[v];
		for (k = 0; k < n; k++)
			indicator_column(v, k);
		break;
	}
}

void
numerical_column(int v)
{
	int i;
	if (npar == MAXPAR)
		stop("Try increasing MAXPAR");
	for (i = 0; i < nobs; i++)
		X[i][npar] = data[i][v];
	npar++;
	if (compute_G() == -1) {
		npar--;		// X'X is singular, remove column
		compute_G();	// restore previous G
	}
}

void
indicator_column(int v, int k)
{
	int i;
	if (npar == MAXPAR)
		stop("Try increasing MAXPAR");
	for (i = 0; i < nobs; i++)
		if (data[i][v] == k)
			X[i][npar] = 1;
		else
			X[i][npar] = 0;
	npar++;
	if (compute_G() == -1) {
		npar--;		// X'X is singular, remove column
		compute_G();	// restore previous G
	}
}

// G = inv X'X

int
compute_G()
{
	int d, i, j, k;
	double a, b, m, t;

	a = 0;
	b = 0;

	// G = I

	for (i = 0; i < npar; i++)
		for (j = 0; j < npar; j++)
			if (i == j)
				G[i][j] = 1;
			else
				G[i][j] = 0;

	// T = X'X

	for (i = 0; i < npar; i++)
		for (j = 0; j < npar; j++) {
			t = 0;
			for (k = 0; k < nobs; k++)
				t += X[k][i] * X[k][j];
			T[i][j] = t;
		}

	// G = inv T

	for (d = 0; d < npar; d++) {

		// find the best pivot row

		k = d;
		for (i = d + 1; i < npar; i++)
			if (fabs(T[i][d]) > fabs(T[k][d]))
				k = i;

		// exchange rows if necessary

		if (k != d) {
			for (j = d; j < npar; j++) { // skip zeroes, start at d
				t = T[d][j];
				T[d][j] = T[k][j];
				T[k][j] = t;
			}
			for (j = 0; j < npar; j++) {
				t = G[d][j];
				G[d][j] = G[k][j];
				G[k][j] = t;
			}
		}

		// multiply the pivot row by 1 / pivot

		m = T[d][d];

		if (m == 0)
			return -1; // singular

		if (fabs(m) > a)
			a = fabs(m);

		m = 1 / m;

		if (fabs(m) > b)
			b = fabs(m);

		for (j = d; j < npar; j++) // skip zeroes, start at d
			T[d][j] *= m;
		for (j = 0; j < npar; j++)
			G[d][j] *= m;

		// clear out column below d

		for (i = d + 1; i < npar; i++) {
			m = -T[i][d];
			for (j = d; j < npar; j++) // skip zeroes, start at d
				T[i][j] += m * T[d][j];
			for (j = 0; j < npar; j++)
				G[i][j] += m * G[d][j];
		}
	}

	// clear out columns above diagonal

	for (d = npar - 1; d > 0; d--)
		for (i = 0; i < d; i++) {
			m = -T[i][d];
			for (j = 0; j < npar; j++)
				G[i][j] += m * G[d][j];
		}

	// check ratio of biggest divisor to smallest divisor

	// domain of ab is [1, inf)

	// printf("cond = %g\n", a * b);

	if (a * b < 1e10)
		return 0;	// ok
	else
		return -1;	// singular
}

// B = G * X' * Y

void
compute_B()
{
	int i, j;
	double t;

	// T = X' * Y

	for (i = 0; i < npar; i++) {
		t = 0;
		for (j = 0; j < nobs; j++)
			t += X[j][i] * Y[j];
		T[0][i] = t;
	}

	// B = G * T

	for (i = 0; i < npar; i++) {
		t = 0;
		for (j = 0; j < npar; j++)
			t += G[i][j] * T[0][j];
		B[i] = t;
	}
}

// Yhat = X * B

void
compute_Yhat()
{
	int i, j;
	double t;
	for (i = 0; i < nobs; i++) {
		t = 0;
		for (j = 0; j < npar; j++)
			t += X[i][j] * B[j];
		Yhat[i] = t;
	}
}

void
compute_ss()
{
	int i;
	ssr = 0;
	sse = 0;
	for (i = 0; i < nobs; i++) {
		ssr += (Yhat[i] - ybar) * (Yhat[i] - ybar);
		sse += (Y[i] - Yhat[i]) * (Y[i] - Yhat[i]);
	}
}

void
compute_param_table()
{
	int i;
	double df, est, se, tval, pval;
	df = nobs - npar;
	for (i = 0; i < npar; i++) {
		est = B[i];
		se = sqrt(mse * G[i][i]);
		tval = est / se;
		pval = 2.0 * (1.0 - tdist(fabs(tval), df));
		param[i].est = est;
		param[i].se = se;
		param[i].tval = tval;
		param[i].pval = pval;
	}
}

void
print_results()
{
	int i, v;
	emit_line_center("Analysis of Variance");
	emit_line("");
	print_anova_table_part1();
	print_anova_table_part2();
	emit_line_center("Parameter Estimates");
	emit_line("");
	print_param_table();
}

// Example:
//
//    Source     DF     Sum of Squares     Mean Square     F Value     Pr > F
//    Model       3        91.71428571     30.57142857       15.29     0.0253
//    Error       3         6.00000000      2.00000000                       
//    Total       6        97.71428571                                       

void
print_anova_table_part1()
{
	nrow = 4;
	ncol = 6;

	// 1st row

	A(0, 0) = strdup("Source");
	A(0, 1) = strdup("DF");
	A(0, 2) = strdup("Sum of Squares");
	A(0, 3) = strdup("Mean Square");
	A(0, 4) = strdup("F Value");
	A(0, 5) = strdup("Pr > F");

	// 2nd row

	A(1, 0) = strdup("Model");

	sprintf(buf, "%d", npar - 1);
	A(1, 1) = strdup(buf);

	sprintf(buf, "%0.8f", ssr);
	A(1, 2) = strdup(buf);

	sprintf(buf, "%0.8f", msr);
	A(1, 3) = strdup(buf);

	sprintf(buf, "%0.2f", fval);
	A(1, 4) = strdup(buf);

	sprintf(buf, "%0.4f", pval);
	A(1, 5) = strdup(buf);

	// 3rd row

	A(2, 0) = strdup("Error");

	sprintf(buf, "%d", nobs - npar);
	A(2, 1) = strdup(buf);

	sprintf(buf, "%0.8f", sse);
	A(2, 2) = strdup(buf);

	sprintf(buf, "%0.8f", mse);
	A(2, 3) = strdup(buf);

	A(2, 4) = strdup("");
	A(2, 5) = strdup("");

	// 4th row

	A(3, 0) = strdup("Total");

	sprintf(buf, "%d", nobs - 1);
	A(3, 1) = strdup(buf);

	sprintf(buf, "%0.8f", css);
	A(3, 2) = strdup(buf);

	A(3, 3) = strdup("");
	A(3, 4) = strdup("");
	A(3, 5) = strdup("");

	fmt[0] = 1; // left justify
	fmt[1] = 0;
	fmt[2] = 0;
	fmt[3] = 0;
	fmt[4] = 0;
	fmt[5] = 0;

	print_table_and_free();
}

// Example:
//
//    Root MSE            0.64801     R-Square     0.3606
//    Dependent Mean      5.63602     Adj R-Sq     0.3561
//    Coeff Var          11.49767                        

void
print_anova_table_part2()
{
	nrow = 3;
	ncol = 4;

	A(0, 0) = strdup("Root MSE");
	A(1, 0) = strdup("Dependent Mean");
	A(2, 0) = strdup("Coeff Var");

	A(0, 2) = strdup("R-Square");
	A(1, 2) = strdup("Adj R-Sq");
	A(2, 2) = strdup("");

	sprintf(buf, "%0.5f", sqrtmse);
	A(0, 1) = strdup(buf);

	sprintf(buf, "%0.5f", ybar);
	A(1, 1) = strdup(buf);

	sprintf(buf, "%0.5f", cv);
	A(2, 1) = strdup(buf);

	sprintf(buf, "%0.4f", rsq);
	A(0, 3) = strdup(buf);

	sprintf(buf, "%0.4f", adjrsq);
	A(1, 3) = strdup(buf);

	A(2, 3) = strdup("");

	fmt[0] = 1; // left justify
	fmt[1] = 0;
	fmt[2] = 1;
	fmt[3] = 0;

	print_table_and_free();
}

void
print_param_table()
{
	int i, j, k, m, n, v;

	// determine number of rows in table

	nrow = 1;
	ncol = 5;

	if (intercept)
		nrow++;

	for (i = 0; i < nexp; i++) {
		v = explanatory_index[i];
		switch (var_type[v]) {
		case 0:
			nrow++;
			break;
		case 1:
			nrow += level_count[v];
			break;
		}
	}

	if (nrow >= MAXROW)
		stop("Try increasing MAXROW");

	// column names

	A(0, 0) = strdup("Parameter");
	A(0, 1) = strdup("Estimate");
	A(0, 2) = strdup("Std Err");
	A(0, 3) = strdup("t Value");
	A(0, 4) = strdup("Pr > |t|");

	// variable names

	j = 1;

	if (intercept)
		A(j++, 0) = strdup("(Intercept)");

	for (i = 0; i < nexp; i++) {
		v = explanatory_index[i];
		switch (var_type[v]) {
		case 0:
			A(j++, 0) = strdup(var_name[v]);
			break;
		case 1:
			n = level_count[v];
			for (k = 0; k < n; k++) {
				sprintf(buf, "%s %s", var_name[v], level_name[v][k]);
				A(j++, 0) = strdup(buf);
			}
		}
	}

	// estimates

	j = 1;
	m = 0;

	if (intercept) {
		sprintf(buf, "%0.5f", param[m].est);
		A(j, 1) = strdup(buf);
		sprintf(buf, "%0.5f", param[m].se);
		A(j, 2) = strdup(buf);
		sprintf(buf, "%0.2f", param[m].tval);
		A(j, 3) = strdup(buf);
		sprintf(buf, "%0.4f", param[m].pval);
		A(j, 4) = strdup(buf);
		j++;
		m++;
	}

	for (i = 0; i < nexp; i++) {

		v = explanatory_index[i];

		switch (var_type[v]) {
		case 0:
			if (df[i] == 0) {
				A(j, 1) = strdup(".");
				A(j, 2) = strdup(".");
				A(j, 3) = strdup(".");
				A(j, 4) = strdup(".");
			} else {
				sprintf(buf, "%0.5f", param[m].est);
				A(j, 1) = strdup(buf);
				sprintf(buf, "%0.5f", param[m].se);
				A(j, 2) = strdup(buf);
				sprintf(buf, "%0.2f", param[m].tval);
				A(j, 3) = strdup(buf);
				sprintf(buf, "%0.4f", param[m].pval);
				A(j, 4) = strdup(buf);
				m++;
			}
			j++;
			break;

		case 1:
			n = level_count[v];
			for (k = 0; k < n; k++) {
				if (df[i] <= k) {
					A(j, 1) = strdup(".");
					A(j, 2) = strdup(".");
					A(j, 3) = strdup(".");
					A(j, 4) = strdup(".");
				} else {
					sprintf(buf, "%0.5f", param[m].est);
					A(j, 1) = strdup(buf);
					sprintf(buf, "%0.5f", param[m].se);
					A(j, 2) = strdup(buf);
					sprintf(buf, "%0.2f", param[m].tval);
					A(j, 3) = strdup(buf);
					sprintf(buf, "%0.4f", param[m].pval);
					A(j, 4) = strdup(buf);
					m++;
				}
				j++;
			}
			break;
		}
	}

	fmt[0] = 1; // right justify
	fmt[1] = 0;
	fmt[2] = 0;
	fmt[3] = 0;
	fmt[4] = 0;

	print_table_and_free();
}

void
print_table_and_free()
{
	int i, j;
	print_table();
	for (i = 0; i < nrow; i++)
		for (j = 0; j < ncol; j++)
			if (A(i, j)) {
				free(A(i, j));
				A(i, j) = NULL;
			}
}

void
print_table()
{
	int c, i, j, k, n, nsp = 0, t;
	char *b, *s;

	static int w[MAXCOL];

	// measure column widths

	t = 0;

	for (j = 0; j < ncol; j++) {
		w[j] = 0;
		for (i = 0; i < nrow; i++) {
			s = A(i, j);
			if (s == NULL)
				continue;
			n = (int) strlen(s);
			if (n > w[j])
				w[j] = n;
		}
		t += w[j];
	}

	// spaces between columns

	if (ncol > 1) {
		nsp = (80 - t) / (ncol - 1);
		if (nsp < 2)
			nsp = 2;
		if (nsp > 5)
			nsp = 5;
		t += (ncol - 1) * nsp;
	}

	// number of spaces to center the line

	if (t < 80)
		c = (80 - t) / 2;
	else
		c = 0;

	// leading spaces

	memset(buf, ' ', c);

	// for each row

	for (i = 0; i < nrow; i++) {

		b = buf + c;

		// for each column

		for (j = 0; j < ncol; j++) {

			// space between columns

			if (j > 0) {
				memset(b, ' ', nsp);
				b += nsp;
			}

			s = A(i, j);

			if (s == NULL) {
				memset(b, ' ', w[j]);
				b += w[j];
			} else {

				n = (int) strlen(s);

				k = w[j] - n;

				switch (fmt[j]) {

				// right aligned

				case 0:
					memset(b, ' ', k);
					b += k;
					strcpy(b, s);
					b += n;
					break;

				// left aligned

				case 1:
					strcpy(b, s);
					b += n;
					memset(b, ' ', k);
					b += k;
					break;
				}
			}
		}

		*b = 0;

		emit_line(buf);
	}

	emit_line("");
}

void
emit_line(char *s)
{
	printf("%s\n", s);
}

void
emit_line_center(char *s)
{
	int i, n;
	n = strlen(s);
	n = (80 - n) / 2;
	for (i = 0; i < n; i++)
		printf(" ");
	printf("%s\n", s);
}

// t-distribution quantile function, like qt() in R

double
qt(double p, double df)
{
	int i;
	double a, t, t1, t2;
	if (isnan(p) || isnan(df) || df < 1.0)
		return NAN;
	t1 = -1000.0;
	t2 = 1000.0;
	for (i = 0; i < 100; i++) {
		t = 0.5 * (t1 + t2);
		a = tdist(t, df);
		if (fabs(a - p) < 1e-10)
			break;
		if (a < p)
			t1 = t;
		else
			t2 = t;
	}
	return t;
}

// t-distribution cdf, like pt() in R

double
tdist(double t, double df)
{
	double a;
	if (isnan(t) || isnan(df) || df < 1.0)
		return NAN;
	a = 0.5 * betai(0.5 * df, 0.5, df / (df + t * t));
	if (t > 0.0)
		a = 1.0 - a;
	return a;
}

// F-distribution quantile function, like qf() in R

double
qf(double p, double df1, double df2)
{
	int i;
	double a, t, t1, t2;
	if (isnan(p) || isnan(df1) || isnan(df2) || df1 < 1.0 || df2 < 1.0)
		return NAN;
	t1 = 0.0;
	t2 = 1000.0;
	for (i = 0; i < 100; i++) {
		t = 0.5 * (t1 + t2);
		a = fdist(t, df1, df2);
		if (fabs(a - p) < 1e-10)
			break;
		if (a < p)
			t1 = t;
		else
			t2 = t;
	}
	return t;
}

// F-distribution cdf, like pf() in R

double
fdist(double t, double df1, double df2)
{
	if (isnan(t) || isnan(df1) || isnan(df2) || df1 < 1.0 || df2 < 1.0)
		return NAN;
	if (t < 0.0)
		return 0.0;
	else
		return betai(0.5 * df1, 0.5 * df2, t / (t + df2 / df1));
}

// From Numerical Recipes in C, p. 214

double
gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// From Numerical Recipes in C, p. 227

double
betai(double a, double b, double x)
{
	double bt;
	if (x < 0.0 || x > 1.0) return NAN;
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

// From Numerical Recipes in C, p. 227

#define FPMIN 1.0e-30

double
betacf(double a, double b, double x)
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=100;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < 1e-10) break;
	}
	return h;
}
