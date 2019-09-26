/* This is a plain and simple anova solver.

(Recall that for anova models the explanatory variables are all categorical.)

To build and run:

     gcc anova.c -lm
     ./a.out infile

Infile has the following format:

1. The first line has a model statement such as Y = A B A*B
   An asterisk indicates an interaction term.

2. The second line has a list of variable names in the same order as the data.
   A dollar sign is needed after each explanatory variable.

3. All of the remaining lines are data lines.

4. Blank lines and lines that start with an asterisk are skipped.

For example, the infile

Y = A B A*B
A$ B$ Y
A1 B1 12
A1 B1 14
A1 B2 11
A1 B2 9
A2 B1 20
A2 B1 18
A2 B2 17

produces the following output:

                              Analysis of Variance

    Source     DF     Sum of Squares     Mean Square     F Value     Pr > F
    Model       3        91.71428571     30.57142857       15.29     0.0253
    Error       3         6.00000000      2.00000000                       
    Total       6        97.71428571                                       

               R-Square     Coeff Var     Root MSE        Y Mean
               0.938596      9.801480     1.414214     14.428571

      Source     DF        Anova SS     Mean Square     F Value     Pr > F
      A           1     80.04761905     80.04761905       40.02     0.0080
      B           1     11.26666667     11.26666667        5.63     0.0982
      A*B         1      0.40000000      0.40000000        0.20     0.6850

                                 Mean Response

              A      N        Mean Y     95% CI MIN     95% CI MAX
              A1     4     11.500000       9.249671      13.750329
              A2     3     18.333333      15.734877      20.931790

                                 Mean Response

              B      N        Mean Y     95% CI MIN     95% CI MAX
              B1     4     16.000000      13.749671      18.250329
              B2     3     12.333333       9.734877      14.931790

Tables of least significant difference tests are also printed. */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#define __USE_ISOC99
#include <math.h>

#ifndef NAN
#define NAN nan("0")
#endif

#define MAXOBS 1000
#define MAXVAR 10
#define MAXLVL 20
#define MAXPAR 100
#define MAXINT 4

//	data	Input data array
//
//	nobs	Number of observations, rows in data[][]
//
//	nvar	Number of variables, columns in data[][]
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

// Example:
//                                 Mean Response
//
//              A      N        Mean Y     95% CI MIN     95% CI MAX
//              A1     4     11.500000       9.249671      13.750329
//              A2     3     18.333333      15.734877      20.931790
//
// For the above table, the following values are computed in compute_means()
//
//	count		Values for N column
//
//	mean		Values for Mean column
//
//	ci_min		Values for CI MIN column
//
//	ci_max		Values for CI MAX column

int count[MAXLVL];
double mean[MAXLVL];
double variance[MAXLVL];
double ci_min[MAXLVL];
double ci_max[MAXLVL];

// Example:
//                  Least Significant Difference Test
//
// A    A      Delta Y   95% CI MIN   95% CI MAX   t Value   Pr > |t|  
// A1   A2   -6.833333   -10.270768    -3.395898     -6.33     0.0080 *
// A2   A1    6.833333     3.395898    10.270768      6.33     0.0080 *
//
// For the above table, the following values are computed in compute_lsd()
//
//	delta		Values for Delta column
//
//	min		Values for CI MIN column
//
//	max		Values for CI MAX column
//
//	tval		Values for t Value column
//
//	pval		Values for Pr > |t| column

struct {
	double delta;
	double min;
	double max;
	double tval;
	double pval;
} lsd[MAXLVL][MAXLVL];

char *var_name[MAXVAR];
int var_type[MAXVAR];
char *level_name[MAXVAR][MAXLVL];
int level_count[MAXVAR];
#define MAXTOK 100
char *token[MAXTOK];
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
int tokenize(char *s);
void scan_obs();
void scan_number(char *s, int v);
void scan_level_name(char *s, int v);
void parse_model();
int parse_interaction_term(int k, int m);
void add_interaction_column(int vtab[], int n);
int get_var_index(char *s);
void stop(char *s);
void stop_and_print_buf(char *s);
void compute();
void prelim();
void fit();
void fit1(int v, int k);
int compute_G();
void compute_B();
void compute_Yhat();
void compute_ss();
void compute_means(int v);
void compute_lsd(int v);
void print_results();
void print_anova_table_part1();
void print_anova_table_part2();
void print_anova_table_part3();
void print_means(int v);
void print_lsd(int v);
void print_ttest(int v);
void check_atab();
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
	int i = 0, n;
	n = tokenize(buf);
	while (i < n) {
		if (nvar == MAXVAR)
			stop("Try increasing MAXVAR");
		var_name[nvar] = token[i++];
		if (i < n && strcmp(token[i], "$") == 0) {
			var_type[nvar] = 1;
			i++;
		}
		nvar++;
	}
}

int
tokenize(char *s)
{
	int c, n = 0;
	char *t;
	for (;;) {
		while (*s && *s <= ' ')
			s++;
		if (*s == 0)
			break;
		if (n == MAXTOK)
			stop("Try increasing MAXTOK");
		t = s;
		if (isalnum(*t))
			while (isalnum(*++t));
		else
			t++;
		c = *t;
		*t = 0;
		token[n++] = strdup(s);
		*t = c;
		s = t;
	}
	return n;
}

// Scan one line of observational data from buf[] and put in data[][]

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
	int k, n, v;
	n = tokenize(model);
	if (n < 3 || strcmp(token[1], "=") != 0)
		stop("Model?");
	response_index = get_var_index(token[0]);
	if (var_type[response_index])
		stop("Type of variable?");
	k = 2;
	while (k < n) {
		if (k < n - 1 && strcmp(token[k + 1], "*") == 0) {
			k = parse_interaction_term(k, n);
			v = nvar++;
		} else {
			v = get_var_index(token[k++]);
			if (var_type[v] == 0)
				stop("Type of variable?");
		}
		explanatory_index[nexp++] = v;
	}
}

int
parse_interaction_term(int k, int m)
{
	int i, n, vtab[MAXINT];

	n = 0;
	*buf = 0;

	// stitch together the name of the interaction variable

	strcat(buf, token[k]);
	vtab[n++] = get_var_index(token[k]);
	k++;

	while (k < m - 1 && strcmp(token[k], "*") == 0) {
		if (n == MAXINT)
			stop("Try increasing MAXINT");
		k++;
		strcat(buf, "*");
		strcat(buf, token[k]);
		vtab[n++] = get_var_index(token[k]);
		k++;
	}

	var_name[nvar] = strdup(buf);
	var_type[nvar] = 1;

	add_interaction_column(vtab, n);

	return k;
}

// Add interaction column
//
//	vtab		Table of indexes into var_name[]
//
//	n		Number of entries in vtab[]

void
add_interaction_column(int vtab[], int n)
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
	stop("Model?");
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
	fit();
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
fit()
{
	int i, k, n, v;

	// put in intercept

	for (i = 0; i < nobs; i++)
		X[i][0] = 1;

	npar = 1;

	// fit explanatory variables

	for (i = 0; i < nexp; i++) {
		v = explanatory_index[i];
		n = level_count[v];
		for (k = 0; k < n; k++)
			fit1(v, k); // fit next column
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

	df[0]--; // subtract 1 for intercept
}

// Add a column in X for treatment level k

void
fit1(int v, int k)
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

// Compute the mean response for each treatment level of variable v

void
compute_means(int v)
{
	int i, level, n;
	double c, se, t, y;

	// n is the number of treatment levels

	n = level_count[v];

	// compute mean and variance for each level

	for (i = 0; i < n; i++) {
		count[i] = 0;
		mean[i] = 0;
		variance[i] = 0;
	}

	for (i = 0; i < nobs; i++) {
		y = data[i][response_index];
		level = data[i][v];
		count[level]++;
		t = mean[level];
		mean[level] += (y - t) / count[level];
		variance[level] += (y - t) * (y - mean[level]);
	}

	// By the way, level variance is not really needed.
	//
	// Instead, MSE is used for confidence intervals and the LSD test.
	//
	// MSE and the corresponding DF have more power than a simple t-test.
	//
	// That is the whole point of linear modeling.
	//
	// But let's go ahead and finish computing level variance anyway.

	for (i = 0; i < n; i++)
		if (count[i] < 2)
			variance[i] = NAN;
		else
			variance[i] /= (count[i] - 1);

	// Confidence interval for the true mean response per treatment level
	//
	//	c	Critical value
	//
	//	se	Standard error

	c = qt(1 - alpha / 2, dfe);

	for (i = 0; i < n; i++) {
		se = sqrt(mse / count[i]);
		ci_min[i] = mean[i] - c * se;
		ci_max[i] = mean[i] + c * se;
	}
}

// Least significant difference test for treatment v

void
compute_lsd(int v)
{
	int i, j, n;
	double c, d, max, min, pval, se, tval;

	//	n	Number of treatment levels
	//
	//	c	Critical value
	//
	//	se	Standard error

	n = level_count[v];

	c = qt(1 - alpha / 2, dfe);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {

			if (i == j)
				continue;

			d = mean[i] - mean[j];
			se = sqrt(mse * (1.0 / count[i] + 1.0 / count[j]));
			min = d - c * se;
			max = d + c * se;
			tval = d / se;
			pval = 2 * (1 - tdist(fabs(tval), dfe));

			lsd[i][j].delta = d;
			lsd[i][j].min = min;
			lsd[i][j].max = max;
			lsd[i][j].tval = tval;
			lsd[i][j].pval = pval;
		}
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
	print_anova_table_part3();
	for (i = 0; i < nexp; i++) {
		v = explanatory_index[i];
		// if (strstr(var_name[v], "*")) continue; // skip interaction terms
		compute_means(v);
		print_means(v);
		compute_lsd(v);
		print_lsd(v);
	}
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

	check_atab();

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
//               R-Square     Coeff Var     Root MSE        Y Mean
//               0.938596      9.801480     1.414214     14.428571

void
print_anova_table_part2()
{
	nrow = 2;
	ncol = 4;

	check_atab();

	A(0, 0) = strdup("R-Square");
	A(0, 1) = strdup("Coeff Var");
	A(0, 2) = strdup("Root MSE");

	sprintf(buf, "%s Mean", var_name[response_index]);
	A(0, 3) = strdup(buf);

	sprintf(buf, "%0.6f", rsq);
	A(1, 0) = strdup(buf);

	sprintf(buf, "%0.6f", cv);
	A(1, 1) = strdup(buf);

	sprintf(buf, "%0.6f", sqrtmse);
	A(1, 2) = strdup(buf);

	sprintf(buf, "%0.6f", ybar);
	A(1, 3) = strdup(buf);

	fmt[0] = 0; // right justified
	fmt[1] = 0;
	fmt[2] = 0;
	fmt[3] = 0;

	print_table_and_free();
}

// Example:
//
//      Source     DF        Anova SS     Mean Square     F Value     Pr > F
//      A           1     80.04761905     80.04761905       40.02     0.0080
//      B           1     11.26666667     11.26666667        5.63     0.0982
//      A*B         1      0.40000000      0.40000000        0.20     0.6850

void
print_anova_table_part3()
{
	int i, x;
	double msq, fval, pval;

	nrow = nexp + 1;
	ncol = 6;

	check_atab();

	A(0, 0) = strdup("Source");
	A(0, 1) = strdup("DF");
	A(0, 2) = strdup("Anova SS");
	A(0, 3) = strdup("Mean Square");
	A(0, 4) = strdup("F Value");
	A(0, 5) = strdup("Pr > F");

	for (i = 0; i < nexp; i++) {

		x = explanatory_index[i];

		// Source

		A(i + 1, 0) = strdup(var_name[x]);

		// DF

		sprintf(buf, "%d", df[i]);
		A(i + 1, 1) = strdup(buf);

		if (df[i] == 0) {
			A(i + 1, 2) = strdup(".");
			A(i + 1, 3) = strdup(".");
			A(i + 1, 4) = strdup(".");
			A(i + 1, 5) = strdup(".");
			continue;
		}

		// Anova SS

		sprintf(buf, "%0.8f", ss[i]);
		A(i + 1, 2) = strdup(buf);

		// Mean Square

		msq = ss[i] / df[i];

		sprintf(buf, "%0.8f", msq);
		A(i + 1, 3) = strdup(buf);

		// F Value

		fval = msq / mse;

		sprintf(buf, "%0.2f", fval);
		A(i + 1, 4) = strdup(buf);

		// Pr > F

		pval = 1 - fdist(fval, df[i], nobs - npar);

		sprintf(buf, "%0.4f", pval);
		A(i + 1, 5) = strdup(buf);
	}

	fmt[0] = 1; // left justify
	fmt[1] = 0;
	fmt[2] = 0;
	fmt[3] = 0;
	fmt[4] = 0;
	fmt[5] = 0;

	print_table_and_free();
}

// Print the mean response for each treatment level in variable v

void
print_means(int v)
{
	int i, j, m, n;
	char *b, *s;

	emit_line_center("Mean Response");
	emit_line("");

	// count the number of interactions

	m = 0;
	s = var_name[v];
	while (*s)
		if (*s++ == '*')
			m++;

	n = level_count[v];

	nrow = n + 1;
	ncol = m + 5;

	check_atab();

	// split up interaction string

	s = var_name[v];
	for (i = 0; i < m + 1; i++) {
		b = buf;
		while (*s && *s != '*')
			*b++ = *s++;
		if (*s)
			s++;
		*b = 0;
		A(0, i) = strdup(buf);
	}

	A(0, m + 1) = strdup("N");

	s = var_name[response_index];
	sprintf(buf, "Mean %s", s);
	A(0, m + 2) = strdup(buf);

	sprintf(buf, "%g%% CI MIN", 100 * (1 - alpha));
	A(0, m + 3) = strdup(buf);

	sprintf(buf, "%g%% CI MAX", 100 * (1 - alpha));
	A(0, m + 4) = strdup(buf);

	for (i = 0; i < n; i++) {

		// Level (decompose interactions)

		s = level_name[v][i];
		for (j = 0; j < m + 1; j++) {
			b = buf;
			while (*s && *s != '*')
				*b++ = *s++;
			if (*s)
				s++;
			*b = 0;
			A(i + 1, j) = strdup(buf);
		}

		// N

		sprintf(buf, "%d", count[i]);
		A(i + 1, m + 1) = strdup(buf);

		if (count[i] < 1) {
			A(i + 1, m + 2) = strdup(".");
			A(i + 1, m + 3) = strdup(".");
			A(i + 1, m + 4) = strdup(".");
			continue;
		}

		// Mean

		sprintf(buf, "%0.6f", mean[i]);
		A(i + 1, m + 2) = strdup(buf);

		// Confidence Interval

		sprintf(buf, "%0.6f", ci_min[i]);
		A(i + 1, m + 3) = strdup(buf);

		sprintf(buf, "%0.6f", ci_max[i]);
		A(i + 1, m + 4) = strdup(buf);
	}

	for (i = 0; i < m + 1; i++)
		fmt[i] = 1; // left justify level name
	fmt[m + 1] = 0;
	fmt[m + 2] = 0;
	fmt[m + 3] = 0;
	fmt[m + 4] = 0;

	print_table_and_free();
}

void
print_lsd(int v)
{
	int i, j, k, n;
	char *s;

	emit_line_center("Least Significant Difference Test");
	emit_line("");

	n = level_count[v];

	nrow = 1 + n * (n - 1);
	ncol = 7;

	check_atab();

	s = var_name[v];
	A(0, 0) = strdup(s);
	A(0, 1) = strdup(s);

	s = var_name[response_index];
	sprintf(buf, "Delta %s", s);
	A(0, 2) = strdup(buf);

	sprintf(buf, "%g%% CI MIN", 100 * (1 - alpha));
	A(0, 3) = strdup(buf);

	sprintf(buf, "%g%% CI MAX", 100 * (1 - alpha));
	A(0, 4) = strdup(buf);

	A(0, 5) = strdup("t Value");
	A(0, 6) = strdup("Pr > |t|  ");

	k = 1;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {

			if (i == j)
				continue;

			s = level_name[v][i];
			A(k, 0) = strdup(s);

			s = level_name[v][j];
			A(k, 1) = strdup(s);

			sprintf(buf, "%0.6f", lsd[i][j].delta);
			A(k, 2) = strdup(buf);

			sprintf(buf, "%0.6f", lsd[i][j].min);
			A(k, 3) = strdup(buf);

			sprintf(buf, "%0.6f", lsd[i][j].max);
			A(k, 4) = strdup(buf);

			sprintf(buf, "%0.2f", lsd[i][j].tval);
			A(k, 5) = strdup(buf);

			if (lsd[i][j].pval > alpha)
				sprintf(buf, "%0.4f  ", lsd[i][j].pval);
			else
				sprintf(buf, "%0.4f *", lsd[i][j].pval);
			A(k, 6) = strdup(buf);

			k++;
		}
	}

	fmt[0] = 1; // left justify
	fmt[1] = 1;
	fmt[2] = 0;
	fmt[3] = 0;
	fmt[4] = 0;
	fmt[5] = 0;
	fmt[6] = 0;

	print_table_and_free();
}

void
print_ttest(int v)
{
	int dfe, i, j, k, n;
	char **a, *s;
	double d, mse, pval, se, sse, t, tval;

	emit_line_center("Two Sample t-Test");
	emit_line("");

	n = level_count[v];

	nrow = 1 + n * (n - 1);
	ncol = 7;

	check_atab();

	s = var_name[v];
	A(0, 0) = strdup(s);
	A(0, 1) = strdup(s);

	s = var_name[response_index];
	sprintf(buf, "Delta %s", s);
	A(0, 2) = strdup(buf);

	sprintf(buf, "%g%% CI MIN", 100 * (1 - alpha));
	A(0, 3) = strdup(buf);

	sprintf(buf, "%g%% CI MAX", 100 * (1 - alpha));
	A(0, 4) = strdup(buf);

	A(0, 5) = strdup("t Value");
	A(0, 6) = strdup("Pr > |t|  ");

	k = 0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {

			if (i == j)
				continue;

			k++;

			// Level

			s = level_name[v][i];
			A(k, 0) = strdup(s);

			// Level

			s = level_name[v][i];
			A(k, 1) = strdup(s);

			// Sanity check

			if (count[i] == 0 || count[j] == 0) {
				A(k, 2) = strdup(".");
				A(k, 3) = strdup(".");
				A(k, 4) = strdup(".");
				A(k, 5) = strdup(".");
				A(k, 6) = strdup(".  ");
				continue;
			}

			// Difference

			d = mean[i] - mean[j];
			sprintf(buf, "%0.6f", d);
			A(k, 2) = strdup(buf);

			// confidence interval

			if (count[i] + count[j] < 3) {
				A(k, 3) = strdup(".");
				A(k, 4) = strdup(".");
				A(k, 5) = strdup(".");
				A(k, 6) = strdup(".  ");
				continue;
			}

			sse = variance[i] * (count[i] - 1) + variance[j] * (count[j] - 1);

			dfe = count[i] + count[j] - 2;

			mse = sse / dfe;

			se = sqrt(mse * (1.0 / count[i] + 1.0 / count[j]));

			tval = d / se;

			pval = 2 * (1 - tdist(fabs(tval), dfe));

			t = qt(1 - alpha / 2, dfe) * se;

			sprintf(buf, "%0.6f", d - t);
			A(k, 3) = strdup(buf);

			sprintf(buf, "%0.6f", d + t);
			A(k, 4) = strdup(buf);

			sprintf(buf, "%0.2f", tval);
			A(k, 5) = strdup(buf);

			if (pval > alpha)
				sprintf(buf, "%0.4f  ", pval);
			else
				sprintf(buf, "%0.4f *", pval);
			A(k, 6) = strdup(buf);
		}
	}

	fmt[0] = 1; // left justify
	fmt[1] = 1;
	fmt[2] = 0;
	fmt[3] = 0;
	fmt[4] = 0;
	fmt[5] = 0;
	fmt[6] = 0;

	print_table_and_free();
}

void
check_atab()
{
	if (nrow >= MAXROW)
		stop("Try increasing MAXROW");
	if (ncol >= MAXCOL)
		stop("Try increasing MAXCOL");
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
