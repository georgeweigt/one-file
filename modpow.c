/* This is a realization of RSA encryption for any key length.

The main function encrypts a certificate signature and checks the result.

To build and run:

	gcc modpow.c
	./a.out
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// certificate values in hex ascii

char *exponent = "010001";

char *modulus =
	"00e182d4da3f6d450048539c7740586bded4444257d023d6e43b988b4a59662f"
	"cee615db4cee92b15ab5623fd66ef30336dad357270f1ee6ae6cb27f3761c1d7"
	"9963da03dae1114b2cc8f9f70f89b36fb6d296e6edf3fdb126f8a34f790dbb65"
	"d0f9449e7b7c5cd3aeafe8484dab0732b16f9e26db80cbc2487fc89fbfb48395"
	"a3bcc5bf79470f37adb11b9681bac4b07a0e5735955a64d88bd5aff11f1ad3ea"
	"5a1ed036632fa929da3412510373e812f998a18759062bde00e45619468c44f4"
	"ea0bb2dd841b2aad022c3705efe373d754551050f4e1f5894b7de6c75faafe9b"
	"c649cfe69a92f341324a9681df17c9010ca4f574e848e7fdd11040894b09a0ac"
	"8d";

char *signature =
	"6b39366e3214e777c886b038c336c7e1a38afb363f39ecde7a6b90abba72e0c8"
	"33bb4183753ded93c910d7f4fa32b5db63461e5b8f1294bb27a76b43aafcb435"
	"4ddee36352b1c4246eeb2a5baa71f0d46c23e4bafe6a15d2cce524f22156b799"
	"e9af641b3edf2433ef1a1e66b6dece964abf8b0b04e9dc5431742163e03450d7"
	"65e2e232eefff9164d7f35b058f7bd7082a57120d38c089bb41fd242362632de"
	"2cd3bc81032153e7f3580f2dff0720530aebb0cddb09808ad925f6ecb5adb9e6"
	"2149e36e83155f0778b4a7336756dfc78eaa00e9572914f1a22ee91fd9a98c92"
	"ab429d7cb83f5917aaa70cd80e60f3a4af199f6839525ac710c91852ccc96cdc";

char *result =
	"0001ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffff00302130"
	"0906052b0e03021a05000414998a88ec8c6ecfe97f035b309614c0524c8ce94c";

uint32_t * modpow(uint32_t *a, uint32_t *b, uint32_t *c);
void mshr(uint32_t *u);
uint32_t * madd(uint32_t *u, uint32_t *v);
uint32_t * msub(uint32_t *u, uint32_t *v);
uint32_t * mmul(uint32_t *u, uint32_t *v);
uint32_t * mdiv(uint32_t *u, uint32_t *v);
void mmod(uint32_t *u, uint32_t *v);
int mcmp(uint32_t *u, uint32_t *v);
uint32_t * mnew(int n);
void mfree(uint32_t *u);
uint32_t * mcopy(uint32_t *u);
void mnorm(uint32_t *u);
uint32_t * str_to_int(char *s);
char * int_to_str(uint32_t *a);

// compute (signature ** exponent) mod modulus, check the result

int
main()
{
	char *s;
	uint32_t *a, *b, *c, *y;

	a = str_to_int(signature);
	b = str_to_int(exponent);
	c = str_to_int(modulus);

	y = modpow(a, b, c);

	s = int_to_str(y);

	mfree(a);
	mfree(b);
	mfree(c);
	mfree(y);

	if (strcmp(s, result) == 0)
		printf("ok\n");
	else
		printf("not ok\n");

	free(s);
}

#define MLENGTH(p) (p)[-1]

// returns (a ** b) mod c

uint32_t *
modpow(uint32_t *a, uint32_t *b, uint32_t *c)
{
	uint32_t *t, *y;

	a = mcopy(a);
	b = mcopy(b);

	// y = 1

	y = mnew(1);
	y[0] = 1;

	for (;;) {

		if (b[0] & 1) {

			// y = (y * a) mod c

			t = mmul(y, a);
			mfree(y);
			y = t;
			mmod(y, c);
		}

		// b = b >> 1

		mshr(b);

		if (MLENGTH(b) == 1 && b[0] == 0)
			break; // b == 0

		// a = (a * a) mod c

		t = mmul(a, a);
		mfree(a);
		a = t;
		mmod(a, c);
	}

	mfree(a);
	mfree(b);

	return y;
}

void
mshr(uint32_t *u)
{
	int i;
	for (i = 0; i < MLENGTH(u) - 1; i++) {
		u[i] >>= 1;
		if (u[i + 1] & 1)
			u[i] |= 0x80000000;
	}
	u[i] >>= 1;
	mnorm(u);
}

// returns u + v

uint32_t *
madd(uint32_t *u, uint32_t *v)
{
	int i, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	if (nu > nv)
		nw = nu + 1;
	else
		nw = nv + 1;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = u[i];
	for (i = nu; i < nw; i++)
		w[i] = 0;
	t = 0;
	for (i = 0; i < nv; i++) {
		t += (uint64_t) w[i] + v[i];
		w[i] = t;
		t >>= 32;
	}
	for (i = nv; i < nw; i++) {
		t += w[i];
		w[i] = t;
		t >>= 32;
	}
	mnorm(w);
	return w;
}

// returns u - v

uint32_t *
msub(uint32_t *u, uint32_t *v)
{
	int i, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	if (nu > nv)
		nw = nu;
	else
		nw = nv;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = u[i];
	for (i = nu; i < nw; i++)
		w[i] = 0;
	t = 0;
	for (i = 0; i < nv; i++) {
		t += (uint64_t) w[i] - v[i];
		w[i] = t;
		t = (int64_t) t >> 32; // cast to extend sign
	}
	for (i = nv; i < nw; i++) {
		t += w[i];
		w[i] = t;
		t = (int64_t) t >> 32; // cast to extend sign
	}
	mnorm(w);
	return w;
}

// returns u * v

uint32_t *
mmul(uint32_t *u, uint32_t *v)
{
	int i, j, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	nw = nu + nv;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = 0;
	for (j = 0; j < nv; j++) {
		t = 0;
		for (i = 0; i < nu; i++) {
			t += (uint64_t) u[i] * v[j] + w[i + j];
			w[i + j] = t;
			t >>= 32;
		}
		w[i + j] = t;
	}
	mnorm(w);
	return w;
}

// returns floor(u / v)

uint32_t *
mdiv(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t *q, qhat, *w;
	uint64_t a, b, t;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(v) == 1 && v[0] == 0)
		return NULL; // v = 0
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	k = nu - nv;
	if (k < 0) {
		q = mnew(1);
		q[0] = 0;
		return q; // u < v, return zero
	}
	u = mcopy(u);
	q = mnew(k + 1);
	w = mnew(nv + 1);
	b = v[nv - 1];
	do {
		q[k] = 0;
		while (nu >= nv + k) {
			// estimate 32-bit partial quotient
			a = u[nu - 1];
			if (nu > nv + k)
				a = a << 32 | u[nu - 2];
			if (a < b)
				break;
			qhat = a / (b + 1);
			if (qhat == 0)
				qhat = 1;
			// w = qhat * v
			t = 0;
			for (i = 0; i < nv; i++) {
				t += (uint64_t) qhat * v[i];
				w[i] = t;
				t >>= 32;
			}
			w[nv] = t;
			// u = u - w
			t = 0;
			for (i = k; i < nu; i++) {
				t += (uint64_t) u[i] - w[i - k];
				u[i] = t;
				t = (int64_t) t >> 32; // cast to extend sign
			}
			if (t) {
				// u is negative, restore u
				t = 0;
				for (i = k; i < nu; i++) {
					t += (uint64_t) u[i] + w[i - k];
					u[i] = t;
					t >>= 32;
				}
				break;
			}
			q[k] += qhat;
			mnorm(u);
			nu = MLENGTH(u);
		}
	} while (--k >= 0);
	mnorm(q);
	mfree(u);
	mfree(w);
	return q;
}

// u = u mod v

void
mmod(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t qhat, *w;
	uint64_t a, b, t;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(v) == 1 && v[0] == 0)
		return; // v = 0
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	k = nu - nv;
	if (k < 0)
		return; // u < v
	w = mnew(nv + 1);
	b = v[nv - 1];
	do {
		while (nu >= nv + k) {
			// estimate 32-bit partial quotient
			a = u[nu - 1];
			if (nu > nv + k)
				a = a << 32 | u[nu - 2];
			if (a < b)
				break;
			qhat = a / (b + 1);
			if (qhat == 0)
				qhat = 1;
			// w = qhat * v
			t = 0;
			for (i = 0; i < nv; i++) {
				t += (uint64_t) qhat * v[i];
				w[i] = t;
				t >>= 32;
			}
			w[nv] = t;
			// u = u - w
			t = 0;
			for (i = k; i < nu; i++) {
				t += (uint64_t) u[i] - w[i - k];
				u[i] = t;
				t = (int64_t) t >> 32; // cast to extend sign
			}
			if (t) {
				// u is negative, restore u
				t = 0;
				for (i = k; i < nu; i++) {
					t += (uint64_t) u[i] + w[i - k];
					u[i] = t;
					t >>= 32;
				}
				break;
			}
			mnorm(u);
			nu = MLENGTH(u);
		}
	} while (--k >= 0);
	mfree(w);
}

// compare u and v

int
mcmp(uint32_t *u, uint32_t *v)
{
	int i;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(u) < MLENGTH(v))
		return -1;
	if (MLENGTH(u) > MLENGTH(v))
		return 1;
	for (i = MLENGTH(u) - 1; i >= 0; i--) {
		if (u[i] < v[i])
			return -1;
		if (u[i] > v[i])
			return 1;
	}
	return 0; // u = v
}

uint32_t *
mnew(int n)
{
	uint32_t *u;
	u = (uint32_t *) malloc((n + 1) * sizeof (uint32_t));
	if (u == NULL) {
		printf("malloc kaput\n");
		exit(1);
	}
	*u = n;
	return u + 1;
}

void
mfree(uint32_t *u)
{
	free(u - 1);
}

uint32_t *
mcopy(uint32_t *u)
{
	int i;
	uint32_t *v;
	v = mnew(MLENGTH(u));
	for (i = 0; i < MLENGTH(u); i++)
		v[i] = u[i];
	return v;
}

// remove leading zeroes

void
mnorm(uint32_t *u)
{
	while (MLENGTH(u) > 1 && u[MLENGTH(u) - 1] == 0)
		MLENGTH(u)--;
}

uint32_t *
str_to_int(char *s)
{
	int d, i, len, n;
	uint32_t *u;
	len = strlen(s);
	n = (len + 7) / 8; // convert len to number of uint32_t ints
	u = mnew(n);
	for (i = 0; i < n; i++)
		u[i] = 0;
	for (i = 0; i < len; i++) {
		d = s[len - i - 1];
		if ('0' <= d && d <= '9')
			d -= '0';
		else if ('A' <= d && d <= 'F')
			d -= 'A' - 10;
		else if ('a' <= d && d <= 'f')
			d -= 'a' - 10;
		else {
			mfree(u);
			return NULL;
		}
		u[i / 8] |= d << (4 * (i % 8));
	}
	mnorm(u);
	return u;
}

char *
int_to_str(uint32_t *u)
{
	int i;
	char *s;
	s = malloc(8 * MLENGTH(u) + 1);
	if (s == NULL)
		return NULL;
	for (i = 0; i < MLENGTH(u); i++)
		sprintf(s + 8 * i, "%08x", u[MLENGTH(u) - i - 1]);
	return s;
}
