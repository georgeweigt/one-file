/* This is a realization of 'Mathematical routines for the NIST prime elliptic curves'

The main function runs a self test.

To build and run:

	gcc ec-nist.c
	./a.out
*/

#include <stdio.h> // printf
#include <stdlib.h> // abs, exit, free, malloc
#include <string.h> // strlen
#include <stdint.h>

#define len(p) (p)[-1]

int malloc_count;

struct point {
	uint32_t *x, *y, *z;
};

uint32_t * ec_modinv(uint32_t *c, uint32_t *p);
void ec_projectify(struct point *S);
void ec_affinify(struct point *S, uint32_t *p);
void ec_double(struct point *R, struct point *S, uint32_t *p);
void ec_add_xyz(struct point *R, struct point *S, struct point *T, uint32_t *p);
void ec_full_add(struct point *R, struct point *S, struct point *T, uint32_t *p);
void ec_full_sub(struct point *R, struct point *S, struct point *T, uint32_t *p);
void ec_mult(struct point *R, uint32_t *d, struct point *S, uint32_t *p);
int ec_get_msbit_index(uint32_t *u);
int ec_get_bit(uint32_t *u, int k);
int ec_F(int t);
void ec_twin_mult(struct point *R, uint32_t *d0, struct point *S, uint32_t *d1, struct point *T, uint32_t *p);
void ec_free_xyz(struct point *u);
uint32_t * ec_add(uint32_t *u, uint32_t *v);
uint32_t * ec_sub(uint32_t *u, uint32_t *v);
uint32_t * ec_mul(uint32_t *u, uint32_t *v);
uint32_t * ec_div(uint32_t *u, uint32_t *v);
void ec_mod(uint32_t *u, uint32_t *v);
uint32_t * ec_pow(uint32_t *u, uint32_t *v);
void ec_shr(uint32_t *u);
int ec_cmp(uint32_t *u, uint32_t *v);
int ec_equal(uint32_t *u, uint32_t v);
uint32_t * ec_int(int k);
uint32_t * ec_new(int n);
void ec_free(uint32_t *u);
uint32_t * ec_dup(uint32_t *u);
void ec_norm(uint32_t *u);
uint32_t * ec_hexstr_to_bignum(char *s);
void ec_test();
void ec_test_full_add();
void ec_test_full_sub();
void ec_test_double();
void ec_test_mult();
void ec_test_twin_mult();

int
main()
{
	ec_test();
	printf("malloc_count %d\n", malloc_count); // should be zero (no memory leaks)
}

// returns 1/c mod p

uint32_t *
ec_modinv(uint32_t *c, uint32_t *p)
{
	uint32_t *k, *r, *u, *v, *t, *x1, *x2;
	u = ec_dup(c);
	v = ec_dup(p);
	x1 = ec_int(1);
	x2 = ec_int(0);
	while (!ec_equal(u, 1) && !ec_equal(v, 1)) {
		while ((u[0] & 1) == 0) {
			ec_shr(u);
			if (x1[0] & 1) {
				t = ec_add(x1, p);
				ec_free(x1);
				x1 = t;
			}
			ec_shr(x1);
		}
		while ((v[0] & 1) == 0) {
			ec_shr(v);
			if (x2[0] & 1) {
				t = ec_add(x2, p);
				ec_free(x2);
				x2 = t;
			}
			ec_shr(x2);
		}
		if (ec_cmp(u, v) >= 0) {
			t = ec_sub(u, v);
			ec_free(u);
			u = t;
			// x1 = x1 - x2
			k = ec_sub(p, x2);
			t = ec_add(x1, k);
			ec_free(x1);
			x1 = t;
			ec_mod(x1, p);
			ec_free(k);
		} else {
			t = ec_sub(v, u);
			ec_free(v);
			v = t;
			// x2 = x2 - x1
			k = ec_sub(p, x1);
			t = ec_add(x2, k);
			ec_free(x2);
			x2 = t;
			ec_mod(x2, p);
			ec_free(k);
		}
	}
	if (ec_equal(u, 1)) {
		r = x1;
		ec_free(x2);
	} else {
		r = x2;
		ec_free(x1);
	}
	ec_free(u);
	ec_free(v);
	return r;
}

void
ec_projectify(struct point *S)
{
	ec_free(S->z);
	S->z = ec_int(1);
}

void
ec_affinify(struct point *S, uint32_t *p)
{
	uint32_t *lambda, *lambda2, *lambda3, *x, *y;

	if (ec_equal(S->z, 0)) {
		printf("cannot affinify\n");
		return;
	}

	lambda = ec_modinv(S->z, p);

	lambda2 = ec_mul(lambda, lambda);
	ec_mod(lambda2, p);

	lambda3 = ec_mul(lambda2, lambda);
	ec_mod(lambda3, p);

	x = ec_mul(lambda2, S->x);
	ec_mod(x, p);

	y = ec_mul(lambda3, S->y);
	ec_mod(y, p);

	ec_free_xyz(S);

	S->x = x;
	S->y = y;

	ec_free(lambda);
	ec_free(lambda2);
	ec_free(lambda3);
}

void
ec_double(struct point *R, struct point *S, uint32_t *p)
{
	uint32_t *k, *t, *t1, *t2, *t3, *t4, *t5;

	// take care to handle the case when R and S are the same pointer

	t1 = ec_dup(S->x);
	t2 = ec_dup(S->y);
	t3 = ec_dup(S->z);

	ec_free_xyz(R);

	if (ec_equal(t3, 0)) {
		R->x = ec_int(1);
		R->y = ec_int(1);
		R->z = ec_int(0);
		ec_free(t1);
		ec_free(t2);
		ec_free(t3);
		return;
	}

	// 7: t4 = t3 * t3

	t4 = ec_mul(t3, t3);
	ec_mod(t4, p);

	// 8: t5 = t1 - t4

	t = ec_sub(p, t4);
	t5 = ec_add(t1, t);
	ec_free(t);
	ec_mod(t5, p);

	// 9: t4 = t1 + t4

	t = ec_add(t1, t4);
	ec_free(t4);
	t4 = t;
	ec_mod(t4, p);

	// 10: t5 = t4 * t5

	t = ec_mul(t4, t5);
	ec_free(t5);
	t5 = t;
	ec_mod(t5, p);

	// 11: t4 = 3 * t5

	k = ec_int(3);
	ec_free(t4);
	t4 = ec_mul(k, t5);
	ec_free(k);
	ec_mod(t4, p);

	// 12: t3 = t3 * t2

	t = ec_mul(t3, t2);
	ec_free(t3);
	t3 = t;
	ec_mod(t3, p);

	// 13: t3 = 2 * t3

	t = ec_add(t3, t3);
	ec_free(t3);
	t3 = t;
	ec_mod(t3, p);

	// 14: t2 = t2 * t2

	t = ec_mul(t2, t2);
	ec_free(t2);
	t2 = t;
	ec_mod(t2, p);

	// 15: t5 = t1 * t2

	t = ec_mul(t1, t2);
	ec_free(t5);
	t5 = t;
	ec_mod(t5, p);

	// 16: t5 = 4 * t5

	k = ec_int(4);
	t = ec_mul(k, t5);
	ec_free(t5);
	t5 = t;
	ec_free(k);
	ec_mod(t5, p);

	// 17: t1 = t4 * t4

	ec_free(t1);
	t1 = ec_mul(t4, t4);
	ec_mod(t1, p);

	// 18: t1 = t1 - 2 * t5

	k = ec_sub(p, t5);
	t = ec_add(t1, k);
	ec_free(t1);
	t1 = ec_add(t, k);
	ec_free(k);
	ec_free(t);
	ec_mod(t1, p);

	// 19: t2 = t2 * t2

	t = ec_mul(t2, t2);
	ec_free(t2);
	t2 = t;
	ec_mod(t2, p);

	// 20: t2 = 8 * t2

	k = ec_int(8);
	t = ec_mul(k, t2);
	ec_free(t2);
	t2 = t;
	ec_free(k);
	ec_mod(t2, p);

	// 21: t5 = t5 - t1

	k = ec_sub(p, t1);
	t = ec_add(t5, k);
	ec_free(t5);
	t5 = t;
	ec_free(k);
	ec_mod(t5, p);

	// 22: t5 = t4 * t5

	t = ec_mul(t4, t5);
	ec_free(t5);
	t5 = t;
	ec_mod(t5, p);

	// 23: t2 = t5 - t2

	t = ec_sub(p, t2);
	ec_free(t2);
	t2 = ec_add(t5, t);
	ec_free(t);
	ec_mod(t2, p);

	R->x = t1;
	R->y = t2;
	R->z = t3;

	ec_free(t4);
	ec_free(t5);
}

void
ec_add_xyz(struct point *R, struct point *S, struct point *T, uint32_t *p)
{
	uint32_t *k, *t, *t1, *t2, *t3, *t4, *t5, *t6, *t7;

	t1 = ec_dup(S->x);
	t2 = ec_dup(S->y);
	t3 = ec_dup(S->z);

	t4 = ec_dup(T->x);
	t5 = ec_dup(T->y);
	t6 = ec_dup(T->z);

	ec_free_xyz(R);

	if (!ec_equal(t6, 1)) {

		// 4: t7 = t6 * t6

		t7 = ec_mul(t6, t6);
		ec_mod(t7, p);

		// 5: t1 = t1 * t7

		t = ec_mul(t1, t7);
		ec_free(t1);
		t1 = t;
		ec_mod(t1, p);

		// 6: t7 = t6 * t7

		t = ec_mul(t6, t7);
		ec_free(t7);
		t7 = t;
		ec_mod(t7, p);

		// 7: t2 = t2 * t7

		t = ec_mul(t2, t7);
		ec_free(t2);
		t2 = t;
		ec_mod(t2, p);

		ec_free(t7);
	}

	// 9: t7 = t3 * t3

	t7 = ec_mul(t3, t3);
	ec_mod(t7, p);

	// 10: t4 = t4 * t7

	t = ec_mul(t4, t7);
	ec_free(t4);
	t4 = t;
	ec_mod(t4, p);

	// 11: t7 = t3 * t7

	t = ec_mul(t3, t7);
	ec_free(t7);
	t7 = t;
	ec_mod(t7, p);

	// 12: t5 = t5 * t7

	t = ec_mul(t5, t7);
	ec_free(t5);
	t5 = t;
	ec_mod(t5, p);

	// 13: t4 = t1 - t4

	t = ec_sub(p, t4);
	ec_free(t4);
	t4 = ec_add(t1, t);
	ec_free(t);
	ec_mod(t4, p);

	// 14: t5 = t2 - t5

	t = ec_sub(p, t5);
	ec_free(t5);
	t5 = ec_add(t2, t);
	ec_free(t);
	ec_mod(t5, p);

	if (ec_equal(t4, 0)) {
		if (ec_equal(t5, 0)) {
			R->x = ec_int(0);
			R->y = ec_int(0);
			R->z = ec_int(0);
		} else {
			R->x = ec_int(1);
			R->y = ec_int(1);
			R->z = ec_int(0);
		}
		ec_free(t1);
		ec_free(t2);
		ec_free(t3);
		ec_free(t4);
		ec_free(t5);
		ec_free(t6);
		ec_free(t7);
		return;
	}

	// 22: t1 = 2 * t1 - t4

	t = ec_add(t1, t1);
	ec_free(t1);
	t1 = t;
	k = ec_sub(p, t4);
	t = ec_add(t1, k);
	ec_free(t1);
	t1 = t;
	ec_free(k);
	ec_mod(t1, p);

	// 23: t2 = 2 * t2 - t5

	t = ec_add(t2, t2);
	ec_free(t2);
	t2 = t;
	k = ec_sub(p, t5);
	t = ec_add(t2, k);
	ec_free(t2);
	t2 = t;
	ec_free(k);
	ec_mod(t2, p);

	if (!ec_equal(t6, 1)) {

		// 25: t3 = t3 * t6

		t = ec_mul(t3, t6);
		ec_free(t3);
		t3 = t;
		ec_mod(t3, p);
	}

	// 27: t3 = t3 * t4

	t = ec_mul(t3, t4);
	ec_free(t3);
	t3 = t;
	ec_mod(t3, p);

	// 28: t7 = t4 * t4

	t = ec_mul(t4, t4);
	ec_free(t7);
	t7 = t;
	ec_mod(t7, p);

	// 29: t4 = t4 * t7

	t = ec_mul(t4, t7);
	ec_free(t4);
	t4 = t;
	ec_mod(t4, p);

	// 30: t7 = t1 * t7

	t = ec_mul(t1, t7);
	ec_free(t7);
	t7 = t;
	ec_mod(t7, p);

	// 31: t1 = t5 * t5

	ec_free(t1);
	t1 = ec_mul(t5, t5);
	ec_mod(t1, p);

	// 32: t1 = t1 - t7

	k = ec_sub(p, t7);
	t = ec_add(t1, k);
	ec_free(t1);
	t1 = t;
	ec_free(k);
	ec_mod(t1, p);

	// 33: t7 = t7 - 2 * t1

	k = ec_sub(p, t1);
	t = ec_add(t7, k);
	ec_free(t7);
	t7 = ec_add(t, k);
	ec_free(k);
	ec_free(t);
	ec_mod(t7, p);

	// 34: t5 = t5 * t7

	t = ec_mul(t5, t7);
	ec_free(t5);
	t5 = t;
	ec_mod(t5, p);

	// 35: t4 = t2 * t4

	t = ec_mul(t2, t4);
	ec_free(t4);
	t4 = t;
	ec_mod(t4, p);

	// 36: t2 = t5 - t4

	t = ec_sub(p, t4);
	ec_free(t2);
	t2 = ec_add(t5, t);
	ec_free(t);
	ec_mod(t2, p);

	// 37: t2 = t2 / 2

	if (t2[0] & 1) {
		t = ec_add(t2, p);
		ec_free(t2);
		t2 = t;
	}
	ec_shr(t2);

	R->x = t1;
	R->y = t2;
	R->z = t3;

	ec_free(t4);
	ec_free(t5);
	ec_free(t6);
	ec_free(t7);
}

void
ec_full_add(struct point *R, struct point *S, struct point *T, uint32_t *p)
{
	uint32_t *x, *y, *z;
	struct point U;

	if (ec_equal(S->z, 0)) {
		x = ec_dup(T->x);
		y = ec_dup(T->y);
		z = ec_dup(T->z);
		ec_free_xyz(R);
		R->x = x;
		R->y = y;
		R->z = z;
		return;
	}

	if (ec_equal(T->z, 0)) {
		x = ec_dup(S->x);
		y = ec_dup(S->y);
		z = ec_dup(S->z);
		ec_free_xyz(R);
		R->x = x;
		R->y = y;
		R->z = z;
		return;
	}

	U.x = NULL;
	U.y = NULL;
	U.z = NULL;

	ec_add_xyz(&U, S, T, p);

	if (ec_equal(U.x, 0) && ec_equal(U.y, 0) && ec_equal(U.z, 0)) {
		ec_free_xyz(&U);
		ec_double(&U, S, p);
	}

	ec_free_xyz(R);

	R->x = U.x;
	R->y = U.y;
	R->z = U.z;
}

void
ec_full_sub(struct point *R, struct point *S, struct point *T, uint32_t *p)
{
	struct point U;

	U.x = ec_dup(T->x);
	U.y = ec_sub(p, T->y);
	U.z = ec_dup(T->z);

	ec_full_add(R, S, &U, p);

	ec_free_xyz(&U);
}

void
ec_mult(struct point *R, uint32_t *d, struct point *S, uint32_t *p)
{
	int h, i, k, l;
	uint32_t *t, *u, *x, *y, *z;
	struct point U;

	if (ec_equal(d, 0)) {
		ec_free_xyz(R);
		R->x = ec_int(1);
		R->y = ec_int(1);
		R->z = ec_int(0);
		return;
	}

	if (ec_equal(d, 1)) {
		x = ec_dup(S->x);
		y = ec_dup(S->y);
		z = ec_dup(S->z);
		ec_free_xyz(R);
		R->x = x;
		R->y = y;
		R->z = z;
		return;
	}

	if (ec_equal(S->z, 0)) {
		ec_free_xyz(R);
		R->x = ec_int(1);
		R->y = ec_int(1);
		R->z = ec_int(0);
		return;
	}

	if (!ec_equal(S->z, 1)) {
		ec_affinify(S, p);
		ec_projectify(S);
	}

	x = ec_dup(S->x);
	y = ec_dup(S->y);
	z = ec_dup(S->z);

	ec_free_xyz(R);

	R->x = x;
	R->y = y;
	R->z = z;

	u = ec_int(3);
	t = ec_mul(u, d);
	ec_free(u);

	l = ec_get_msbit_index(t);

	for (i = l - 1; i > 0; i--) {

		U.x = NULL;
		U.y = NULL;
		U.z = NULL;

		ec_double(R, R, p);

		h = ec_get_bit(t, i);
		k = ec_get_bit(d, i);

		if (h == 1 && k == 0)
			ec_full_add(&U, R, S, p);

		if (h == 0 && k == 1)
			ec_full_sub(&U, R, S, p);

		if (h != k) {
			ec_free_xyz(R);
			R->x = U.x;
			R->y = U.y;
			R->z = U.z;
		}
	}

	ec_free(t);
}

int
ec_get_msbit_index(uint32_t *u)
{
	int k, n;
	uint32_t m;
	m = 0x80000000;
	n = len(u);
	k = 32 * n - 1;
	while (m > 1) {
		if (u[n - 1] & m)
			break;
		m >>= 1;
		k--;
	}
	return k;
}

int
ec_get_bit(uint32_t *u, int k)
{
	int j;
	uint32_t m;
	if (k < 0)
		return 0;
	j = k / 32;
	if (j >= len(u))
		return 0;
	m = 1 << (k % 32);
	if (u[j] & m)
		return 1;
	else
		return 0;
}

int
ec_F(int t)
{
	if (18 <= t && t < 22)
		return 9;

	if (14 <= t && t < 18)
		return 10;

	if (22 <= t && t < 24)
		return 11;

	if (4 <= t && t < 12)
		return 14;

	return 12;
}

// R cannot point to S or T

void
ec_twin_mult(struct point *R, uint32_t *d0, struct point *S, uint32_t *d1, struct point *T, uint32_t *p)
{
	int c[2][6], h[2], i, k, m, m0, m1, u[2];
	struct point SpT, SmT;

	if (R == S || R == T) {
		printf("arg error\n");
		return;
	}

	SpT.x = NULL;
	SpT.y = NULL;
	SpT.z = NULL;

	SmT.x = NULL;
	SmT.y = NULL;
	SmT.z = NULL;

	ec_full_add(&SpT, S, T, p);
	ec_full_sub(&SmT, S, T, p);

	m0 = ec_get_msbit_index(d0) + 1;
	m1 = ec_get_msbit_index(d1) + 1;

	if (m0 > m1)
		m = m0;
	else
		m = m1;

	c[0][0] = 0;
	c[0][1] = 0;
	c[0][2] = ec_get_bit(d0, m - 1);
	c[0][3] = ec_get_bit(d0, m - 2);
	c[0][4] = ec_get_bit(d0, m - 3);
	c[0][5] = ec_get_bit(d0, m - 4);

	c[1][0] = 0;
	c[1][1] = 0;
	c[1][2] = ec_get_bit(d1, m - 1);
	c[1][3] = ec_get_bit(d1, m - 2);
	c[1][4] = ec_get_bit(d1, m - 3);
	c[1][5] = ec_get_bit(d1, m - 4);

	R->x = ec_int(1);
	R->y = ec_int(1);
	R->z = ec_int(0);

	for (k = m; k > -1; k--) {

		for (i = 0; i < 2; i++) {
			h[i] = 16 * c[i][1] + 8 * c[i][2] + 4 * c[i][3] + 2 * c[i][4] + c[i][5];
			if (c[i][0] == 1)
				h[i] = 31 - h[i];
		}

		for (i = 0; i < 2; i++) {
			if (h[i] < ec_F(h[1 - i]))
				u[i] = 0;
			else {
				if (c[i][0] & 1)
					u[i] = -1;
				else
					u[i] = 1;
			}
		}

		c[0][0] = abs(u[0]) ^ c[0][1];
		c[0][1] = c[0][2];
		c[0][2] = c[0][3];
		c[0][3] = c[0][4];
		c[0][4] = c[0][5];
		c[0][5] = ec_get_bit(d0, k - 5);

		c[1][0] = abs(u[1]) ^ c[1][1];
		c[1][1] = c[1][2];
		c[1][2] = c[1][3];
		c[1][3] = c[1][4];
		c[1][4] = c[1][5];
		c[1][5] = ec_get_bit(d1, k - 5);

		ec_double(R, R, p);

		if (u[0] == -1 && u[1] == -1)
			ec_full_sub(R, R, &SpT, p);

		if (u[0] == -1 && u[1] == 0)
			ec_full_sub(R, R, S, p);

		if (u[0] == -1 && u[1] == 1)
			ec_full_sub(R, R, &SmT, p);

		if (u[0] == 0 && u[1] == -1)
			ec_full_sub(R, R, T, p);

		if (u[0] == 0 && u[1] == 1)
			ec_full_add(R, R, T, p);

		if (u[0] == 1 && u[1] == -1)
			ec_full_add(R, R, &SmT, p);

		if (u[0] == 1 && u[1] == 0)
			ec_full_add(R, R, S, p);

		if (u[0] == 1 && u[1] == 1)
			ec_full_add(R, R, &SpT, p);
	}

	ec_free_xyz(&SpT);
	ec_free_xyz(&SmT);
}

void
ec_free_xyz(struct point *u)
{
	ec_free(u->x);
	ec_free(u->y);
	ec_free(u->z);
	u->x = NULL;
	u->y = NULL;
	u->z = NULL;
}

// returns u + v

uint32_t *
ec_add(uint32_t *u, uint32_t *v)
{
	int i, n, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = len(u);
	nv = len(v);
	if (nu > nv)
		nw = nu + 1;
	else
		nw = nv + 1;
	w = ec_new(nw);
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
	ec_norm(w);
	return w;
}

// returns u - v

uint32_t *
ec_sub(uint32_t *u, uint32_t *v)
{
	int i, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = len(u);
	nv = len(v);
	if (nu > nv)
		nw = nu;
	else
		nw = nv;
	w = ec_new(nw);
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
	ec_norm(w);
	return w;
}

// returns u * v

uint32_t *
ec_mul(uint32_t *u, uint32_t *v)
{
	int i, j, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = len(u);
	nv = len(v);
	nw = nu + nv;
	w = ec_new(nw);
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
	ec_norm(w);
	return w;
}

// returns floor(u / v)

uint32_t *
ec_div(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t *q, qhat, *w;
	uint64_t a, b, t;
	ec_norm(u);
	ec_norm(v);
	if (len(v) == 1 && v[0] == 0)
		return NULL; // v = 0
	nu = len(u);
	nv = len(v);
	k = nu - nv;
	if (k < 0) {
		q = ec_new(1);
		q[0] = 0;
		return q; // u < v, return zero
	}
	u = ec_dup(u);
	q = ec_new(k + 1);
	w = ec_new(nv + 1);
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
			ec_norm(u);
			nu = len(u);
		}
	} while (--k >= 0);
	ec_norm(q);
	ec_free(u);
	ec_free(w);
	return q;
}

// u = u mod v

void
ec_mod(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t qhat, *w;
	uint64_t a, b, t;
	ec_norm(u);
	ec_norm(v);
	if (len(v) == 1 && v[0] == 0)
		return; // v = 0
	nu = len(u);
	nv = len(v);
	k = nu - nv;
	if (k < 0)
		return; // u < v
	w = ec_new(nv + 1);
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
			ec_norm(u);
			nu = len(u);
		}
	} while (--k >= 0);
	ec_free(w);
}

// returns u ** v

uint32_t *
ec_pow(uint32_t *u, uint32_t *v)
{
	uint32_t *t, *w;
	u = ec_dup(u);
	v = ec_dup(v);
	// w = 1
	w = ec_new(1);
	w[0] = 1;
	for (;;) {
		if (v[0] & 1) {
			// w = w * u
			t = ec_mul(w, u);
			ec_free(w);
			w = t;
		}
		// v = v >> 1
		ec_shr(v);
		// v = 0?
		if (len(v) == 1 && v[0] == 0)
			break;
		// u = u * u
		t = ec_mul(u, u);
		ec_free(u);
		u = t;
	}
	ec_free(u);
	ec_free(v);
	return w;
}

// u = u >> 1

void
ec_shr(uint32_t *u)
{
	int i;
	for (i = 0; i < len(u) - 1; i++) {
		u[i] >>= 1;
		if (u[i + 1] & 1)
			u[i] |= 0x80000000;
	}
	u[i] >>= 1;
	ec_norm(u);
}

// compare u and v

int
ec_cmp(uint32_t *u, uint32_t *v)
{
	int i;
	ec_norm(u);
	ec_norm(v);
	if (len(u) < len(v))
		return -1;
	if (len(u) > len(v))
		return 1;
	for (i = len(u) - 1; i >= 0; i--) {
		if (u[i] < v[i])
			return -1;
		if (u[i] > v[i])
			return 1;
	}
	return 0; // u = v
}

int
ec_equal(uint32_t *u, uint32_t v)
{
	if (len(u) == 1 && u[0] == v)
		return 1;
	else
		return 0;
}

uint32_t *
ec_int(int k)
{
	uint32_t *u;
	u = ec_new(1);
	u[0] = k;
	return u;
}

uint32_t *
ec_new(int n)
{
	uint32_t *u;
	u = (uint32_t *) malloc((n + 1) * sizeof (uint32_t));
	if (u == NULL) {
		printf("malloc error\n");
		exit(1);
	}
	malloc_count++;
	u[0] = n;
	return u + 1;
}

void
ec_free(uint32_t *u)
{
	if (u) {
		free(u - 1);
		malloc_count--;
	}
}

uint32_t *
ec_dup(uint32_t *u)
{
	int i;
	uint32_t *v;
	v = ec_new(len(u));
	for (i = 0; i < len(u); i++)
		v[i] = u[i];
	return v;
}

// remove leading zeroes

void
ec_norm(uint32_t *u)
{
	while (len(u) > 1 && u[len(u) - 1] == 0)
		len(u)--;
}

uint32_t *
ec_hexstr_to_bignum(char *s)
{
	int d, i, len, n;
	uint32_t *u;
	len = strlen(s);
	n = (len + 7) / 8; // convert len to number of uint32_t ints
	u = ec_new(n);
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
			ec_free(u);
			return NULL;
		}
		u[i / 8] |= d << (4 * (i % 8));
	}
	ec_norm(u);
	return u;
}

void
ec_test()
{
	ec_test_full_add();
	ec_test_full_sub();
	ec_test_double();
	ec_test_mult();
	ec_test_twin_mult();
}

void
ec_test_full_add()
{
	uint32_t *p, *x, *y;
	struct point R, S, T;

	char *str_p384 =
		"ffffffffffffffffffffffffffffffffffffffffffffffff"
		"fffffffffffffffeffffffff0000000000000000ffffffff";

	char *str_xs =
		"fba203b81bbd23f2b3be971cc23997e1ae4d89e69cb6f923"
		"85dda82768ada415ebab4167459da98e62b1332d1e73cb0e";

	char *str_ys =
		"5ffedbaefdeba603e7923e06cdb5d0c65b22301429293376"
		"d5c6944e3fa6259f162b4788de6987fd59aed5e4b5285e45";

	char *str_xt =
		"aacc05202e7fda6fc73d82f0a66220527da8117ee8f8330e"
		"ad7d20ee6f255f582d8bd38c5a7f2b40bcdb68ba13d81051";

	char *str_yt =
		"84009a263fefba7c2c57cffa5db3634d286131afc0fca8d2"
		"5afa22a7b5dce0d9470da89233cee178592f49b6fecb5092";

	char *str_xr =
		"12dc5ce7acdfc5844d939f40b4df012e68f865b89c3213ba"
		"97090a247a2fc009075cf471cd2e85c489979b65ee0b5eed";

	char *str_yr =
		"167312e58fe0c0afa248f2854e3cddcb557f983b3189b67f"
		"21eee01341e7e9fe67f6ee81b36988efa406945c8804a4b0";

	p = ec_hexstr_to_bignum(str_p384);

	S.x = ec_hexstr_to_bignum(str_xs);
	S.y = ec_hexstr_to_bignum(str_ys);
	S.z = ec_int(1);

	T.x = ec_hexstr_to_bignum(str_xt);
	T.y = ec_hexstr_to_bignum(str_yt);
	T.z = ec_int(1);

	R.x = NULL;
	R.y = NULL;
	R.z = NULL;

	ec_full_add(&R, &S, &T, p);

	ec_affinify(&R, p);

	x = ec_hexstr_to_bignum(str_xr);
	y = ec_hexstr_to_bignum(str_yr);

	if (ec_cmp(R.x, x) == 0 && ec_cmp(R.y, y) == 0)
		printf("ok: ec_full_add\n");
	else
		printf("err: ec_full_add\n");

	ec_free(p);
	ec_free(x);
	ec_free(y);

	ec_free_xyz(&R);
	ec_free_xyz(&S);
	ec_free_xyz(&T);
}

void
ec_test_full_sub()
{
	uint32_t *p, *x, *y;
	struct point R, S, T;

	char *str_p384 =
		"ffffffffffffffffffffffffffffffffffffffffffffffff"
		"fffffffffffffffeffffffff0000000000000000ffffffff";

	char *str_xs =
		"fba203b81bbd23f2b3be971cc23997e1ae4d89e69cb6f923"
		"85dda82768ada415ebab4167459da98e62b1332d1e73cb0e";

	char *str_ys =
		"5ffedbaefdeba603e7923e06cdb5d0c65b22301429293376"
		"d5c6944e3fa6259f162b4788de6987fd59aed5e4b5285e45";

	char *str_xt =
		"aacc05202e7fda6fc73d82f0a66220527da8117ee8f8330e"
		"ad7d20ee6f255f582d8bd38c5a7f2b40bcdb68ba13d81051";

	char *str_yt =
		"84009a263fefba7c2c57cffa5db3634d286131afc0fca8d2"
		"5afa22a7b5dce0d9470da89233cee178592f49b6fecb5092";

	char *str_xr =
		"6afdaf8da8b11c984cf177e551cee542cda4ac2f25cd522d"
		"0cd710f88059c6565aef78f6b5ed6cc05a6666def2a2fb59";

	char *str_yr =
		"7bed0e158ae8cc70e847a60347ca1548c348decc6309f48b"
		"59bd5afc9a9b804e7f7876178cb5a7eb4f6940a9c73e8e5e";

	p = ec_hexstr_to_bignum(str_p384);

	S.x = ec_hexstr_to_bignum(str_xs);
	S.y = ec_hexstr_to_bignum(str_ys);
	S.z = ec_int(1);

	T.x = ec_hexstr_to_bignum(str_xt);
	T.y = ec_hexstr_to_bignum(str_yt);
	T.z = ec_int(1);

	R.x = NULL;
	R.y = NULL;
	R.z = NULL;

	ec_full_sub(&R, &S, &T, p);

	ec_affinify(&R, p);

	x = ec_hexstr_to_bignum(str_xr);
	y = ec_hexstr_to_bignum(str_yr);

	if (ec_cmp(R.x, x) == 0 && ec_cmp(R.y, y) == 0)
		printf("ok: ec_full_sub\n");
	else
		printf("err: ec_full_sub\n");

	ec_free(p);
	ec_free(x);
	ec_free(y);

	ec_free_xyz(&R);
	ec_free_xyz(&S);
	ec_free_xyz(&T);
}

void
ec_test_double()
{
	uint32_t *p, *x, *y;
	struct point R, S;

	char *str_p384 =
		"ffffffffffffffffffffffffffffffffffffffffffffffff"
		"fffffffffffffffeffffffff0000000000000000ffffffff";

	char *str_xs =
		"fba203b81bbd23f2b3be971cc23997e1ae4d89e69cb6f923"
		"85dda82768ada415ebab4167459da98e62b1332d1e73cb0e";

	char *str_ys =
		"5ffedbaefdeba603e7923e06cdb5d0c65b22301429293376"
		"d5c6944e3fa6259f162b4788de6987fd59aed5e4b5285e45";

	char *str_xr =
		"2a2111b1e0aa8b2fc5a1975516bc4d58017ff96b25e1bdff"
		"3c229d5fac3bacc319dcbec29f9478f42dee597b4641504c";

	char *str_yr =
		"fa2e3d9dc84db8954ce8085ef28d7184fddfd1344b4d4797"
		"343af9b5f9d837520b450f726443e4114bd4e5bdb2f65ddd";

	p = ec_hexstr_to_bignum(str_p384);

	S.x = ec_hexstr_to_bignum(str_xs);
	S.y = ec_hexstr_to_bignum(str_ys);
	S.z = ec_int(1);

	R.x = NULL;
	R.y = NULL;
	R.z = NULL;

	ec_double(&R, &S, p);

	ec_affinify(&R, p);

	x = ec_hexstr_to_bignum(str_xr);
	y = ec_hexstr_to_bignum(str_yr);

	if (ec_cmp(R.x, x) == 0 && ec_cmp(R.y, y) == 0)
		printf("ok: ec_double\n");
	else
		printf("err: ec_double\n");

	ec_free(p);
	ec_free(x);
	ec_free(y);

	ec_free_xyz(&R);
	ec_free_xyz(&S);
}

void
ec_test_mult()
{
	uint32_t *d, *p, *x, *y;
	struct point R, S;

	char *str_p384 =
		"ffffffffffffffffffffffffffffffffffffffffffffffff"
		"fffffffffffffffeffffffff0000000000000000ffffffff";

	char *str_xs =
		"fba203b81bbd23f2b3be971cc23997e1ae4d89e69cb6f923"
		"85dda82768ada415ebab4167459da98e62b1332d1e73cb0e";

	char *str_ys =
		"5ffedbaefdeba603e7923e06cdb5d0c65b22301429293376"
		"d5c6944e3fa6259f162b4788de6987fd59aed5e4b5285e45";

	char *str_d =
		"a4ebcae5a665983493ab3e626085a24c104311a761b5a8fd"
		"ac052ed1f111a5c44f76f45659d2d111a61b5fdd97583480";

	char *str_xr =
		"e4f77e7ffeb7f0958910e3a680d677a477191df166160ff7"
		"ef6bb5261f791aa7b45e3e653d151b95dad3d93ca0290ef2";

	char *str_yr =
		"ac7dee41d8c5f4a7d5836960a773cfc1376289d3373f8cf7"
		"417b0c6207ac32e913856612fc9ff2e357eb2ee05cf9667f";

	p = ec_hexstr_to_bignum(str_p384);

	S.x = ec_hexstr_to_bignum(str_xs);
	S.y = ec_hexstr_to_bignum(str_ys);
	S.z = ec_int(1);

	d = ec_hexstr_to_bignum(str_d);

	R.x = NULL;
	R.y = NULL;
	R.z = NULL;

	ec_mult(&R, d, &S, p);

	ec_affinify(&R, p);

	x = ec_hexstr_to_bignum(str_xr);
	y = ec_hexstr_to_bignum(str_yr);

	if (ec_cmp(R.x, x) == 0 && ec_cmp(R.y, y) == 0)
		printf("ok: ec_mult\n");
	else
		printf("err: ec_mult\n");

	ec_free(p);
	ec_free(d);
	ec_free(x);
	ec_free(y);

	ec_free_xyz(&R);
	ec_free_xyz(&S);
}

void
ec_test_twin_mult()
{
	uint32_t *d, *e, *p, *x, *y;
	struct point R, S, T;

	char *str_p384 =
		"ffffffffffffffffffffffffffffffffffffffffffffffff"
		"fffffffffffffffeffffffff0000000000000000ffffffff";

	char *str_xs =
		"fba203b81bbd23f2b3be971cc23997e1ae4d89e69cb6f923"
		"85dda82768ada415ebab4167459da98e62b1332d1e73cb0e";

	char *str_ys =
		"5ffedbaefdeba603e7923e06cdb5d0c65b22301429293376"
		"d5c6944e3fa6259f162b4788de6987fd59aed5e4b5285e45";

	char *str_xt =
		"aacc05202e7fda6fc73d82f0a66220527da8117ee8f8330e"
		"ad7d20ee6f255f582d8bd38c5a7f2b40bcdb68ba13d81051";

	char *str_yt =
		"84009a263fefba7c2c57cffa5db3634d286131afc0fca8d2"
		"5afa22a7b5dce0d9470da89233cee178592f49b6fecb5092";

	char *str_d =
		"a4ebcae5a665983493ab3e626085a24c104311a761b5a8fd"
		"ac052ed1f111a5c44f76f45659d2d111a61b5fdd97583480";

	char *str_e =
		"afcf88119a3a76c87acbd6008e1349b29f4ba9aa0e12ce89"
		"bcfcae2180b38d81ab8cf15095301a182afbc6893e75385d";

	char *str_xr =
		"917ea28bcd641741ae5d18c2f1bd917ba68d34f0f0577387"
		"dc81260462aea60e2417b8bdc5d954fc729d211db23a02dc";

	char *str_yr =
		"1a29f7ce6d074654d77b40888c73e92546c8f16a5ff6bcbd"
		"307f758d4aee684beff26f6742f597e2585c86da908f7186";

	p = ec_hexstr_to_bignum(str_p384);

	S.x = ec_hexstr_to_bignum(str_xs);
	S.y = ec_hexstr_to_bignum(str_ys);
	S.z = ec_int(1);

	T.x = ec_hexstr_to_bignum(str_xt);
	T.y = ec_hexstr_to_bignum(str_yt);
	T.z = ec_int(1);

	d = ec_hexstr_to_bignum(str_d);
	e = ec_hexstr_to_bignum(str_e);

	ec_twin_mult(&R, d, &S, e, &T, p);

	ec_affinify(&R, p);

	x = ec_hexstr_to_bignum(str_xr);
	y = ec_hexstr_to_bignum(str_yr);

	if (ec_cmp(R.x, x) == 0 && ec_cmp(R.y, y) == 0)
		printf("ok: ec_twin_mult\n");
	else
		printf("err: ec_twin_mult\n");

	ec_free(p);
	ec_free(d);
	ec_free(e);
	ec_free(x);
	ec_free(y);

	ec_free_xyz(&R);
	ec_free_xyz(&S);
	ec_free_xyz(&T);
}
