/* This is a realization of the cipher DES.

The main function runs a self test.

To build and run:

	gcc des.c
	./a.out
*/

#include <stdio.h>
#include <string.h>
#include <stdint.h>

void des_compute_key_schedule(uint8_t *K, uint8_t *KEY);
void des_permuted_choice_1(uint8_t *C, uint8_t *D, uint8_t *KEY);
void des_permuted_choice_2(uint8_t *K, uint8_t *C, uint8_t *D);
void des_left_shift(uint8_t *a);
void des_encrypt(uint8_t *a, uint8_t *K);
void des_decrypt(uint8_t *a, uint8_t *K);
void des_initial_permutation(uint8_t *L, uint8_t *R, uint8_t *a);
void des_inverse_initial_permutation(uint8_t *a, uint8_t *R, uint8_t *L);
void des_func(uint8_t *R, uint8_t *L, uint8_t *K);

static uint8_t bit[8] = {0x80, 0x40, 0x20, 0x10, 8, 4, 2, 1};

// K is 96 bytes, KEY is 8 bytes

void
des_compute_key_schedule(uint8_t *K, uint8_t *KEY)
{
	int i;
	uint8_t C[4], D[4];

	static int s[16] = {1,1,2,2,2,2,2,2,1,2,2,2,2,2,2,1};

	des_permuted_choice_1(C, D, KEY);

	for (i = 0; i < 16; i++) {
		des_left_shift(C);
		des_left_shift(D);
		if (s[i] == 2) {
			des_left_shift(C);
			des_left_shift(D);
		}
		des_permuted_choice_2(K + 6 * i, C, D);
	}
}

void
des_permuted_choice_1(uint8_t *C, uint8_t *D, uint8_t *KEY)
{
	int i, j;

	static int x[28] = {
	57, 49, 41, 33, 25, 17,  9,
	 1, 58, 50, 42, 34, 26, 18,
	10,  2, 59, 51, 43, 35, 27,
	19, 11,  3, 60, 52, 44, 36};

	static int y[28] = {
	63, 55, 47, 39, 31, 23, 15,
	 7, 62, 54, 46, 38, 30, 22,
	14,  6, 61, 53, 45, 37, 29,
	21, 13,  5, 28, 20, 12,  4};

	bzero(C, 4);
	bzero(D, 4);

	for (i = 0; i < 28; i++) {
		j = x[i] - 1;
		if (KEY[j / 8] & bit[j % 8])
			C[i / 8] |= bit[i % 8];
		j = y[i] - 1;
		if (KEY[j / 8] & bit[j % 8])
			D[i / 8] |= bit[i % 8];
	}
}

void
des_permuted_choice_2(uint8_t *K, uint8_t *C, uint8_t *D)
{
	int i, j;

	static int x[48] = {
	14, 17, 11, 24,  1,  5,
	 3, 28, 15,  6, 21, 10,
	23, 19, 12,  4, 26,  8,
	16,  7, 27, 20, 13,  2,
	41, 52, 31, 37, 47, 55,
	30, 40, 51, 45, 33, 48,
	44, 49, 39, 56, 34, 53,
	46, 42, 50, 36, 29, 32};

	bzero(K, 6);

	for (i = 0; i < 48; i++) {
		j = x[i] - 1;
		if (j < 28) {
			if (C[j / 8] & bit[j % 8])
				K[i / 8] |= bit[i % 8];
		} else {
			j -= 28;
			if (D[j / 8] & bit[j % 8])
				K[i / 8] |= bit[i % 8];
		}
	}
}

void
des_left_shift(uint8_t *a)
{
	uint8_t t;

	t = a[0];

	a[0] <<= 1;
	if (a[1] & 0x80)
		a[0] |= 1;

	a[1] <<= 1;
	if (a[2] & 0x80)
		a[1] |= 1;

	a[2] <<= 1;
	if (a[3] & 0x80)
		a[2] |= 1;

	a[3] <<= 1;
	if (t & 0x80)
		a[3] |= 0x10;
}

void
des_encrypt(uint8_t *a, uint8_t *K)
{
	int i;
	uint8_t L[4], R[4], t[4];
	des_initial_permutation(L, R, a);
	for (i = 0; i < 16; i++) {
		memcpy(t, R, 4);
		des_func(R, L, K + 6 * i);
		memcpy(L, t, 4);
	}
	des_inverse_initial_permutation(a, R, L);
}

void
des_decrypt(uint8_t *a, uint8_t *K)
{
	int i;
	uint8_t L[4], R[4], t[4];
	des_initial_permutation(R, L, a);
	for (i = 0; i < 16; i++) {
		memcpy(t, L, 4);
		des_func(L, R, K + 6 * (15 - i));
		memcpy(R, t, 4);
	}
	des_inverse_initial_permutation(a, L, R);
}

void
des_initial_permutation(uint8_t *L, uint8_t *R, uint8_t *a)
{
	int i, j;

	static int x[64] = {
	58, 50, 42, 34, 26, 18, 10,  2,
	60, 52, 44, 36, 28, 20, 12,  4,
	62, 54, 46, 38, 30, 22, 14,  6,
	64, 56, 48, 40, 32, 24, 16,  8,
	57, 49, 41, 33, 25, 17,  9,  1,
	59, 51, 43, 35, 27, 19, 11,  3,
	61, 53, 45, 37, 29, 21, 13,  5,
	63, 55, 47, 39, 31, 23, 15,  7};

	bzero(L, 4);
	bzero(R, 4);

	for (i = 0; i < 32; i++) {
		j = x[i] - 1;
		if (a[j / 8] & bit[j % 8])
			L[i / 8] |= bit[i % 8];
		j = x[i + 32] - 1;
		if (a[j / 8] & bit[j % 8])
			R[i / 8] |= bit[i % 8];
	}
}

void
des_inverse_initial_permutation(uint8_t *a, uint8_t *R, uint8_t *L)
{
	int i, j;

	static int x[64] = {
	40,  8, 48, 16, 56, 24, 64, 32,
	39,  7, 47, 15, 55, 23, 63, 31,
	38,  6, 46, 14, 54, 22, 62, 30,
	37,  5, 45, 13, 53, 21, 61, 29,
	36,  4, 44, 12, 52, 20, 60, 28,
	35,  3, 43, 11, 51, 19, 59, 27,
	34,  2, 42, 10, 50, 18, 58, 26,
	33,  1, 41,  9, 49, 17, 57, 25};

	bzero(a, 8);

	for (i = 0; i < 64; i++) {
		j = x[i] - 1;
		if (j < 32) {
			if (R[j / 8] & bit[j % 8])
				a[i / 8] |= bit[i % 8];
		} else {
			j -= 32;
			if (L[j / 8] & bit[j % 8])
				a[i / 8] |= bit[i % 8];
		}
	}
}

static int E[48] = {
	32,  1,  2,  3,  4,  5,
	 4,  5,  6,  7,  8,  9,
	 8,  9, 10, 11, 12, 13,
	12, 13, 14, 15, 16, 17,
	16, 17, 18, 19, 20, 21,
	20, 21, 22, 23, 24, 25,
	24, 25, 26, 27, 28, 29,
	28, 29, 30, 31, 32,  1,
};

static int S1[4][16] = {
	14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7,
	0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8,
	4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0,
	15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13,
};

static int S2[4][16] = {
	15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10,
	3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5,
	0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15,
	13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9,
};

static int S3[4][16] = {
	10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8,
	13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1,
	13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7,
	1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12,
};

static int S4[4][16] = {
	7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15,
	13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9,
	10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4,
	3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14,
};

static int S5[4][16] = {
	2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9,
	14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6,
	4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14,
	11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3,
};

static int S6[4][16] = {
	12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11,
	10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8,
	9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6,
	4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13,
};

static int S7[4][16] = {
	4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1,
	13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6,
	1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2,
	6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12,
};

static int S8[4][16] = {
	13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7,
	1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2,
	7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8,
	2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11,
};

static int P[32] = {
	16,  7, 20, 21,
	29, 12, 28, 17,
	 1, 15, 23, 26,
	 5, 18, 31, 10,
	 2,  8, 24, 14,
	32, 27,  3,  9,
	19, 13, 30,  6,
	22, 11,  4, 25,
};

#define BLOCK1(p) ((p)[0] >> 2)
#define BLOCK2(p) ((((p)[0] << 4) & 0x30) | ((p)[1] >> 4))
#define BLOCK3(p) ((((p)[1] << 2) & 0x3c) | ((p)[2] >> 6))
#define BLOCK4(p) ((p)[2] & 0x3f)
#define LOOKUP(s, k) s[(((k) >> 4) & 2) | ((k) & 1)][((k) >> 1) & 0xf]

// R = L ^ f(R, K)

void
des_func(uint8_t *R, uint8_t *L, uint8_t *K)
{
	int i, j, k;
	uint8_t s[4], t[6];

	bzero(t, 6);

	for (i = 0; i < 48; i++) {
		j = E[i] - 1;
		if (R[j / 8] & bit[j % 8])
			t[i / 8] |= bit[i % 8];
	}

	t[0] ^= K[0];
	t[1] ^= K[1];
	t[2] ^= K[2];
	t[3] ^= K[3];
	t[4] ^= K[4];
	t[5] ^= K[5];

	k = BLOCK1(t);				// S1
	s[0] = LOOKUP(S1, k) << 4;

	k = BLOCK2(t);				// S2
	s[0] |= LOOKUP(S2, k);

	k = BLOCK3(t);				// S3
	s[1] = LOOKUP(S3, k) << 4;

	k = BLOCK4(t);				// S4
	s[1] |= LOOKUP(S4, k);

	k = BLOCK1(t + 3);			// S5
	s[2] = LOOKUP(S5, k) << 4;

	k = BLOCK2(t + 3);			// S6
	s[2] |= LOOKUP(S6, k);

	k = BLOCK3(t + 3);			// S7
	s[3] = LOOKUP(S7, k) << 4;

	k = BLOCK4(t + 3);			// S8
	s[3] |= LOOKUP(S8, k);

	bzero(t, 4);

	for (i = 0; i < 32; i++) {
		j = P[i] - 1;
		if (s[j / 8] & bit[j % 8])
			t[i / 8] |= bit[i % 8];
	}

	R[0] = L[0] ^ t[0];
	R[1] = L[1] ^ t[1];
	R[2] = L[2] ^ t[2];
	R[3] = L[3] ^ t[3];
}

int
main()
{
	uint8_t K[96];
	uint8_t a[8], key[8];

	bzero(a, 8);

	key[0] = 0x3b;
	key[1] = 0x38;
	key[2] = 0x98;
	key[3] = 0x37;
	key[4] = 0x15;
	key[5] = 0x20;
	key[6] = 0xf7;
	key[7] = 0x5e;

	des_compute_key_schedule(K, key);

	des_encrypt(a, K);

//	printf("%02x%02x%02x%02x%02x%02x%02x%02x\n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);

	if (a[0] == 0x83 && a[1] == 0xa1 && a[2] == 0xe8 && a[3] == 0x14 && a[4] == 0x88 && a[5] == 0x92 && a[6] == 0x53 && a[7] == 0xe0)
		printf("encryption test passed\n");
	else
		printf("encryption test failed\n");

	des_decrypt(a, K);

//	printf("%02x%02x%02x%02x%02x%02x%02x%02x\n", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]);

	if (a[0] == 0x00 && a[1] == 0x00 && a[2] == 0x00 && a[3] == 0x00 && a[4] == 0x00 && a[5] == 0x00 && a[6] == 0x00 && a[7] == 0x00)
		printf("decryption test passed\n");
	else
		printf("decryption test failed\n");
}
