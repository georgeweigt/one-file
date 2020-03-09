/* Print the contents of an RSA key file.

To build and run:

	gcc print-rsa-keys.c
	./a.out

To generate the key file:

	openssl genrsa -out key.pem 2048
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define INTEGER 2
#define BIT_STRING 3
#define OCTET_STRING 4
#define OID 6
#define UTF8STRING 12
#define PRINTABLE_STRING 19
#define IA5STRING 22
#define UTCTIME 0x17
#define GENERALIZEDTIME 0x18
#define SEQUENCE 0x30
#define SET 0x31

struct keyinfo {

	int modulus_offset;
	int modulus_length;

	int public_exponent_offset;
	int public_exponent_length;

	int private_exponent_offset;
	int private_exponent_length;

	int prime1_offset;
	int prime1_length;

	int prime2_offset;
	int prime2_length;

	int exponent1_offset;
	int exponent1_length;

	int exponent2_offset;
	int exponent2_length;

	int coefficient_offset;
	int coefficient_length;

	int key_data_length;
	unsigned char key_data[10000];
};

int read_key_file(char *filename, char *passphrase, struct keyinfo *ki);
int decode(int c);
int get_rsa_keys_from_keyinfo(struct keyinfo *ki);
int get_type_and_length(uint8_t *buf, int end, int *offset, int *type, int *length);

char strbuf[100];
struct keyinfo provisioned_ki;

int
read_key_file(char *filename, char *passphrase, struct keyinfo *ki)
{
	int a, b, c, d, k, n, w;
	char *s;
	FILE *f;
	f = fopen(filename, "r");
	if (f == NULL)
		return -1;
	for (;;) {
		if (fgets(strbuf, sizeof strbuf, f) == NULL)
			goto err;
		if (strcmp(strbuf, "-----BEGIN RSA PRIVATE KEY-----\n") == 0)
			break;
	}
#if 0
	int i;
	uint8_t iv[16], key[16], expanded_key[272];
	if (passphrase) { // encrypted key file
		if (fgets(strbuf, sizeof strbuf, f) == NULL)
			goto err;
		if (strcmp(strbuf, "Proc-Type: 4,ENCRYPTED\n") != 0)
			goto err;
		if (fgets(strbuf, sizeof strbuf, f) == NULL)
			goto err;
		if (strlen(strbuf) != 55)
			goto err;
		if (strncmp(strbuf, "DEK-Info: AES-128-CBC,", 22) != 0)
			goto err;
		for (i = 0; i < 16; i++) { // parse iv
			a = strbuf[22 + 2 * i];
			b = strbuf[23 + 2 * i];
			if ('0' <= a && a <= '9')
				a = a - '0';
			else if ('A' <= a && a <= 'F')
				a = a - 'A' + 10;
			else if ('a' <= a && a <= 'f')
				a = a - 'a' + 10;
			else
				goto err;
			if ('0' <= b && b <= '9')
				b = b - '0';
			else if ('A' <= b && b <= 'F')
				b = b - 'A' + 10;
			else if ('a' <= b && b <= 'f')
				b = b - 'a' + 10;
			else
				goto err;
			iv[i] = a << 4 | b;
		}
		n = strlen(passphrase);
		strncpy(strbuf, passphrase, n);
		memcpy(strbuf + n, iv, 8);
		md5(strbuf, n + 8, key);
		if (fgets(strbuf, sizeof strbuf, f) == NULL) // skip blank line
			goto err;
	}
#endif
	k = 0;
	for (;;) {
		if (fgets(strbuf, sizeof strbuf, f) == NULL)
			goto err;
		if (strcmp(strbuf, "-----END RSA PRIVATE KEY-----\n") == 0)
			goto done;
		n = (int) strlen(strbuf) - 1;
		if (n % 4)
			goto err;
		if (k + 3 * n / 4 > sizeof ki->key_data)
			goto err;
		s = strbuf;
		while (s < strbuf + n) {
			a = decode(s[0]);
			b = decode(s[1]);
			c = decode(s[2]);
			d = decode(s[3]);
			// check for padding
			if (s[2] == '=' || s[3] == '=')
				goto pad;
			if (a == -1 || b == -1 || c == -1 || d == -1)
				goto err;
			w = a << 18 | b << 12 | c << 6 | d;
			ki->key_data[k++] = w >> 16;
			ki->key_data[k++] = w >> 8;
			ki->key_data[k++] = w;
			s += 4;
		}
	}
pad:	if (s[2] == '=') {
		if (a == -1 || b == -1 || s[3] != '=')
			goto err;
		w = a << 18 | b << 12;
		ki->key_data[k++] = w >> 16;
	} else {
		if (a == -1 || b == -1 || c == -1)
			goto err;
		w = a << 18 | b << 12 | c << 6;
		ki->key_data[k++] = w >> 16;
		ki->key_data[k++] = w >> 8;
	}
	// must be end of string and file
	if (s + 4 != strbuf + n)
		goto err;
	if (fgets(strbuf, sizeof strbuf, f) == NULL)
		goto err;
	if (strcmp(strbuf, "-----END RSA PRIVATE KEY-----\n") != 0)
		goto err;
done:	ki->key_data_length = k;
	fclose(f);
#if 0
	if (passphrase) {
		if (ki->key_data_length & 0xf)
			return -1;
		aes128_decrypt_direct(key, iv, ki->key_data, ki->key_data_length / 16);
		n = ki->key_data[ki->key_data_length - 1]; // pad format is different from rsa cipher suites
		ki->key_data_length -= n; // undo pad
		if (ki->key_data_length < 0)
			return -1;
	}
#endif
	return 0;
err:	fclose(f);
	return -1;
}

// convert base-64 digit to binary

int
decode(int c)
{
	if ('A' <= c && c <= 'Z')
		return c - 'A';
	if ('a' <= c && c <= 'z')
		return c - 'a' + 26;
	if ('0' <= c && c <= '9')
		return c - '0' + 52;
	if (c == '+')
		return 62;
	if (c == '/')
		return 63;
	return -1;
}

int
get_rsa_keys_from_keyinfo(struct keyinfo *ki)
{
	int end, err, length, offset, type;

	end = ki->key_data_length;
	offset = 0;

	// SEQUENCE

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != SEQUENCE)
		return -1;
	end = offset + length; // cannot go past end of this sequence

	// version

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	offset += length;

	// modulus

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->modulus_offset = offset;
	ki->modulus_length = length;
	offset += length;

	// public exponent

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->public_exponent_offset = offset;
	ki->public_exponent_length = length;
	offset += length;

	// private exponent

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->private_exponent_offset = offset;
	ki->private_exponent_length = length;
	offset += length;

	// prime1

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->prime1_offset = offset;
	ki->prime1_length = length;
	offset += length;

	// prime2

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->prime2_offset = offset;
	ki->prime2_length = length;
	offset += length;

	// exponent1

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->exponent1_offset = offset;
	ki->exponent1_length = length;
	offset += length;

	// exponent2

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->exponent2_offset = offset;
	ki->exponent2_length = length;
	offset += length;

	// coefficient

	err = get_type_and_length(ki->key_data, end, &offset, &type, &length);
	if (err || type != INTEGER)
		return -1;
	ki->coefficient_offset = offset;
	ki->coefficient_length = length;
	offset += length;

	return 0;
}

// this function advances offset to beginning of V in TLV

int
get_type_and_length(uint8_t *buf, int end, int *offset, int *type, int *length)
{
	int i, k, l, n, t;
	k = *offset;
	if (k + 2 > end)
		return -1;
	t = buf[k++];
	l = buf[k++];
	if (l < 128) {
		if (k + l > end)
			return -1;
		*offset = k;
		*type = t;
		*length = l;
		return 0;
	}
	n = l & 0x7f;
	if (k + n > end)
		return -1;
	l = 0;
	for (i = 0; i < n; i++) {
		if (l & 0xff000000)
			return -1;
		l = l << 8 | buf[k++];
	}
	if (k + l > end)
		return -1;
	*offset = k;
	*type = t;
	*length = l;
	return 0;
}

void print_field(char *s, uint8_t *p, int length);

void
print_key_info(struct keyinfo *p)
{
	print_field("modulus", p->key_data + p->modulus_offset, p->modulus_length);

	print_field("public exponent", p->key_data + p->public_exponent_offset, p->public_exponent_length);
	print_field("private exponent", p->key_data + p->private_exponent_offset, p->private_exponent_length);

	print_field("prime1", p->key_data + p->prime1_offset, p->prime1_length);
	print_field("prime2", p->key_data + p->prime2_offset, p->prime2_length);

	print_field("exponent1", p->key_data + p->exponent1_offset, p->exponent1_length);
	print_field("exponent2", p->key_data + p->exponent2_offset, p->exponent2_length);

	print_field("coefficient", p->key_data + p->coefficient_offset, p->coefficient_length);
}

void
print_field(char *s, uint8_t *p, int length)
{
	int i;

	printf("%s (%d bytes)\n", s, length);

	for (i = 0; i < length; i++) {
		printf("%02x", p[i]);
		if (i % 16 == 15)
			printf("\n");
	}

	if (length % 16)
		printf("\n");

	printf("\n");
}

int
main(int argc, char *argv[])
{
	int err;

	err = read_key_file("key.pem", NULL, &provisioned_ki);

	if (err) {
		printf("key.pem error\n");
		exit(1);
	}

	err = get_rsa_keys_from_keyinfo(&provisioned_ki);

	if (err) {
		printf("key.pem error\n");
		exit(1);
	}

	print_key_info(&provisioned_ki);
}
