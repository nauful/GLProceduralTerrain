#include "lexer.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

const char* lex_source_buffer = 0;
int lex_source_buffer_len = 0;

char* file_read(const char* path)
{
	assert(path);

	FILE* f = fopen(path, "rb");
	assert(f);

	fseek(f, 0, SEEK_END);
	size_t size = ftell(f);
	fseek(f, 0, SEEK_SET);

	char* buf = new char[size + 1];
	assert(buf);
	fread(buf, 1, size, f);
	buf[size] = '\0';

	fclose(f);
	return buf;
}

bool is_whitespace(int c)
{
	return strchr(" \t\r\n", c) != 0;
}

void lex_set_source_buffer(const char* src)
{
	lex_source_buffer = src;

	if (src)
	{
		lex_source_buffer_len = strlen(src);
	}
}

void lex_err_trace(const char* src)
{
	const static int err_source_len = 50;
	char buf[err_source_len + 1];
	int srclen;
	int i;
	int pos_line, pos_column;

	if (!src)
	{
		printf("src is null\n");
		return;
	}

	if (lex_source_buffer && src >= lex_source_buffer && src < lex_source_buffer + lex_source_buffer_len)
	{
		pos_column = 1;
		pos_line = 1;
		for (i = 0; i < (int)(src - lex_source_buffer); i++)
		{
			if (lex_source_buffer[i] == '\n')
			{
				pos_column = 1;
				pos_line++;
			}
			else
			{
				pos_column++;
			}
		}
		printf("lexer position: line %d, column %d\n", pos_line, pos_column);
	}
	
	srclen = 0;
	src += skip_whitespace(src);
	for (i = 0; src[i] && !is_whitespace(src[i]) && srclen <= err_source_len; i++)
	{
		srclen++;
	}

	if (srclen >= err_source_len)
	{
		memcpy(buf, src, err_source_len - 4);
		buf[err_source_len - 4] = ' ';
		buf[err_source_len - 3] = '.';
		buf[err_source_len - 2] = '.';
		buf[err_source_len - 1] = '.';
		buf[err_source_len - 0] = '\0';
	}
	else
	{
		memcpy(buf, src, srclen);
		buf[srclen] = '\0';
	}

	printf("lexer error line: %s\n", buf);
}

int skip_whitespace(const char* str)
{
	const char* start = str;
	while (*str)
	{
		if (is_whitespace(*str)) { str++; }
		// Single-line comment
		else if (str[0] == '/' && str[1] == '/')
		{
			while (*str && *str != '\n') { str++; }
			if (*str == '\n') { str++; }
		}
		// Multi-line comment
		else if (str[0] == '/' && str[1] == '*')
		{
			while (*str && !(str[0] == '*' && str[1] == '/')) { str++; }
			if (str[0] == '*' && str[1] == '/') { str += 2; }
		}
		else { break; }
	}
	return (int)(str - start);
}

bool strmatch(const char* s1, const char* s2, int n)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if (*s1 != *s2) { return false; }
		if (!*s1 || !*s2) { return false; }
		s1++;
		s2++;
	}
	return true;
}

int lex_token(const char* src, char* dest, int dest_len)
{
	// dest must have dest_len + 1 bytes
	const char* cur;
	// Initial and terminal flag
	char flag = 0;

	if (!dest) { dest_len = 0; }
	cur = src;
	cur += skip_whitespace(cur);

	if (*cur == '\"')
	{
		flag = *cur;
		cur++;
		cur += skip_whitespace(cur);
	}

	while (*cur && (!is_whitespace(*cur) || flag) && (!flag || *cur != flag))
	{
		if (dest)
		{
			if (dest_len < 1)
			{
				printf("lex_token: Buffer overflow\n");
				lex_err_trace(src);
				exit(EXIT_FAILURE);
			}
			dest_len--;
			*dest++ = *cur++;
		}
		else
		{
			cur++;
		}
	}

	if (flag && *cur == flag) { cur++; }
	if (dest) { *dest = '\0'; }

	return (int)(cur - src);
}

bool lex_is_end(const char* src)
{
	const char* cur;

	cur = src;
	cur += skip_whitespace(cur);
	return *cur == '\0';
}

int lex_int(const char* src, int* res)
{
	const char* cur;
	char buf[128];
	char* dest;
	int dest_len;
	bool read_digits;

	dest = buf;
	dest_len = sizeof(buf);
	*dest = '\0';
	read_digits = false;

	cur = src;
	cur += skip_whitespace(cur);

	if (*cur && strchr("+-", *cur))
	{
		--dest_len;
		*dest++ = *cur++;
	}
	while (*cur && strchr("0123456789", *cur))
	{
		if (dest_len < 1)
		{
			printf("lex_int: Buffer overflow\n");
			lex_err_trace(src);
			exit(EXIT_FAILURE);
		}
		--dest_len;
		*dest++ = *cur++;
		read_digits = true;
	}
	*dest = '\0';

	if (!read_digits)
	{
		printf("lex_int: Expected int\n");
		lex_err_trace(src);
		exit(EXIT_FAILURE);
	}

	if (res) { *res = *buf ? atoi(buf) : 0; }
	return (int)(cur - src);
}

int lex_float(const char* src, float* res, bool must_read)
{
	const char* cur;
	char buf[128];
	char* dest;
	int dest_len;
	bool read_digits;
	int period_count;
	const char* exponent_pos;

	dest = buf;
	dest_len = sizeof(buf);
	*dest = '\0';
	read_digits = false;
	period_count = 0;
	exponent_pos = 0;

	cur = src;
	cur += skip_whitespace(cur);

	if (*cur && strchr("+-", *cur))
	{
		--dest_len;
		*dest++ = *cur++;
	}
	while (*cur)
	{
		if (period_count == 0 && *cur == '.') { ++period_count; }
		else if (strchr("0123456789", *cur)) { }
		else if (!exponent_pos && strchr("eE", *cur)) { exponent_pos = cur; }
		else if (exponent_pos == cur - 1 && strchr("+-", *cur)) { }
		else if (strchr("fF", *cur))
		{
			cur++;
			break;
		}
		else { break; }

		if (dest_len < 1)
		{
			printf("lex_float: Buffer overflow\n");
			lex_err_trace(src);
			exit(EXIT_FAILURE);
		}
		--dest_len;
		*dest++ = *cur++;
		read_digits = true;
	}
	*dest = '\0';

	if (!read_digits)
	{
		if (must_read)
		{
			printf("lex_float: Expected float\n");
			lex_err_trace(src);
			exit(EXIT_FAILURE);
		}
		else { return 0; }
	}

	if (res) { *res = *buf ? atof(buf) : 0; }
	return (int)(cur - src);
}

lexer_punc_t* lex_get_punc(const char* src, lexer_punc_t* table)
{
	const char* cur;
	int best_match_len;
	lexer_punc_t* best_match;

	src += skip_whitespace(src);
	best_match = 0;
	while (table && table->str)
	{
		cur = src;
		if (*cur && strmatch(cur, table->str, table->str_len))
		{
			if (!best_match || table->str_len > best_match_len)
			{
				best_match_len = table->str_len;
				best_match = table;
			}
		}

		table++;
	}

	return best_match;
}

int lex_consume_punc(const char* src, lexer_punc_t* punc)
{
	const char* cur;

	cur = src;
	cur += skip_whitespace(cur);
	if (*cur && strmatch(cur, punc->str, punc->str_len))
	{
		cur += punc->str_len;
		return (int)(cur - src);
	}
	else
	{
		printf("lex_consume_punc: Unexpected punctuation, expected %s\n", punc->str);
		lex_err_trace(src);
		exit(EXIT_FAILURE);
	}

	return 0;
}

int lex_end_line(const char* src)
{
	const char* cur;

	cur = src;
	while (*cur && *cur != '\n') { cur++; }
	if (*cur == '\n') { cur++; }

	return (int)(cur - src);
}

char lex_peek_symbol(const char* src)
{
	const char* cur;

	cur = src;
	cur += skip_whitespace(cur);
	return *cur;
}

int lex_consume_symbol(const char* src, char sym)
{
	const char* cur;

	cur = src;
	cur += skip_whitespace(cur);
	if (*cur && *cur == sym)
	{
		cur++;
		return (int)(cur - src);
	}
	else
	{
		printf("lex_consume_symbol: Unexpected symbol, expected %c\n", sym);
		lex_err_trace(src);
		exit(EXIT_FAILURE);
	}

	return 0;
}

const char* lex_strchr(const char* src, char c)
{
	const char* cur;

	cur = src;
	while (*cur)
	{
		if (*cur == c) { return cur; }
		cur++;
	}

	return 0;
}