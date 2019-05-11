#ifndef LEXER_H
#define LEXER_H

char* file_read(const char* path);

typedef struct
{
	const char* str;
	int str_len;
	int value;
} lexer_punc_t;

bool is_whitespace(int c);
int skip_whitespace(const char* str);
bool strmatch(const char* s1, const char* s2, int n);

bool lex_is_end(const char* src);
int lex_token(const char* src, char* dest, int dest_len);
int lex_int(const char* src, int* res);
int lex_float(const char* src, float* res, bool must_read = true);
int lex_end_line(const char* src);

lexer_punc_t* lex_get_punc(const char* src, lexer_punc_t* table);
int lex_consume_punc(const char* src, lexer_punc_t* punc);
int lex_consume_symbol(const char* src, char sym);
char lex_peek_symbol(const char* src);

const char* lex_strchr(const char* src, char c);

void lex_err_trace(const char* src);
void lex_set_source_buffer(const char* source_buffer);

#endif