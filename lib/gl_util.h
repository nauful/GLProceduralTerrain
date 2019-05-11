#ifndef GL_UTIL_H
#define GL_UTIL_H

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>

#include "image.h"

struct GLCharacter
{
	int value;
	int x, y, w, h;
	int ext_x, ext_y, ext_w, ext_h;
	int advance_x;

	int horz_kerning[256];
};

struct GLFont
{
	GLuint tex, shader;

	int tex_w, tex_h;
	GLCharacter chars[128];
};

GLFont gl_font_load(const char* file_path, const char* shader_path);
void gl_font_draw(GLFont* font, const char* str, int x, int y, float size, int screen_w, int screen_h);

GLuint gl_shader_create(const char* file_path, GLenum type, const char* defines);
GLuint gl_program_create(GLuint vs, GLuint fs);
GLuint gl_buffer_create(GLenum target, const void* data, GLsizei buffer_size, GLenum usage);
GLuint gl_image_load(const char* file_path);

void gl_text_draw(const char* str, int x, int y, int screen_w, int screen_h);

#endif