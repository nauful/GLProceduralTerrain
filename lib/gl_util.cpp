#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

#include "gl_util.h"
#include "math_util.h"
#include "lexer.h"

GLuint gl_shader_create(const char* file_path, GLenum type, const char* defines)
{
	char* src = file_read(file_path);
	if (!src)
	{
		printf("gl_shader_create: Failed to read %s\n", file_path);
		exit(EXIT_FAILURE);
	}

	if (defines)
	{
		int new_len = strlen(src) + strlen(defines);
		char* new_src = new char[new_len + 1];
		memcpy(new_src, defines, strlen(defines));
		memcpy(new_src + strlen(defines), src, strlen(src));
		new_src[new_len] = 0;
		delete[] src;
		src = new_src;
	}

	GLint src_len = strlen(src);
	GLuint shader = glCreateShader(type);
	glShaderSource(shader, 1, (const GLchar**)&src, &src_len);
	glCompileShader(shader);

	delete[] src;

	GLint shader_good = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_good);

	if (!shader_good)
	{
		GLint log_length = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &log_length);
		if (log_length > 0)
		{
			char* log = new char[log_length + 1];
			glGetShaderInfoLog(shader, log_length, NULL, log);
			printf("gl_shader_create: Failed to compile %s\n%s\n", file_path, log);
			delete[] log;
		}
		else
		{
			printf("gl_shader_create: Failed to compile %s, reason unknown\n", file_path);
		}
		exit(EXIT_FAILURE);
	}

	return shader;
}

GLuint gl_program_create(GLuint vs, GLuint fs)
{
	GLuint program = glCreateProgram();
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);

	GLint program_good = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &program_good);

	if (!program_good)
	{
		GLint log_length = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &log_length);
		if (log_length > 0)
		{
			char* log = new char[log_length + 1];
			glGetProgramInfoLog(program, log_length, NULL, log);
			printf("gl_program_create: Failed to link\n%s\n", log);
			delete[] log;
		}
		else
		{
			printf("gl_program_create: Failed to link, reason unknown\n");
		}
		exit(EXIT_FAILURE);
	}

	return program;
}

GLuint gl_buffer_create(GLenum target, const void* data, GLsizei buffer_size, GLenum usage)
{
    GLuint buffer = 0;
    glGenBuffers(1, &buffer);
    if (data)
	{
		glBindBuffer(target, buffer);
		glBufferData(target, buffer_size, data, usage);
	}
    return buffer;
}

GLuint gl_image_load(const char* file_path)
{
	assert(file_path);
	FILE* f = fopen(file_path, "rb");
	if (!f)
	{
		printf("gl_image_load: Unable to open %s\n", file_path);
		assert(f);
		exit(EXIT_FAILURE);
	}

	int tex_w, tex_h;
	unsigned char* image_data = image_load(f, &tex_w, &tex_h);
	fclose(f);
	if (!image_data) { return 0; }

	GLuint img_tex;
	glGenTextures(1, &img_tex);
	glBindTexture(GL_TEXTURE_2D, img_tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_w, tex_h, 0, GL_RGB, GL_UNSIGNED_BYTE, image_data);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	//delete yuv_conv;
	glGenerateMipmap(GL_TEXTURE_2D);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	GLfloat max_anisotropy = 0;
	glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &max_anisotropy);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, max_anisotropy);
	delete[] image_data;

	return img_tex;
}

void gl_text_draw(const char* str, int x, int y, int screen_w, int screen_h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glScalef(1.0f / screen_w, -1.0f / screen_h, 1.0f);
	glTranslatef(-screen_w, -screen_h, 0.0f);

	glRasterPos2i(x, y);	
	while (*str)
	{
		if (*str == '\n')
		{
			y += 30;
			glRasterPos2i(x, y);
		}
		else
		{
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *str);
		}
		str++;
	}
}

GLFont gl_font_load(const char* file_path, const char* shader_path)
{
	GLFont font;
	FILE* file = fopen(file_path, "rb");
	if (!file)
	{
		printf("gl_font_load: Unable to open %s\n", file_path);
		exit(EXIT_FAILURE);
	}

	int version = 0;
	fread(&version, sizeof(int), 1, file);
	fread(&font.tex_w, sizeof(int), 1, file);
	fread(&font.tex_h, sizeof(int), 1, file);
	fread(&font.chars, sizeof(GLCharacter), 128, file);
	unsigned char* uncompressed_distance_field = (unsigned char*)malloc(font.tex_w * font.tex_h);

	const int word_length = 8;
	const int rle_code = (1 << word_length) - 1;

	int write_pos = 0;
	while (!feof(file))
	{
		int symbol = getc(file);
		if (symbol == rle_code)
		{
			symbol = getc(file);
			int count = getc(file) + 1;
			if (symbol != rle_code) { count += 3; }
			while (count-- > 0) { uncompressed_distance_field[write_pos++] = symbol; }
		}
		else
		{
			uncompressed_distance_field[write_pos++] = symbol;
		}
	}

	if (write_pos < font.tex_h * font.tex_w)
	{
		printf("gl_font_load: Couldn't read entire distance field texture from %s\n", file_path);
		exit(EXIT_FAILURE);
	}

	glGenTextures(1, &font.tex);
	glBindTexture(GL_TEXTURE_2D, font.tex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, font.tex_w, font.tex_h, 0, GL_RED, GL_UNSIGNED_BYTE, uncompressed_distance_field);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	font.shader = gl_program_create(
		gl_shader_create(shader_path, GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n"),
		gl_shader_create(shader_path, GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n"));

	free(uncompressed_distance_field);
	fclose(file);
	return font;
}

void gl_font_draw(GLFont* font, const char* str, int x, int y, float size, int screen_w, int screen_h)
{
	const static int default_line_spacing = 64;
	const static int default_line_space_spacing = 16;

	if (!str || !*str) { return; }

	struct FontVertex
	{
		float x, y;
		float u, v;
		float r, g, b;
	};

	Matrix screen_transform, screen_scale, screen_vp;
	screen_transform.Translate(Vector3(-1.0f, 1.0f, 0.0f));
	screen_scale.Scale(2.0f / screen_w, -2.0f / screen_h, 0);
	screen_transform = screen_scale * screen_transform;

	size_t str_len = strlen(str);
	FontVertex* vertices = (FontVertex*)malloc(4 * str_len * sizeof(FontVertex));
	GLushort* indices = (GLushort*)malloc(6 * str_len * sizeof(GLushort));
	int cur_vert = 0;
	int cur_idx = 0;

	float cur_x = x;
	float cur_y = y;

	float sdf_min_t = size > 1.0f ? 0.5f : 0.5f * powf(size, 1.0f / 12.0f);
	float antialias_t = 0.03f;
	if (size < 1.0f) { antialias_t /= max(0.3f, 3.0f * size); }
	if (size > 1.0f) { antialias_t /= size; }
	float sdf_boost_t = 1.0f;

	glUseProgram(font->shader);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, font->shader);
	GLuint uniform_wvp_loc = glGetUniformLocation(font->shader, "wvp");
	GLuint uniform_basetex_loc = glGetUniformLocation(font->shader, "basetex");
	GLuint uniform_sdf_min_t_loc = glGetUniformLocation(font->shader, "sdf_min_t");
	GLuint uniform_antialias_t_loc = glGetUniformLocation(font->shader, "antialias_t");
	GLuint uniform_sdf_boost_t_loc = glGetUniformLocation(font->shader, "sdf_boost_t");

	GLuint attrib_pos_loc = glGetAttribLocation(font->shader, "position");
	GLuint attrib_texcoord_loc = glGetAttribLocation(font->shader, "texcoord");
	GLuint attrib_color_loc = glGetAttribLocation(font->shader, "vertex_color");

	float cur_color[3] = { 1, 1, 1 };
	const char* color_table = "0123456789ABCDEF";

	for (const char* c = str; *c; c++)
	{
		if (c[0] == '$' && c[1] == 'c')
		{
			const char* red_value = lex_strchr(color_table, c[2]);
			const char* green_value = lex_strchr(color_table, c[3]);
			const char* blue_value = lex_strchr(color_table, c[4]);
			if (red_value && green_value && blue_value)
			{
				cur_color[0] = (red_value - color_table) / 15.0f;
				cur_color[1] = (green_value - color_table) / 15.0f;
				cur_color[2] = (blue_value - color_table) / 15.0f;
				c += 4;
				continue;
			}
		}

		if (*c == ' ')
		{
			cur_x += default_line_space_spacing * size;
			continue;
		}
		if (*c == '\n')
		{
			cur_x = x;
			cur_y += default_line_spacing * size;
			continue;
		}
		if (*c < ' ') { continue; }

		const GLCharacter& cv = font->chars[*c];
		float cur_w = cv.w;
		float cur_h = cv.h;
		cur_w *= size;
		cur_h *= size;

		float min_u = (float)cv.x / (float)font->tex_w;
		float max_u = (float)(cv.x + cv.w) / (float)font->tex_w;
		float min_v = (float)cv.y / (float)font->tex_h;
		float max_v = (float)(cv.y + cv.h) / (float)font->tex_h;

		float min_x = cur_x;
		float max_x = cur_x + cur_w;
		float min_y = cur_y;
		float max_y = cur_y + cur_h;

		vertices[cur_vert].x = min_x;
		vertices[cur_vert].y = max_y;
		vertices[cur_vert].u = min_u;
		vertices[cur_vert].v = max_v;
		vertices[cur_vert].r = cur_color[0];
		vertices[cur_vert].g = cur_color[1];
		vertices[cur_vert].b = cur_color[2];

		vertices[cur_vert + 1].x = max_x;
		vertices[cur_vert + 1].y = max_y;
		vertices[cur_vert + 1].u = max_u;
		vertices[cur_vert + 1].v = max_v;
		vertices[cur_vert + 1].r = cur_color[0];
		vertices[cur_vert + 1].g = cur_color[1];
		vertices[cur_vert + 1].b = cur_color[2];

		vertices[cur_vert + 2].x = max_x;
		vertices[cur_vert + 2].y = min_y;
		vertices[cur_vert + 2].u = max_u;
		vertices[cur_vert + 2].v = min_v;
		vertices[cur_vert + 2].r = cur_color[0];
		vertices[cur_vert + 2].g = cur_color[1];
		vertices[cur_vert + 2].b = cur_color[2];

		vertices[cur_vert + 3].x = min_x;
		vertices[cur_vert + 3].y = min_y;
		vertices[cur_vert + 3].u = min_u;
		vertices[cur_vert + 3].v = min_v;
		vertices[cur_vert + 3].r = cur_color[0];
		vertices[cur_vert + 3].g = cur_color[1];
		vertices[cur_vert + 3].b = cur_color[2];

		cur_x += cur_w + 2;

		indices[cur_idx++] = cur_vert;
		indices[cur_idx++] = cur_vert + 1;
		indices[cur_idx++] = cur_vert + 2;
		indices[cur_idx++] = cur_vert;
		indices[cur_idx++] = cur_vert + 2;
		indices[cur_idx++] = cur_vert + 3;
		cur_vert += 4;
	}

	if (cur_idx > 0)
	{
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, font->tex);

		if (uniform_basetex_loc != -1) { glUniform1i(uniform_basetex_loc, 0); }
		if (uniform_sdf_min_t_loc != -1) { glUniform1f(uniform_sdf_min_t_loc, sdf_min_t); }
		if (uniform_antialias_t_loc != -1) { glUniform1f(uniform_antialias_t_loc, antialias_t); }
		if (uniform_sdf_boost_t_loc != -1) { glUniform1f(uniform_sdf_boost_t_loc, sdf_boost_t); }

		if (uniform_wvp_loc != -1) { glUniformMatrix4fv(uniform_wvp_loc, 1, true, (const GLfloat*)screen_transform.cell); }

		GLuint vb = gl_buffer_create(GL_ARRAY_BUFFER, vertices, sizeof(FontVertex) * cur_vert, GL_STATIC_DRAW);
		GLuint ib = gl_buffer_create(GL_ELEMENT_ARRAY_BUFFER, indices, sizeof(GLushort) * cur_idx, GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, vb);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib);

		glEnableVertexAttribArray(attrib_pos_loc);
		glEnableVertexAttribArray(attrib_texcoord_loc);
		glEnableVertexAttribArray(attrib_color_loc);

		glVertexAttribPointer(attrib_pos_loc, 2, GL_FLOAT, GL_FALSE, sizeof(FontVertex), (void*)0);
		glVertexAttribPointer(attrib_texcoord_loc, 2, GL_FLOAT, GL_FALSE, sizeof(FontVertex), (void*)(sizeof(float) * 2));
		glVertexAttribPointer(attrib_color_loc, 3, GL_FLOAT, GL_FALSE, sizeof(FontVertex), (void*)(sizeof(float) * 4));

		glDrawElements(GL_TRIANGLES, cur_idx, GL_UNSIGNED_SHORT, 0);

		glDisableVertexAttribArray(attrib_pos_loc);
		glDisableVertexAttribArray(attrib_texcoord_loc);
		glDisableVertexAttribArray(attrib_color_loc);
	}

	free(vertices);
	free(indices);
}