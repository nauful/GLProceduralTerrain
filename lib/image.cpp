#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "image.h"
#include "math_util.h"

typedef unsigned char byte;

const static byte image_magic_header[] = { 'F', 'W', 'I', 'M', 'G', '\0' };
const static int image_magic_header_len = 6;

const static int CODING_MAX_BITS = 64;

const static int BLOCK_END_FLAG = (1 << 17) - 1;
const static int RLE_BITS = 10;
const static int RLE_MAX_DIST = (1 << RLE_BITS) - 1;

struct CodingNode
{
	int val, freq;
	int parent, left, right;

	int code[CODING_MAX_BITS];
	int code_len;
};

struct FreqNode
{
	int val, freq;
};

int npow2(int x)
{
	if (x < 0) { return x; }
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x + 1;
}

bool ispow2(int x) { return x > 0 && (x & (x - 1)) == 0; }

int partition(FreqNode* table, int p, int r)
{
	FreqNode pivot = table[p];
	int left = p;

	for (int i = p + 1; i <= r; i++)
	{
		if (table[i].freq > pivot.freq)
		{
			++left;
			FreqNode tmp = table[i];
			table[i] = table[left];
			table[left] = tmp;
		}
	}

	if (p != left)
	{
		FreqNode tmp = table[p];
		table[p] = table[left];
		table[left] = tmp;
	}
	return left;
}

void qsort(FreqNode* table, int p, int r)
{
	if (p < r)
	{
		int q = partition(table, p, r);
		qsort(table, p, q - 1);
		qsort(table, q + 1, r);
	}
}

int coding_node_lowest_unused(CodingNode* nodes, int start, int num, int exclude_index = -1)
{
	int best = -1;
	int best_freq = 0;

	for (int i = start; i < start + num; i++)
	{
		if (nodes[i].parent == -1 && i != exclude_index && (best == -1 || nodes[i].freq <= best_freq))
		{
			best = i;
			best_freq = nodes[i].freq;
		}
	}

	return best;
}

void coding_node_assign(CodingNode* nodes, int num)
{
	assert(num > 0);
	nodes[0].code_len = 0;

	for (int i = 1; i < num; i++)
	{
		CodingNode& n = nodes[i];
		const CodingNode& parent = nodes[n.parent];

		if (parent.code_len > 0)
		{
			n.code_len = parent.code_len + 1;
			for (int j = 0; j < parent.code_len; j++) { n.code[j] = parent.code[j]; }
		}
		else
		{
			n.code_len = 1;
		}

		assert(n.code_len < CODING_MAX_BITS);
		assert(i == parent.left || i == parent.right);
		n.code[n.code_len - 1] = i == parent.left ? 0 : 1;
	}
}

int bit_read(FILE* file, byte& cur_byte, int& cur_bits, int num_bits)
{
	int value = 0;
	int value_bits = 0;
	
	assert(!feof(file) || num_bits > 0);

	while (value_bits < num_bits)
	{
		if (cur_bits == 0)
		{
			if (feof(file)) { cur_byte = 0; }
			else { fread(&cur_byte, 1, 1, file); }
		}

		int get = 8 - cur_bits;
		if (get > (num_bits - value_bits)) { get = num_bits - value_bits; }
		int fraction = cur_byte;
		fraction >>= cur_bits;
		fraction &= (1 << get) - 1;
		value |= fraction << value_bits;
		value_bits += get;
		cur_bits = (cur_bits + get) & 7;
	}

	return value;
}

void decode_huffman_buffer(FILE* file, int* buf, int buf_len)
{
	int num_symbols;
	fread(&num_symbols, sizeof(int), 1, file);
	int* buf_sym = new int[num_symbols];
	int* buf_freq = new int[num_symbols];
	fread(buf_sym, sizeof(int), num_symbols, file);
	fread(buf_freq, sizeof(int), num_symbols, file);

	int num_nodes = num_symbols * 2 - 1;
	CodingNode* nodes = new CodingNode[num_nodes];
	FreqNode* freq_table = new FreqNode[num_symbols];

	for (int i = 0; i < num_symbols; i++)
	{
		freq_table[i].val = buf_sym[i];
		freq_table[i].freq = buf_freq[i];
	}
	qsort(freq_table, 0, num_symbols - 1);

	for (int i = 0; i < num_nodes; i++)
	{
		nodes[i].val = -1;
		nodes[i].freq = -1;
		nodes[i].parent = -1;
		nodes[i].left = nodes[i].right = -1;
	}

	int max_node_used = 0;
	for (int i = 0; i < num_symbols; i++)
	{
		nodes[num_symbols - 1 + i].freq = freq_table[i].freq;
		nodes[num_symbols - 1 + i].val = freq_table[i].val;
		if (nodes[num_symbols - 1 + i].freq > 0) { max_node_used = num_symbols - 1 + i; }
	}

	int cur_head = num_symbols - 2;
	while (true)
	{
		int lowest1 = coding_node_lowest_unused(nodes, cur_head + 1, max_node_used - cur_head);
		int lowest2 = coding_node_lowest_unused(nodes, cur_head + 1, max_node_used - cur_head, lowest1);

		assert(lowest1 != -1 || lowest2 != -1);

		// Last node is the head
		if (lowest1 == -1 || lowest2 == -1)
		{
			assert(lowest1 > -1 && lowest2 == -1);
			cur_head = lowest1;
			break;
		}

		assert(cur_head >= 0);
		nodes[cur_head].freq = nodes[lowest1].freq + nodes[lowest2].freq;
		nodes[lowest1].parent = nodes[lowest2].parent = cur_head;
		nodes[cur_head].left = lowest1;
		nodes[cur_head].right = lowest2;
		cur_head--;
	}

	assert(cur_head == 0);
	coding_node_assign(nodes, num_nodes);

	int io_cur_bits = 0;
	byte io_cur_byte = 0;

	int buf_pos = 0;
	int num_compressed_symbols;
	fread(&num_compressed_symbols, sizeof(int), 1, file);
	memset(buf, 0, sizeof(int) * buf_len);
	while (buf_pos < num_compressed_symbols && !feof(file))
	{
		int block_flag = bit_read(file, io_cur_byte, io_cur_bits, 1);
		int block_length = bit_read(file, io_cur_byte, io_cur_bits, RLE_BITS);
		assert(block_length > 0);
		if (block_flag)
		{
			int cur_node = 0;
			int cur_symbol = 0;
			while (cur_symbol < block_length)
			{
				int node_flag = bit_read(file, io_cur_byte, io_cur_bits, 1);
				
				if (node_flag) { cur_node = nodes[cur_node].right; }
				else { cur_node = nodes[cur_node].left; }

				assert(cur_node >= 0 && cur_node < num_nodes);

				if (nodes[cur_node].right == -1 && nodes[cur_node].left == -1)
				{
					assert(buf_pos < buf_len);
					buf[buf_pos++] = nodes[cur_node].val;
					++cur_symbol;
					cur_node = 0;
				}
			}
		}
		else
		{
			int prev_pos = buf_pos;
			assert(buf_pos < buf_len);
			buf_pos += block_length;
			assert(buf_pos <= buf_len);
		}
	}
	assert(buf_pos <= buf_len);
	delete[] nodes;
	delete[] freq_table;
}

void icdf97_2D_pass(float* buf, float* tmp, int dim, int pitch)
{
	const float a1 = 1.586134342f;
	const float a2 = 0.05298011854f;
	const float a3 = -0.8829110762f;
	const float a4 = -0.4435068522f;
	const float a5 = 1.149604398f;

	for (int r = 0; r < dim; r++)
	{
		float* s = buf + r * pitch;

		memcpy(tmp, s, sizeof(float) * dim);
		for (int i = 0; i < dim / 2; i++)
		{
			tmp[i * 2] = s[i];
			tmp[i * 2 + 1] = s[i + dim / 2];
		}
		memcpy(s, tmp, sizeof(float) * dim);

		for (int i = 0; i < dim; i++)
		{
			if (i & 1) { s[i] *= a5; }
			else { s[i] /= a5; }
		}
	
		for (int i = 2; i < dim; i += 2) { s[i] += a4 * (s[i - 1] + s[i + 1]); }
		s[0] += 2 * a4 * s[1];

		for (int i = 1; i < dim - 2; i += 2) { s[i] += a3 * (s[i - 1] + s[i + 1]); }
		s[dim - 1] += 2 * a3 * s[dim - 2];

		for (int i = 2; i < dim; i += 2) { s[i] += a2 * (s[i - 1] + s[i + 1]); }
		s[0] += 2 * a2 * s[1];

		for (int i = 1; i < dim - 2; i += 2) { s[i] += a1 * (s[i - 1] + s[i + 1]); }
		s[dim - 1] += 2 * a1 * s[dim - 2];
	}

	// transpose for repeated operation along the other axis
	memcpy(tmp, buf, sizeof(float) * pitch * pitch);
	for (int r = 0; r < dim; r++)
	{
		for (int c = 0; c < dim; c++) { buf[r * pitch + c] = tmp[c * pitch + r]; }
	}
}

void icdf97_2D_levels(float* s, float* tmp, int dim, int levels, float quantization_mult)
{
	int pitch = dim;
	for (int i = 0; i < levels - 1; i++) { dim /= 2; }

	for (int i = 0; i < levels; i++)
	{
		for (int y = 0; y < dim; y++)
		{
			for (int x = 0; x < dim; x++)
			{
				if (x >= dim / 2 || y >= dim / 2) { s[y * pitch + x] *= quantization_mult; }
			}
		}
		
		icdf97_2D_pass(s, tmp, dim, pitch);
		icdf97_2D_pass(s, tmp, dim, pitch);
		dim *= 2;
	}
}

int decode_block_wavelet(int* coefs, float* channel, int channel_width, int channel_height, int block_dim, float quantization_coef)
{
	assert(channel_width <= block_dim);
	assert(channel_height <= block_dim);
	assert(ispow2(block_dim));

	float* buf = new float[block_dim * block_dim];
	float* tmp = new float[block_dim * block_dim];

	int num_coefs = 0;
	for (int i = 0; i < block_dim * block_dim; i++)
	{
		int val = coefs[num_coefs++];
		if (val == BLOCK_END_FLAG)
		{
			for (int j = i; j < block_dim * block_dim; j++) { buf[j] = 0; }
			break;
		}
		buf[i] = val;
	}

	float log_level = logf(block_dim) / logf(2) + 1e-6f;
	icdf97_2D_levels(buf, tmp, block_dim, max(1, (int)(log_level - 3)), quantization_coef);

	for (int y = 0; y < block_dim; y++)
	{
		if (y >= channel_height) { break; }
		for (int x = 0; x < block_dim; x++)
		{
			if (x >= channel_width) { break; }
			channel[y * channel_width + x] = buf[y * block_dim + x] + 127.0f;
		}
	}

	delete[] buf;
	delete[] tmp;
	return num_coefs;
}

byte* image_load(FILE* f, int* width, int* height)
{
	if (!f) { return 0; }
	if (width) { *width = 0; }
	if (height) { *height = 0; }

	byte magic_str[image_magic_header_len];
	fread(magic_str, 1, image_magic_header_len, f);
	if (memcmp(image_magic_header, magic_str, image_magic_header_len))
	{
		printf("image_load: Invalid image file\n");
		return 0;
	}

	byte header_data[1024];
	int header_size;
	fread(&header_size, sizeof(int), 1, f);
	if (header_size > (int)sizeof(header_data))
	{
		printf("image_load: Header too large\n");
		return 0;
	}
	fread(header_data, header_size, 1, f);

	int* img_dim = (int*)(header_data);
	if (*(int*)(header_data + 8) != 2)
	{
		printf("image_load: Invalid channel specifier\n");
		return 0;
	}
	int bpp = 3;
	int chroma_downstep = *(int*)(header_data + 16);
	float wavelet_quantization_step = *(float*)(header_data + 24);

	if (width) { *width = img_dim[0]; }
	if (height) { *height = img_dim[1]; }

	int chroma_width = img_dim[0] / chroma_downstep;
	int chroma_height = img_dim[1] / chroma_downstep;
	int block_dim = npow2(max(img_dim[0], img_dim[1]));
	int chroma_block_dim = npow2(max(chroma_width, chroma_height));

	int quantized_buf_size = 0;
	fread(&quantized_buf_size, sizeof(int), 1, f);
	if (quantized_buf_size <= 0)
	{
		printf("image_load: Invalid quantized_buf_size\n");
		return 0;
	}

	int* quantized_buf = new int[quantized_buf_size];
	decode_huffman_buffer(f, quantized_buf, quantized_buf_size);

	byte* loaded_data = new byte[img_dim[0] * img_dim[1] * bpp];

	float* channel_lum = new float[img_dim[0] * img_dim[1]];
	int quantized_buf_read_pos = 0;
	quantized_buf_read_pos += decode_block_wavelet(quantized_buf + quantized_buf_read_pos, channel_lum, img_dim[0], img_dim[1], block_dim, wavelet_quantization_step);
	for (int y = 0; y < img_dim[1]; y++)
	{
		for (int x = 0; x < img_dim[0]; x++)
		{
			loaded_data[(y * img_dim[0] + x) * bpp] = min(255, max(0, channel_lum[y * img_dim[0] + x]));
		}
	}
	delete[] channel_lum;

	float* channel_u = new float[(img_dim[0] + chroma_downstep) * (img_dim[1] + chroma_downstep) / (chroma_downstep * chroma_downstep)];
	float* channel_v = new float[(img_dim[0] + chroma_downstep) * (img_dim[1] + chroma_downstep) / (chroma_downstep * chroma_downstep)];
	quantized_buf_read_pos += decode_block_wavelet(quantized_buf + quantized_buf_read_pos, channel_u, chroma_width, chroma_height, chroma_block_dim, wavelet_quantization_step);
	quantized_buf_read_pos += decode_block_wavelet(quantized_buf + quantized_buf_read_pos, channel_v, chroma_width, chroma_height, chroma_block_dim, wavelet_quantization_step);
	for (int y = 0; y < img_dim[1]; y += chroma_downstep)
	{
		for (int x = 0; x < img_dim[0]; x += chroma_downstep)
		{
			float u = channel_u[(y / chroma_downstep) * img_dim[0] / chroma_downstep + (x / chroma_downstep)];
			float v = channel_v[(y / chroma_downstep) * img_dim[0] / chroma_downstep + (x / chroma_downstep)];
			u = min(255, max(0, u));
			v = min(255, max(0, v));
			for (int i = 0; i < chroma_downstep; i++)
			{
				for (int j = 0; j < chroma_downstep; j++)
				{
					loaded_data[((y + i) * img_dim[0] + x + j) * bpp + 1] = u;
					loaded_data[((y + i) * img_dim[0] + x + j) * bpp + 2] = v;
				}
			}
		}
	}

	delete[] channel_u;
	delete[] channel_v;
	delete[] quantized_buf;

	for (int idx = 0; idx < img_dim[0] * img_dim[1]; idx++)
	{
		float y = loaded_data[idx * 3] / 255.0f;
		float cb = loaded_data[idx * 3 + 1] / 255.0f;
		float cr = loaded_data[idx * 3 + 2] / 255.0f;
		float r = -0.701f + 1.402f * cr + y;
		float g = 0.529136f - 0.344136f * cb - 0.714136f * cr + y;
		float b = -0.886f + 1.772f * cb + y;

		loaded_data[idx * 3] = (byte)min(255, max(0, r * 255.0f));
		loaded_data[idx * 3 + 1] = (byte)min(255, max(0, g * 255.0f));
		loaded_data[idx * 3 + 2] = (byte)min(255, max(0, b * 255.0f));
	}

	return loaded_data;
}
