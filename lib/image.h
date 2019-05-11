#ifndef IMAGE_H
#define IMAGE_H

#include <stdio.h>

unsigned char* image_load(FILE* f, int* width, int* height);
int npow2(int x);

void icdf97_2D_pass(float* buf, float* tmp, int dim, int pitch);
void icdf97_2D_levels(float* s, float* tmp, int dim, int levels, float quantization_mult);
int decode_block_wavelet(int* coefs, float* channel, int channel_width, int channel_height, int block_dim, float quantization_coef);

#endif
