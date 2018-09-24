/*                                                                       
 * Copyright 10/12/2017 - Dr. Christopher H. S. Aylett                   
 *                                                                       
 * This program is free software; you can redistribute it and/or modify  
 * it under the terms of version 3 of the GNU General Public License as  
 * published by the Free Software Foundation.                            
 *                                                                       
 * This program is distributed in the hope that it will be useful,       
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
 * GNU General Public License for more details - YOU HAVE BEEN WARNED!   
 *                                                                       
 * Program: LAFTER V1.0                                                  
 *                                                                       
 * Authors: Chris Aylett                                                 
 *                                                                       
 */

// Library header inclusion for linking                                  
#include "lafter.h"

// Definitions
#define THREAD 0
#define STEP   1

// FFT determine voxel
#define FFT_IJK(I, J, K, SIZE) \
  (((K < 0) ? (K + SIZE) : (K)) * SIZE * (SIZE / 2 + 1) + ((J < 0) ? (J + SIZE) : J) * (SIZE / 2 + 1) + (I))


// Soften map
void soften_map(double *in, double *out, int32_t size){
  int32_t i, j, k;
  int32_t neighbours[6] = {1, -1, size, -size, size * size, -size * size};
  int32_t sqr = size * size;
  int32_t cub = size * size * size;
  double mean;
  double *tmp;
  for (k = 0; k < 8; k++){
    for (i = 0; i < cub; i++){
      if (in[i] != 0.0){
	out[i] = in[i];
	continue;
      }
      mean = 0;
      k = 0;
      for (j = 0; j < 6; j++){
	mean += in[(cub + i + neighbours[j]) % (cub)];
      }
      out[i] = mean / 6.0;
    }
    tmp = in;
    in = out;
    out = tmp;
  }
  return;
}

// Write out upsampled mrc file combining in1/2 - sharpening specified by command line args
void write_upsampled(fftw_complex *in, double res, r_mrc *mask, char *name, arguments *args, int32_t size){

  // Allocate memory for maps
  fftw_complex *ks = fftw_alloc_complex((size * 2) * (size * 2) * (size + 1));
  double *rs = fftw_alloc_real((size * 2) * (size * 2) * (size * 2));

  // Make FFTW plans
  fftw_plan fft_ks_rs = fftw_plan_dft_c2r_3d(size * 2, size * 2, size * 2, ks, rs, FFTW_ESTIMATE);

  // Zero fill maps
  memset(ks, 0, (size * 2) * (size * 2) * (size + 1) * sizeof(fftw_complex));
  memset(rs, 0, (size * 2) * (size * 2) * (size * 2) * sizeof(double));

  // Butterworth lowpass with sharpening from in to out
  double hires = res * res;
  double dim = (double) size;
  double norms, kd, jd, id;
  int32_t full = size * 2;
  int32_t half = size + 1;
  int32_t cut1 = size / 2;
  int32_t cut2 = full - cut1;
  int32_t index = 0;
  for(int32_t k = 0; k < full; k++){
    if(k > cut1 && k <= cut2){
      continue;
    }
    kd = ((double) ((k > size) ? k - full: k)) / dim;
    kd = kd * kd;
    for(int32_t j = 0; j < full; j++){
      if(j > cut1 && j <= cut2){
	continue;
      }
      jd = ((double) ((j > size) ? j - full: j)) / dim;
      jd = jd * jd;
      for(int i = THREAD; i < half; i += STEP){
	if (i > cut1){
	  continue;
	}
	id = ((double) i) / dim;
	id = id * id;
	norms = kd + jd + id;
	ks[(k * full * half + j * half + i)] = in[index++] * (pow(2.0, (8.0 * norms * args->sharp)) / (1.0 + pow((norms / hires), 8.0)));
      }
    }
  }
  fftw_execute(fft_ks_rs);
  write_mrc(mask, rs, name, full);
  fftw_free(ks);
  fftw_free(rs);
  return;
}
