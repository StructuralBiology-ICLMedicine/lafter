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
#define STEP 1

// Add FFT in to FFT out
void add_fft(fftw_complex *in, fftw_complex *out, int32_t full){
  int32_t size = full * full * (full / 2 + 1);
  for(int index = THREAD; index < size; index += STEP){
    out[index] += in[index];
  }
  return;
}

// Calculate FSC over map
double calc_fsc(fftw_complex *half1, fftw_complex *half2, int32_t full){
  long double numerator = 0.0;
  long double denomin_1 = 0.0;
  long double denomin_2 = 0.0;
  int32_t size = full * full * (full / 2 + 1);
  for(int index = THREAD; index < size; index += STEP){
    numerator += creal(half1[index] * conj(half2[index]));
    denomin_1 += creal(half1[index] * conj(half1[index]));
    denomin_2 += creal(half2[index] * conj(half2[index]));
  }
  return (double) (numerator / sqrtl(denomin_1 * denomin_2));
}

// Apply bandpass to in and writes to out
void bandpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full){
  double hires = node->res + node->stp;
  double lores = node->res;
  hires = hires * hires;
  lores = lores * lores;
  double dim = (double) full;
  double norms, kd, jd, id;
  int32_t size = (full / 2) + 1;
  int32_t index = 0;
  for(int _k = 0, k = 0; _k < full; _k++, k = (_k < size) ? _k : _k - full){
    kd = ((double) k) / dim;
    for(int _j = 0, j = 0; _j < full; _j++, j = (_j < size) ? _j : _j - full){
      jd = ((double) j) / dim;
      for(int _i = THREAD, i = THREAD; _i < size; _i += STEP, i = _i){
	id = ((double) i) / dim;
	norms = kd * kd + jd * jd + id * id;
	out[index] = in[index] * ((1.0 / (1.0 + pow((norms / hires), 8.0))) - (1.0 / (1.0 + pow((norms / lores), 8.0))));
	index++;
      }
    }
  }
  return;
}

// Butterworth lowpass from in to out
void lowpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t full){
  double hires = node->res + node->stp;
  hires = hires * hires;
  double dim = (double) full;
  double norms, kd, jd, id;
  int32_t size = (full / 2) + 1;
  int32_t index = 0;
  for(int _k = 0, k = 0; _k < full; _k++, k = (_k < size) ? _k : _k - full){
    kd = ((double) k) / dim;
    for(int _j = 0, j = 0; _j < full; _j++, j = (_j < size) ? _j : _j - full){
      jd = ((double) j) / dim;
      for(int _i = THREAD, i = THREAD; _i < size; _i += STEP, i = _i){
	id = ((double) i) / dim;
	norms = kd * kd + jd * jd + id * id;
	out[index] = in[index] * (1.0 / (1.0 + pow((norms / hires), 8.0)));
	index++;
      }
    }
  }
  return;
}
