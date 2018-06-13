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

#define THREAD 0
#define STEP 1

// Suppress noise between in/out
double suppress_noise(double *in1, double *in2, double *out1, double *out2, r_mrc *mask, list *node, int32_t size){
  int32_t i, max = size * size * size;
  double snr, rmsd, cur, cor, count = 0.0;
  long double p = 0.0;
  long double noise = 0.0;
  long double power = 0.0;
  long double norml = 0.0;

  // Calculate mean noise and mean signal
  for (i = THREAD; i < max; i += STEP){
    if (mask->data[i] < 0.99){
      continue;
    }
    count += 1.0;
    cur = in1[i] - in2[i];
    noise += cur * cur;
    cur = in1[i] + in2[i];
    power += cur * cur;
    norml += in1[i] * in1[i] + in2[i] * in2[i];
  }
  noise /= count;
  power /= count;
  snr  = fabsl(1.0 - noise / power);
  norml /= count * 2.0;
  rmsd = sqrtl(norml);
  // Correct according to probability and power
  for (i = THREAD; i < max; i += STEP){
    cor = in1[i] + in2[i];
    cur = snr * (1.0 - exp(-0.5 * ((cor * cor) / noise)));
    out1[i] += (in1[i] / rmsd) * 20.0 * node->stp * cur;
    out2[i] += (in2[i] / rmsd) * 20.0 * node->stp * cur;
    if(mask->data[i] > 0.99){
      p += cur;
    }
  }
  p /= count;
  return (double) p;
}
