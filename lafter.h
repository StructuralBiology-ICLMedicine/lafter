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
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>
#include <pthread.h>
#include <unistd.h>
#include <complex.h>
#include <fftw3.h>

/* Type definitions */

// #define DEBUG

// Arguments
typedef struct {
  char   *vol1;
  char   *vol2;
  char   *mask;
  double sharp;
  double ovfit;
  double   rad;
  double   cut;
} arguments;

// List node
typedef struct list list;
struct list {
  double res;
  double stp;
  double crf;
  double fsc;
  double max;
  list  *prv;
  list  *nxt;
};

// MRC files
typedef struct {
  int32_t n_crs[3];
  int32_t mode;
  int32_t start_crs[3];
  int32_t n_xyz[3];
  float   length_xyz[3];
  float   angle_xyz[3];
  int32_t map_crs[3];
  float   d_min;
  float   d_max;
  float   d_mean;
  int32_t ispg;
  int32_t nsymbt;
  int32_t extra[25];
  int32_t ori_xyz[3];
  char    map[4];
  char    machst[4];
  float   rms;
  int32_t nlabl;
  char    label[800];
  float  *data;
} r_mrc;


/* Function definitions */

arguments *parse_args(int argc, char **argv);
// Read and return arguments structure
// Exit if required args not specified

int get_num_jobs();
// Returns number of processors

list *extend_list(list *node, double p);
// Extend list by one using p-val
// Calculates step size

list *end_list(list *node);
// Finish list to 0.5 for overfit calculation

r_mrc *read_mrc(char* filename);
// Read mrc file and build struct

void write_mrc(r_mrc *mrc, double *map,char* filename, int32_t size);
// Read mrc file and build struct

r_mrc *make_msk(r_mrc *mrc, arguments *args);
// Make mask from radius in voxels

void add_map(r_mrc *in, double *out);
// Add MRC map in to out

void add_fft(fftw_complex *in, fftw_complex *out, int32_t size);
// Add FFT in to out

void apply_mask(r_mrc *in, double *out);
// Multiply out by in elementwise

void bandpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t size);
// Apply bandpass to in and writes to out
// List node specifies resolutions

void lowpass_filter(fftw_complex *in, fftw_complex *out, list *node, int32_t size);
// Butterworth lowpass from in to out
// List node specifies resolution

double calc_fsc(fftw_complex *half1, fftw_complex *half2, int32_t size);
// Calculate FSC over map
// Returns FSC

double suppress_noise(double *in1, double *in2, double *out1, double *out2, r_mrc *mask, list *node, int32_t size);
// Suppress noise between in/out
// Returns mean p-val in mask

double truncate_map(double *in1, double *in2, double *out, r_mrc *mask, list *node, arguments *args, int32_t size);
// Updates out if in1/2 over noise
// Returns fractional recovery

void soften_map(double *in, double *out, int32_t size);
// Soften edges of map

void write_upsampled(fftw_complex *in, double max_res, r_mrc *mask, char *name, arguments *args, int32_t size);
// Write out upsampled mrc file combining in1/2
// Sharpening specified by command line
