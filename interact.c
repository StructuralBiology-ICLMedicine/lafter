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

int get_num_jobs(void){
  // Obtain thread number from environmental variables
  char* thread_number = getenv("OMP_NUM_THREADS");
  int nthreads = 0;
  if (thread_number){
    // If thread number specified by user - use this one
    nthreads = atoi(thread_number);
  }
  if (nthreads < 1){
    // If thread number still not set - try sysconf
    nthreads = sysconf(_SC_NPROCESSORS_ONLN);
  }
  if (nthreads < 1){
    // If variables are both empty - use a single thread
    nthreads = 1;
  }
  return nthreads;
}

arguments *parse_args(int argc, char **argv){
  // Print usage and disclaimer
  if (argc < 7){
    printf("\n    Usage: %s --v1 half_map1.mrc --v2 half_map2.mrc ( --mask mask.mrc || --rad voxel_radius )[ --sharp ][ --fsc cut_off ][ --correct_overfitting ]\n", argv[0]);
  }

  printf("\n\n\n");
  printf("                 _______             _______      _____________  ___________________  ____________  ______________         \n");
  printf("                |\\      \\         _-         --__|\\            \\|\\                  \\|\\           \\|\\             --_      \n");
  printf("                | \\      \\       |       __      \\ \\       _____\\ \\______       _____\\ \\      _____\\ \\       __      \\     \n");
  printf("                \\  \\      \\      |\\      \\_\\      \\ \\         \\ | |     |\\      \\    |  \\         \\|  \\      \\_\\      |    \n");
  printf("                 \\  \\      \\_____| \\       __      \\ \\       __\\|\\|_____| \\      \\___|\\  \\      ___\\__ \\       __     \\_   \n");
  printf("                  \\  \\            \\ \\      \\ \\      \\ \\      \\ |        \\  \\      \\    \\  \\           \\ \\      \\ -_     -_ \n");
  printf("                   \\  \\____________\\ \\______\\ \\______\\ \\______\\|         \\  \\______\\    \\  \\___________\\ \\______\\  \\______\\\n");
  printf("                    \\ |            | |      | |      | |      |           \\ |      |     \\ |           | |      |\\ |      |\n");
  printf("                     \\|____________|\\|______|\\|______|\\|______|            \\|______|      \\|___________|\\|______| \\|______|\n\n\n\n");

  printf("    PLEASE NOTE: LAFTER absolutely requires the unfiltered halfmaps and mask used for processing - anything else is invalid\n");
  printf("                 The output maps from LAFTER are incompatible with atomic coordinate refinement but will aid model building\n");
  printf("                 Junk in = Junk out is one thing we will guarantee. Report any bugs to c.aylett@imperial.ac.uk - good luck!\n\n");
  printf("    LAFTER v1.0: Noise suppression and SNR filtering - header correct on 14-06-2018 - GNU licensed - K Ramlaul & CHS Aylett\n\n");

  // Capture user requested settings
  int i;
  arguments *args = malloc(sizeof(arguments));
  memset(args, 0, sizeof(arguments));
  args->rad = 0.000;
  args->cut = 0.143;
  args->sharp = 0.0;
  args->ovfit = 1.0;
  for (i = 1; i < argc; i++){
    if (!strcmp(argv[i], "--v1") && ((i + 1) < argc)){
      args->vol1 = argv[i + 1];
    } else if (!strcmp(argv[i], "--v2") && ((i + 1) < argc)){
      args->vol2 = argv[i + 1];
    } else if (!strcmp(argv[i], "--mask") && ((i + 1) < argc)){
      args->mask = argv[i + 1];
    } else if (!strcmp(argv[i], "--correct_overfitting")){
      args->ovfit = 0.0;
      printf("    LAFTER will attempt to correct overfitting. This can be used if residual noise is visible, but we recommend re-refining \n\n");
    } else if (!strcmp(argv[i], "--fsc") && ((i + 1) < argc)){
      args->cut = atof(argv[i + 1]);
#ifndef DEBUG
      if (args->cut < 0.1 || args->cut > 1.0){
	args->cut = 0.143;
      }
#endif
      printf("    LAFTER will use this figure as the FSC cut-off at which no new resolution shells are incorporated. The default is 0.143 \n\n");
    } else if (!strcmp(argv[i], "--rad") && ((i + 1) < argc)){
      args->rad = atof(argv[i + 1]);
      if (args->rad < 10.0){
        printf("    Necessary maps not found or not specified LAFTER requires both volumes and mask - set --rad when refined without a mask \n\n");
        printf("    Radius was not specified correctly - a float or integer value in voxels is required. Divide by Å/pixel if mask set in Å \n\n");
        exit(1);
      }
      args->mask = argv[i];
      printf("    LAFTER will soft mask the map at this radius in voxels. Please note that if refinement was masked you MUST use the mask \n\n");
    } else if (!strcmp(argv[i], "--sharp")){
      args->sharp = 1.0;
      printf("    LAFTER will sharpen the output map after the last cycle. This will not affect the FSC or any of the statistics provided \n\n");
    }
  }
  if (args->vol1 == NULL || args->vol2 == NULL || args->mask == NULL){
    printf("    Necessary maps not found or not specified LAFTER requires both volumes and mask - set --rad when refined without a mask \n\n");
    exit(1);
  }
  return args;
}

// Extend list by one using p-val
list *extend_list(list *current, double p){
  list *node = calloc(1, sizeof(list));
  node->res = current->res + current->stp;
  node->stp = p * (node->res / 8.0);
  node->prv = current;
  node->nxt = NULL;
  current->nxt = node;
  return node;
}

// End list for overfitting calculation
list *end_list(list *current){
  list *node = calloc(1, sizeof(list));
  node->res = current->res + current->stp;
  node->stp = 0.475 - (current->res + current->stp);
  if (node->stp < 0.0625){
    node->stp = 0.0625;
  }
  node->prv = current;
  node->nxt = NULL;
  current->nxt = node;
  return node;
}

// Read map header and data and return corresponding data structure
r_mrc *read_mrc(char *filename){

  FILE *f;
  int i;
  f = fopen(filename, "rb");
  if (!f){
    printf("\n\tError reading %s - bad file handle\n\n", filename);
    exit(1);
  }
  r_mrc *header = calloc(1, sizeof(r_mrc));
  fread(&header->n_crs, 4, 3, f);
  fread(&header->mode, 4, 1, f);
  /* We only accept mode 2 - c float32 / FORTRAN real - because implementing
     checks for other data types would require pointer casting everywhere */
  if (header->mode != 2){
    printf("Error reading %s - not 32 bit data\n", filename);
    exit(1);
  }
  fread(&header->start_crs,  4, 3,   f);
  fread(&header->n_xyz,      4, 3,   f);
  fread(&header->length_xyz, 4, 3,   f);
  fread(&header->angle_xyz,  4, 3,   f);
  fread(&header->map_crs,    4, 3,   f);
  fread(&header->d_min,      4, 1,   f);
  fread(&header->d_max,      4, 1,   f);
  fread(&header->d_mean,     4, 1,   f);
  fread(&header->ispg,       4, 1,   f);
  fread(&header->nsymbt,     4, 1,   f);
  fread(&header->extra,      4, 25,  f);
  fread(&header->ori_xyz,    4, 3,   f);
  fread(&header->map,        1, 4,   f);
  fread(&header->machst,     1, 4,   f);
  fread(&header->rms,        4, 1,   f);
  fread(&header->nlabl,      4, 1,   f);
  fread(&header->label,      1, 800, f);
  // Check cube
  if (header->n_crs[0] != header->n_crs[1] || header->n_crs[0] != header->n_crs[2]){
    printf("Error reading %s - map is not a cube\n", filename);
    exit(1);
  }
  /* Assign float array for data, initialize it to zero and read it in - note that the
     endianness is not corrected - architectures may therefore be cross-incompatible*/
  header->data = calloc((header->n_crs[0] * header->n_crs[1] * header->n_crs[2]), sizeof(float));
  if (!header->data){
    printf("Error reading %s - map not allocated\n", filename);
    exit(1);
  }
  fread(header->data, sizeof(float), (header->n_crs[0] * header->n_crs[1] * header->n_crs[2]), f);
  fclose(f);
  return header;
}

// Write MRC file given an mrc structure and corresponding data
void write_mrc(r_mrc* header, double *vol, char* filename, int32_t size){

  int i;
  double total    = size * size * size;
  double current  = 0;
  double tmp, sum = 0;

  float *hold = header->data;
  header->data = calloc(total, sizeof(float));
  if (!header->data){
    printf("Error reading %s - map not allocated\n", filename);
    exit(1);
  }
  for (i = 0; i < total; i++){
    header->data[i] = (float) vol[i];
  }

  // Calculate new min, max and mean figures for header
  header->d_min = header->data[0];
  header->d_max = header->data[0]; 
  for (i = 0; i < total; i++){
    current = (double) header->data[i];
    if (current < header->d_min){
      header->d_min = (float) current;
    }
    if (current > header->d_max){
      header->d_min = (float) current;
    }
    sum += current;
  }
  header->d_mean = (float) (sum / total);

  // Calculate RMSD to fill in the rms field for scaling
  sum = 0;
  for (i = 0; i < total; i++){
    tmp  = (double) (header->d_mean - header->data[i]);
    sum += tmp * tmp;
  }
  header->rms = (float) sqrt((sum / total));

  // Write out 1024 byte header
  FILE * f;
  f = fopen(filename, "wb");
  if (!f){
    printf("Error writing %s - bad file handle\n", filename);
    exit(1);
  }
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&header->mode,       4, 1,   f);
  fwrite(&header->start_crs,  4, 3,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&size,               4, 1,   f);
  fwrite(&header->length_xyz, 4, 3,   f);
  fwrite(&header->angle_xyz,  4, 3,   f);
  fwrite(&header->map_crs,    4, 3,   f);
  fwrite(&header->d_min,      4, 1,   f);
  fwrite(&header->d_max,      4, 1,   f);
  fwrite(&header->d_mean,     4, 1,   f);
  fwrite(&header->ispg,       4, 1,   f);
  fwrite(&header->nsymbt,     4, 1,   f);
  fwrite(&header->extra,      4, 25,  f);
  fwrite(&header->ori_xyz,    4, 3,   f);
  fwrite(&header->map,        1, 4,   f);
  fwrite(&header->machst,     1, 4,   f);
  fwrite(&header->rms,        4, 1,   f);
  fwrite(&header->nlabl,      4, 1,   f);
  fwrite(&header->label,      1, 800, f);

  // Write data to file
  fwrite(header->data, 4, total, f);
  fclose(f);

  free(header->data);
  header->data = hold;
  return;
}

// Make mask from radius in voxels
r_mrc *make_msk(r_mrc *in, arguments *args){
  r_mrc *out = malloc(sizeof(r_mrc));
  int32_t size = in->n_crs[0];
  int32_t ind = 0;
  args->rad *= args->rad;
  double i, j, k, cen = (double) size / 2;
  double norm;
  memcpy(out, in, sizeof(r_mrc));
  out->data = calloc(size * size * size, sizeof(float));
  for(int32_t _k = 0; _k < size; _k++){
    k = (double) _k - cen;
    k = k * k;
    for(int32_t _j = 0; _j < size; _j++){
      j = (double) _j - cen;
      j = j * j;
      for(int32_t _i = 0; _i < size; _i++){
	i = (double) _i - cen;
	i = i * i;
	norm = (double) k + j + i;
	out->data[ind++] = 1.0 / (1.0 + pow((norm / args->rad), 8.0));
      }
    }
  }
  return out;
}

// Add MRC map in to out
void add_map(r_mrc *in, double *out){
  int32_t i, max = in->n_crs[0] * in->n_crs[1] * in->n_crs[2];
  for (i = 0; i < max; i++){
    out[i] += (double) in->data[i];
  }
  return;
}

// Multiply out by in elementwise
void apply_mask(r_mrc *in, double *out){
  int32_t i, max = in->n_crs[0] * in->n_crs[1] * in->n_crs[2];
  for (i = 0; i < max; i++){
    out[i] *= (double) in->data[i];
  }
  return;
}
