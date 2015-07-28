// LICENSE
//
//  Copyright (c) 1999-2015 The European Bioinformatics Institute and
//  Genome Research Limited.  All rights reserved.
//
//  This software is distributed under a modified Apache license.
//  For license details, please see
//
// http://www.ensembl.org/info/about/code_licence.html
//
#ifndef _HDF5_WRAPPER_H_
#define _HDF5_WRAPPER_H_

#ifndef bool
#define bool char
#define true 1
#define false 0
#endif

#include "hdf5.h"

typedef struct string_array_st {
	char * array;
	hsize_t length;
	hsize_t count;
} StringArray;

typedef struct result_table_st {
	hsize_t rows, columns;
	hsize_t * dims;
	hsize_t ** coords;
	double * values;
} ResultTable;

typedef struct string_result_table_st {
	hsize_t rows, columns;
	char ** dims;
	char *** coords;
	double * values;
	StringArray ** dim_labels;
	StringArray * dim_names;
} StringResultTable;

hid_t create_file(char * filename, hsize_t rank, char ** dim_names, hsize_t * dim_sizes, char *** dim_labels, hsize_t * chunk_sizes);
void store_values(hid_t file, hsize_t count, hsize_t ** coords, double * values);
hid_t open_file(char * filename);
StringResultTable * fetch_string_values(hid_t file, bool * set_dims, hsize_t * constraints);
void destroy_string_result_table(StringResultTable * table);
void close_file(hid_t file);

hsize_t get_file_core_rank(hid_t file);
hsize_t get_file_rank(hid_t file);
StringArray * get_dim_names(hid_t file);
StringArray * get_all_dim_labels(hid_t file, hsize_t dim);
char * get_string_in_array(StringArray * sarray, hsize_t index);
void destroy_string_array(StringArray * sarray);
#endif
