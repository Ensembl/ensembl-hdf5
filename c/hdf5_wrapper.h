// Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//      http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
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

hid_t create_file(char * filename, hsize_t rank, char ** dim_names, hsize_t * dim_sizes, hsize_t * dim_label_lengths, hsize_t * chunk_sizes);
void store_dim_labels(hid_t file, char * dim_name, hsize_t dim_size, char ** dim_labels);
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
void set_hdf5_log(int value);
#endif
