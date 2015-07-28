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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5_wrapper.h"
#include "hdf5_wrapper_priv.h"

int main(int argc, char ** argv) {
	int rank = 2;
	char * dim_names[] = {"X", "Y"};
	hsize_t dim_sizes[] = {3, 3};	
	char * xlabels[] = {"Hi","Lo","Middle"};
	char * ylabels[] = {"Left","Right","Center"};
	char ** dim_labels[] = {xlabels, ylabels};

	printf("Resetting big dim cutoff\n");
	set_big_dim_length(0);
	set_hdf5_log(true);

	printf("Creating files\n");
	hid_t file = create_file("TEST.hd5", rank, dim_names, dim_sizes, dim_labels, NULL);

	if (get_file_rank(file) != 2)
		abort();

	if (get_file_core_rank(file) != 2)
		abort();

	puts("Testing dim names");
	StringArray * names = get_dim_names(file);
	if (names->count != 2)
		abort();

	if (strcmp(get_string_in_array(names, 0), "X") && strcmp(get_string_in_array(names, 0), "Y"))
		abort();

	if (strcmp(get_string_in_array(names, 1), "X") && strcmp(get_string_in_array(names, 1), "Y"))
		abort();

	puts("Testing dim labels");
	StringArray * labels;
	if (strcmp(get_string_in_array(names, 0), "X") == 0) {
		labels = get_all_dim_labels(file, 0);
	} else {
		labels = get_all_dim_labels(file, 1);
	}

	puts("Testing dim labels");
	if (labels->count != 3)
		abort();
	if (strcmp(get_string_in_array(labels, 0), "Hi"))
		abort();
	if (strcmp(get_string_in_array(labels, 1), "Lo"))
		abort();
	if (strcmp(get_string_in_array(labels, 2), "Middle"))
		abort();

	hsize_t coord[] = {0,0};
	hsize_t coord2[] = {1,1};
	hsize_t * coord_array[] = {coord, coord2};
	double values[] = {1, 2};
	printf("Storing values\n");
	store_values(file, 2, coord_array, values);

	bool set_dims[] = {1, 0};
	hsize_t constraints[] = {0, 0};

	printf("Fetching values\n");
	StringResultTable * res = fetch_string_values(file, set_dims, constraints);

	printf("Testing output:\n");

	printf("- Column count: %lli\n", res->columns);
	if (res->columns != 1)
		abort();

	printf("- Variable dims:\n");
	if (strcmp(res->dims[0], "Y"))
		abort();

	printf("- Data points: %lli\n", res->rows);
	if (res->rows != 1)
		abort();

	printf("- Value:\n");
	if (res->values[0] != 1)
		abort();

	bool set_dims2[] = {0, 0};

	printf("Fetching values again\n");
	StringResultTable * res2 = fetch_string_values(file, set_dims2, constraints);

	printf("- Column count: %lli\n", res2->columns);
	if (res2->columns != 2)
		abort();

	printf("- Data points: %lli\n", res2->rows);
	if (res2->rows != 2)
		abort();

	printf("Success\n");
	return 0;
}
