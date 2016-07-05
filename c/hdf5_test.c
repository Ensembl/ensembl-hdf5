// Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
// Copyright [2016] EMBL-European Bioinformatics Institute
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5_wrapper.h"
#include "hdf5_wrapper_priv.h"

int main(int argc, char ** argv) {
	int rank = 2;
	char * dim_names[] = {"snp", "gene"};
	hsize_t dim_sizes[] = {2, 2};
	char * xlabels[] = {"A","B"};
	char * ylabels[] = {"rs1"};
	char * ylabels2[] = {"rs2"};
	hsize_t dim_label_lengths[] = {3, 1};

	printf("Resetting big dim cutoff\n");
	set_big_dim_length(1);
	set_hdf5_log(2);

	printf("Creating files\n");
	hid_t file = create_file("TEST.hd5", rank, dim_names, dim_sizes, dim_label_lengths, NULL);

	store_dim_labels(file, "snp", 1, ylabels);

	close_file(file);

	file = open_file("TEST.hd5");
	store_dim_labels(file, "gene", 2, xlabels);
	store_dim_labels(file, "snp", 1, ylabels2);

	if (get_file_rank(file) != 2)
		abort();

	if (get_file_core_rank(file) != 2)
		abort();

	puts("Testing dim names");
	StringArray * names = get_dim_names(file);
	if (names->count != 2)
		abort();

	if (strcmp(get_string_in_array(names, 0), "gene") && strcmp(get_string_in_array(names, 0), "snp"))
		abort();

	if (strcmp(get_string_in_array(names, 1), "gene") && strcmp(get_string_in_array(names, 1), "snp"))
		abort();

	puts("Testing dim labels");
	StringArray * labelsX, * labelsY;
	if (strcmp(get_string_in_array(names, 0), "gene") == 0) {
		labelsX = get_all_dim_labels(file, 0);
		labelsY = get_all_dim_labels(file, 1);
	} else {
		labelsY = get_all_dim_labels(file, 0);
		labelsX = get_all_dim_labels(file, 1);
	}
	destroy_string_array(names);

	puts("Testing dim labels");
	if (labelsX->count != 2)
		abort();
	if (strcmp(get_string_in_array(labelsX, 0), "A"))
		abort();
	if (strcmp(get_string_in_array(labelsX, 1), "B"))
		abort();
	destroy_string_array(labelsX);

	if (labelsY->count != 2)
		abort();
	if (strcmp(get_string_in_array(labelsY, 0), "rs1"))
		abort();
	if (strcmp(get_string_in_array(labelsY, 1), "rs2"))
		abort();
	destroy_string_array(labelsY);

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
	if (strcmp(res->dims[0], "gene"))
		abort();

	printf("- Data points: %lli\n", res->rows);
	if (res->rows != 1)
		abort();

	printf("- Value:\n");
	if (res->values[0] != 1)
		abort();

	destroy_string_result_table(res);

	bool set_dims2[] = {0, 0};

	printf("Fetching values again\n");
	StringResultTable * res2 = fetch_string_values(file, set_dims2, constraints);

	printf("- Column count: %lli\n", res2->columns);
	if (res2->columns != 2)
		abort();

	printf("- Data points: %lli\n", res2->rows);
	if (res2->rows != 2)
		abort();

	destroy_string_result_table(res2);

	bool set_dims3[] = {1, 1};

	printf("Fetching values single\n");
	StringResultTable * res3 = fetch_string_values(file, set_dims3, constraints);

	if (res3->columns)
		abort();

	destroy_string_result_table(res3);

	printf("Success\n");
	return 0;
}
