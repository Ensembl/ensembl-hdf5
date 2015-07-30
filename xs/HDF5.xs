// LICENSE
//  Copyright (c) 1999-2015 The European Bioinformatics Institute and
//  Genome Research Limited.  All rights reserved.
//
//  This software is distributed under a modified Apache license.
//  For license details, please see
//
// http://www.ensembl.org/info/about/code_licence.html
//
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "hdf5_wrapper.h"

#include "ppport.h"

struct hdf5_file_st {
	hid_t file;
	HV * dim_indices;
	int * dim_name_lengths;
};

MODULE = Bio::EnsEMBL::HDF5 PACKAGE = Bio::EnsEMBL::HDF5 PREFIX=hdf5_

void
hdf5_create(filename_sv, dim_labels_hv)
		SV * filename_sv;
		HV * dim_labels_hv;
	PREINIT:
		hsize_t rank;
		char ** dim_names;
		hsize_t * dim_sizes;
		char *** dim_labels;
		int dim;
		hsize_t index;
		char * filename;
	CODE:
		// Allocate memory
		rank = hv_iterinit(dim_labels_hv);
		dim_names = calloc(rank, sizeof(char*));
		dim_sizes = calloc(rank, sizeof(hsize_t));
		dim_labels = calloc(rank, sizeof(char **));

		// Read dimension names and labels from hash ref
		for (dim = 0; dim < rank; dim++) {
			// Iteration
			HE * hash_entry = hv_iternext(dim_labels_hv);

			// Dereference data
			SV * sv_key = hv_iterkeysv(hash_entry);
			SV * labels_sv = hv_iterval(dim_labels_hv, hash_entry);
			AV * labels_av = (AV *) SvRV(labels_sv);

			// Update C arrays
			dim_names[dim] = SvPV(sv_key, PL_na);
			dim_sizes[dim] = av_len(labels_av) + 1;
			dim_labels[dim] = calloc(dim_sizes[dim], sizeof(char **));
			for (index = 0; index < dim_sizes[dim]; index++) {
				SV ** dim_label_sv = av_fetch(labels_av, index, 0); 
				dim_labels[dim][index] = SvPV(*dim_label_sv, PL_na);
			}
		}

		// Create file
		filename = SvPV_nolen(filename_sv);
		create_file(filename, rank, dim_names, dim_sizes, dim_labels, NULL);

		// Clean up memory
		for (dim = 0; dim < rank; dim++) {
			free(dim_labels[dim]);
		}
		free(dim_names);
		free(dim_sizes);
		free(dim_labels);

void * 
hdf5_open(filename_sv)
		SV * filename_sv
	PREINIT:
		struct hdf5_file_st * file;
		char * filename;
		StringArray * dim_names_sa;
		int index, dim, rank;
	CODE:
		// Create struct 
		file = calloc(1, sizeof(struct hdf5_file_st));

		// Open file
		filename = SvPV_nolen(filename_sv);
		file->file = open_file(filename);

		// Allocate storage
		rank = get_file_rank(file->file);
		file->dim_indices = newHV();
		file->dim_name_lengths = calloc(rank, sizeof(int));

		// Read dimension indices from the file store in HV
		dim_names_sa = get_dim_names(file->file);
		for (dim = 0; dim < rank; dim++) {
			char * name = get_string_in_array(dim_names_sa, dim);
			// Store dimension index
			hv_store(file->dim_indices, name, strlen(name), newSViv(dim), 0);
			// Store length of dimension name
			file->dim_name_lengths[dim] = strlen(name);
		}
		destroy_string_array(dim_names_sa);

		RETVAL = file; 
	OUTPUT:
		RETVAL

SV *
hdf5_get_dim_labels(file)
		void * file
	PREINIT:
		struct hdf5_file_st * file_st = (struct hdf5_file_st *) file;
		HV * dim_labels_hv;
		int dim, rank;
		hid_t index;
		StringArray * dim_names_sa;
	CODE:
		dim_labels_hv = newHV();
		dim_names_sa = get_dim_names(file_st->file);
		rank = get_file_rank(file_st->file);
		for (dim = 0; dim < rank; dim++) {
			// Read labels
			StringArray * dim_labels_sa = get_all_dim_labels(file_st->file, dim);

			// Store in array ref
			AV * dim_labels_av = newAV();
			for (index = 0; index < dim_labels_sa->count; index++) {
				char * dim_label = get_string_in_array(dim_labels_sa, index);
				SV * dim_label_sv = newSVpv(dim_label, 0);
				av_push(dim_labels_av, dim_label_sv);
			}

			// Store array ref in hash
			char * dim_name = get_string_in_array(dim_names_sa, dim);
			hv_store(dim_labels_hv, dim_name, strlen(dim_name), newRV_inc((SV *) dim_labels_av), 0);

			// Free names
			destroy_string_array(dim_labels_sa);
		}
		destroy_string_array(dim_names_sa);

		// Return reference
		RETVAL = newRV_inc((SV *) dim_labels_hv);
	OUTPUT:
		RETVAL

void
hdf5_store(file, points_sv)
		void * file
		SV * points_sv
	PREINIT:
		struct hdf5_file_st * file_st = (struct hdf5_file_st *) file;
		AV * points_av;
		hsize_t ** coords;
		double * values;
		hsize_t count, rank, index;
		int coord_index, coord_count;
		I32 length;
	CODE:
		// Defensive coding:
		if (file_st == NULL) {
			puts("Cannot write into null file handle");
			exit(1);
		}

		// Dereference array ref
		points_av = (AV *) SvRV(points_sv);

		// Allocate memory
		count = av_len(points_av) + 1;
		coords = calloc(count, sizeof(hsize_t *));
		values = calloc(count, sizeof(double));

		// Go through data
		rank = get_file_rank(file_st->file);
		for (index = 0; index < count; index++) {
			// Extract datapoint hash ref
			SV ** point_sv = av_fetch(points_av, index, 0);

			// Dereference hash ref
			HV * point_hv = (HV *) SvRV(* point_sv);

			// Allocate C matrix row
			coords[index] = calloc(rank, sizeof(hsize_t));

			// Store coordinates
			coord_count = hv_iterinit(point_hv);
			for (coord_index = 0; coord_index < rank + 1; coord_index++) {
				// Iteration
				HE * hash_entry = hv_iternext(point_hv);

				// Extract info from hash entry
				char * dim_name = hv_iterkey(hash_entry, &length);
				if (strcmp(dim_name, "value") == 0)
					continue;
				SV ** dim_sv = hv_fetch(file_st->dim_indices, dim_name, strlen(dim_name), 0);
				if (dim_sv == NULL) {
					printf("Dimension '%s' unknown!\n", dim_name);
					exit(1);
				}
				int dim = SvIV(*dim_sv);
				hid_t value = SvIV(hv_iterval(point_hv, hash_entry));

				// Updating C data
				coords[index][dim] = value;
			}

			// Store value
			SV ** value = hv_fetch(point_hv, "value", 5, 0);
			values[index] = SvNV(*value);
		}

		// Store into file
		store_values(file_st->file, count, coords, values);

		// Clean up data
		for (index = 0; index < count; index++) {
			free(coords[index]);
		}
		free(coords);
		free(values);

SV * 
hdf5_fetch(file, constraints_hv)
		void * file
		HV * constraints_hv
	PREINIT:
		struct hdf5_file_st * file_st = (struct hdf5_file_st *) file;
		hsize_t rank;
		int index, dim, constraint_count, value;
		bool * set_dims;
		hsize_t * constraints;
		StringResultTable * table;
		AV * results_av;
		char * dim_name;
		int * dim_lengths;
		HV * row_hv;
		I32 length;
	CODE:
		// Allocating dynamic arrays
		rank = get_file_rank(file_st->file);
		set_dims = calloc(rank, sizeof(bool));
		constraints = calloc(rank, sizeof(hsize_t));

		// Reading constraint hash ref to fill two C arrays (set_dims and constraints) with values
		constraint_count = hv_iterinit(constraints_hv);
		for (index = 0; index < constraint_count; index++) {
			// Iteration
			HE * hash_entry = hv_iternext(constraints_hv);

			// Extract info from hash entry
			dim_name = hv_iterkey(hash_entry, &length);
			dim = SvIV(*hv_fetch(file_st->dim_indices, dim_name, strlen(dim_name), 0));
			value = SvIV(hv_iterval(constraints_hv, hash_entry));

			// Updating C data
			set_dims[dim] = 1;
			constraints[dim] = value;
		}

		// Query the file
		table = fetch_string_values(file_st->file, set_dims, constraints);

		// Produce array ref of hash refs for output from C-style StringResultTable object
		results_av = newAV();
		for (index = 0; index < table->rows; index++) {
			// Create hash ref from row of C-style StringResultTable object
			row_hv = newHV();
			for (dim = 0; dim < table->columns; dim++) {
				hv_store(row_hv, table->dims[dim], file_st->dim_name_lengths[dim], newSVpv(table->coords[index][dim], 0), 0);
			}
			hv_store(row_hv, "value", 5, newSVnv(table->values[index]), 0);

			// Append to output array ref
			av_push(results_av, newRV((SV*) row_hv));
		}

		// Cleaning up dynamically allocated arrays
		free(set_dims);
		free(constraints);
		destroy_string_result_table(table);

		RETVAL = newRV((SV *) results_av);
	OUTPUT:
		RETVAL

void
hdf5_close(file)
		void * file
	PREINIT:
		struct hdf5_file_st * file_st = (struct hdf5_file_st *) file;
	CODE:
		close_file(file_st->file);

void
hdf5_set_log(value)
		SV * value;
	CODE:
		set_hdf5_log(SvIV(value));
