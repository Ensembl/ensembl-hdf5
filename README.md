Ensembl HDF5
============

Copyright holder: [EMBL-EBI](http://www.ebi.ac.uk) (Apache 2 License)

Wrapper for the Ensembl Perl API around large numerical arrays in HDF5 format.

Requirements
------------

This package requires an installation of the [HDF5 library](https://www.hdfgroup.org/HDF5/).

Installation
------------

If you did not install the HDF5 library on a standard location (e.g. usr/shared/lib) then you may need to tell the C compiler about it. 

In c/Makefile change these lines:
```
-INC=
+INC=-I/path/to/include/directory/
-LIB_PATHS=-L./
+LIB_PATHS=-L./ -L/path/to/lib/directory
```

In xs/Makefile.PL change these lines:
```
-    LIBS              => ['-L../c -lhdf5_wrapper -lhdf5'],
+    LIBS              => ['-L../c -lhdf5_wrapper -L/path/to/lib/directory -lhdf5'],
-    INC               => '-I../c',
+    INC               => '-I../c -I/path/to/include/directory',
```

Type 'make'. Unit tests are run automatically.

If things are not working, here is a decomposition of the make process:
```
# C library
cd c
make

# XS Library
cd ../xs
perl Makefile.PL
make
mkdir -p  `echo $PERL5LIB | sed -e 's/:.*//'`/auto/Bio/EnsEMBL/HDF5
cp ./blib/arch/auto/Bio/EnsEMBL/HDF5/HDF5.so `echo $PERL5LIB | sed -e 's/:.*//'`/auto/Bio/EnsEMBL/HDF5/HDF5.so
perl t/*

# Perl Library
cd ..
perl t/*
```

Data structures
---------------

All data, whether input or output, is assumed to come as an arrayref containing hashfrefs, where each hashref associates one value to each key. There can be as many keys as you wish, so long as every hashref has the same keys. The only reserved and required keyword is _value_ which points to the numerical value contained in the matrix, e.g:

```
my $original_data = [
  {gene => 'A', snp => 'rs1', value=>.1},
  {gene => 'B', snp => 'rs2', value=>.2}
];
```

Creating a new database
-----------------------

```
my %dim_labels = (
  gene => ['A', 'B'],
  snp => ['rs1', 'rs2']
);

# Compute the number of labels for each dimension
my $dim_sizes = {
  gene => scalar @{$dim_labels{gene}},
  snp => scalar @{$dim_labels{snp}}
}

# Compute the maximum label length for each dimension
my $dim_label_lengths = {
  gene => max map length @{$dim_labels{gene}},
  snp => max map length @{$dim_labels{snp}}
}

# Create the basic file
my $aa = new Bio::EnsEMBL::HDF5::ArrayAdaptor(
  -FILENAME => $filename, 
  -SIZES => $dim_sizes,
  -LABEL_LENGTHS => $dim_label_lengths
);

# Load the labels
$as->store_dim_labels('gene', $dim_labels{gene});
$as->store_dim_labels('snp', ['rs1']);
$as->store_dim_labels('snp', ['rs2']);

$aa->store($original_data);
```

Note how in this example the labels for the 'snp' dimension are added in multiple steps. You can do this for any dimension. This is conveninent if there are so many labels that a maintaining a list in memory would be problematic

Querying a database
-------------------

All the data except for undefined or null values:
```
my $data_points = $aa->fetch();
```
which returns:
```
[
  {gene => 'A', snp => 'rs1', value=>.1},
  {gene => 'B', snp => 'rs2', value=>.2}
]
```

A subset of the data:
```
my $data_points = $aa->fetch({gene => 'A'});
```
which returns:
```
[
  {snp => 'rs1', value=>.1}
]
```

You can add as many constraints as you wish. Each dimension either has a fixed value or none.

Longer example
--------------

See modules/t/ArrayAdaptor.t for a longer example.

Contact us
----------

Please email comments or questions to the public [Ensembl developers' list](http://lists.ensembl.org/mailman/listinfo/dev).

Questions may also be sent to the [Ensembl help desk](http://www.ensembl.org/Help/Contact).
