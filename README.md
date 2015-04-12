Ensembl HDF5
============

Copyright holder: [EMBL-EBI](http://www.ebi.ac.uk) (Apache 2 License)

Wrapper for the Ensembl Perl API around large numerical arrays in HDF5 format.

Requirements
------------

This package depends on [PDL](http://pdl.perl.org/) and [PDL::IO::HDF5](http://search.cpan.org/~chm/PDL-IO-HDF5-0.6501/hdf5.pd). (In my experience, CPAN cannot install them properly, so I did it manually.)

Data structures
---------------

All data, whether input or output, is assumed to come as an arrayref containing hashfrefs, where each hashref associates one value to each key. The only reserved keyword is _value_ which points to the numerical value contained in the matrix, e.g:

```
my $original_data = [
  {gene => 'A', snp => 'rs1', value=>.1},
  {gene => 'B', snp => 'rs2', value=>.2}
];
```

Creating a new database
-----------------------

```
my $aa = new Bio::EnsEMBL::HDF5::ArrayAdaptor($filename);
$aa->set_data($original_data);
```

Querying a database
-------------------

All the data except for undefined or null values:
```
my $data_points = $aa->get();
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
my $data_points = $aa->get({gene => 'A'});
```
which returns:
```
[
  {snp => 'rs1', value=>.1}
]
```

Contact us
----------

Please email comments or questions to the public [Ensembl developers' list](http://lists.ensembl.org/mailman/listinfo/dev).

Questions may also be sent to the [Ensembl help desk](http://www.ensembl.org/Help/Contact).
