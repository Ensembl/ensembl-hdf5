Ensembl HDF5
============

Copyright holder: [EMBL-EBI](http://www.ebi.ac.uk) (Apache 2 License)

Wrapper for the Ensembl Perl API around large numerical arrays in HDF5 format.

Requirements
------------

This package requires an installation of the [HDF5 library](https://www.hdfgroup.org/HDF5/).

Installation
------------

Type 'make'. Unit tests are run automatically.

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
my $dim_labels = {
  gene => ['A', 'B'],
  snp => ['rs1', 'rs2']
};
my $aa = new Bio::EnsEMBL::HDF5::ArrayAdaptor(-FILENAME => $filename, -LABELS => $dim_labels);
$aa->store($original_data);
```

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
