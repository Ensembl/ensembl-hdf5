=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

use strict;
use warnings;

use Test::More;
use File::Temp qw/tempfile/;

BEGIN { use_ok('Bio::EnsEMBL::HDF5_sqlite') };

use Bio::EnsEMBL::HDF5_sqlite;

# For chatty output
#Bio::EnsEMBL::HDF5_sqlite::set_log(1);

my ($fh, $filename) = tempfile();
Bio::EnsEMBL::HDF5_sqlite::create($filename, {gene => 2, snp => 2}, {gene => 1, snp => 3});

ok(my $hdfh = Bio::EnsEMBL::HDF5_sqlite::open($filename));

Bio::EnsEMBL::HDF5_sqlite::store_dim_labels($hdfh, 'gene', ['A','B']);
Bio::EnsEMBL::HDF5_sqlite::store_dim_labels($hdfh, 'snp', ['rs1']);
Bio::EnsEMBL::HDF5_sqlite::store_dim_labels($hdfh, 'snp', ['rs2']);

my $labels = Bio::EnsEMBL::HDF5_sqlite::get_all_dim_labels($hdfh);

my @gene_names = sort { $a cmp $b } @{Bio::EnsEMBL::HDF5_sqlite::get_dim_labels($hdfh, "gene")};

ok($gene_names[0] eq "A");
ok($gene_names[1] eq "B");

my $original_data = [
  {gene => 0, snp => 0, value=>.1},
];

my $original_data2 = [
  {gene => 1, snp => 1, value=>.2}
];

# Inserting data into file
Bio::EnsEMBL::HDF5_sqlite::store($hdfh, $original_data);
Bio::EnsEMBL::HDF5_sqlite::store($hdfh, $original_data2);

# Pulling out all the data
my @output_data = @{Bio::EnsEMBL::HDF5_sqlite::fetch($hdfh, {})};

ok(scalar(@output_data) == 2);

foreach my $data_point (@output_data) {
  ok(defined $data_point->{gene});
  ok(defined $data_point->{snp});
  ok(defined $data_point->{value});
  if ($data_point->{gene} eq 'A') {
    ok($data_point->{snp} eq 'rs1');
    ok($data_point->{value} == .1);
  } elsif ($data_point->{gene} eq 'B') {
    ok($data_point->{snp} eq 'rs2');
    ok($data_point->{value} == .2);
  }
}

# Pulling out a subset of the data
@output_data = @{Bio::EnsEMBL::HDF5_sqlite::fetch($hdfh, {gene => 0})};

ok(scalar(@output_data) == 1);
my $data_point = pop @output_data;
ok(!defined $data_point->{gene});
ok(defined $data_point->{snp});
ok(defined $data_point->{value});
ok($data_point->{snp} eq 'rs1');
ok($data_point->{value} == .1);

# Test whether an error is raised when an unkown gene is requested
#@output_data = @{Bio::EnsEMBL::HDF5_sqlite::fetch($hdfh, {gene => 2})};

done_testing;

unlink $filename;

# Little convenience function for debugging
sub print_dataset {
  my $data_set = shift;
  print "{\n";
  foreach my $data_point (@$data_set) {
    foreach my $key (keys %$data_point) {
      print "\t$key:$data_point->{$key}"; 
    }  
    print "\n";
  }
  print "}\n";
}
