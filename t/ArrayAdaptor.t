=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::HDF5::ArrayAdaptor;
use File::Temp qw/tempfile/;
use Data::Dumper;
use feature qw(say);

# For chatty output
#Bio::EnsEMBL::HDF5::set_log(1);

my ($fh, $filename) = tempfile();
print "Creating file\n";
my $aa = new Bio::EnsEMBL::HDF5::ArrayAdaptor(-FILENAME => $filename, -SIZES => {gene => 2, snp => 2}, -LABEL_LENGTHS => {gene => 1, snp => 3}, -DBNAME => ':memory:');
$aa->store_dim_labels('gene', ['A', 'B']);
$aa->store_dim_labels('snp', ['rs1']);
$aa->store_dim_labels('snp', ['rs2']);


my $gene_names = $aa->get_dim_labels("gene");
ok(ref $gene_names eq 'HASH');
ok( scalar keys(%{$gene_names}) == 2);
ok( defined $gene_names->{A} );
ok( defined $gene_names->{B} );

my $original_data = [
  {gene => 'A', snp => 'rs1', value=>.1},
  {gene => 'B', snp => 'rs2', value=>.2}
];

print "Storing data\n";
$aa->store($original_data);

print "Fetching all data\n";
my @output_data = @{$aa->fetch()};

ok(scalar(@output_data) == 2);

foreach my $data_point (@output_data) {
  ok(defined $data_point->{gene});
  ok(defined $data_point->{snp});
  ok(defined $data_point->{value});
  if ($data_point->{gene} eq 'A') {
      ok($data_point->{snp} eq 'rs1');
      ok(abs($data_point->{value} - .1) < 1e-4);
    } elsif ($data_point->{gene} eq 'B') {
      ok($data_point->{snp} eq 'rs2');
      ok(abs($data_point->{value} - .2) < 1e-4);
    }
}

print "Fetching selected data\n";
@output_data = @{$aa->fetch({gene => 'A'})};

ok(scalar(@output_data) == 1);
my $data_point = pop @output_data;
ok(!defined $data_point->{gene});
ok(defined $data_point->{snp});
ok(defined $data_point->{value});
ok($data_point->{snp} eq 'rs1');
ok(abs($data_point->{value} - .1) < 1e-4);

# Test whether an error is raised when an unkown gene is requested
ok(eval {$aa->fetch({gene => 'C'}); 0;} || 1);

$aa->close;

done_testing;

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
