use strict;
use warnings;

use Test::More;
use Bio::EnsEMBL::HDF5::ArrayAdaptor;
use File::Temp qw/tempfile/;

my ($fh, $filename) = tempfile();
my $aa = new Bio::EnsEMBL::HDF5::ArrayAdaptor($filename);

my $original_data = [
  {gene => 'A', snp => 'rs1', value=>.1},
  {gene => 'B', snp => 'rs2', value=>.2}
];

# Inserting data into file
ok($aa->set_data($original_data));

# Pulling out all the data
my @output_data = @{$aa->get()};

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
@output_data = @{$aa->get({gene => 'A'})};

ok(scalar(@output_data) == 1);
my $data_point = pop @output_data;
ok(!defined $data_point->{gene});
ok(defined $data_point->{snp});
ok(defined $data_point->{value});
ok($data_point->{snp} eq 'rs1');
ok($data_point->{value} == .1);

# Test whether an error is raised when an unkown gene is requested
ok(eval {$aa->get({gene => 'C'}); 0;} || 1);

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
