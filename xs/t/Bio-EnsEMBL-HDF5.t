use strict;
use warnings;

use Test::More;
use File::Temp qw/tempfile/;

BEGIN { use_ok('Bio::EnsEMBL::HDF5') };

use Bio::EnsEMBL::HDF5;

# For chatty output
#Bio::EnsEMBL::HDF5::set_log(1);

my ($fh, $filename) = tempfile();
Bio::EnsEMBL::HDF5::create($filename, {'gene' => ['A','B'], 'snp' => ['rs1','rs2']});

ok(my $hdfh = Bio::EnsEMBL::HDF5::open($filename));

my $labels = Bio::EnsEMBL::HDF5::get_dim_labels($hdfh);

my $original_data = [
  {gene => 0, snp => 0, value=>.1},
];

my $original_data2 = [
  {gene => 1, snp => 1, value=>.2}
];

# Inserting data into file
Bio::EnsEMBL::HDF5::store($hdfh, $original_data);
Bio::EnsEMBL::HDF5::store($hdfh, $original_data2);

# Pulling out all the data
my @output_data = @{Bio::EnsEMBL::HDF5::fetch($hdfh, {})};

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
@output_data = @{Bio::EnsEMBL::HDF5::fetch($hdfh, {gene => 0})};

ok(scalar(@output_data) == 1);
my $data_point = pop @output_data;
ok(!defined $data_point->{gene});
ok(defined $data_point->{snp});
ok(defined $data_point->{value});
ok($data_point->{snp} eq 'rs1');
ok($data_point->{value} == .1);

# Test whether an error is raised when an unkown gene is requested
#@output_data = @{Bio::EnsEMBL::HDF5::fetch($hdfh, {gene => 2})};

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
