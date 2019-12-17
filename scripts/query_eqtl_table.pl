=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

ArrayAdaptor - A generic array adaptor for HDF5 files

=cut

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw/ zip /;
use File::Copy qw/ copy /;
use File::Temp qw/ tempfile /;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::HDF5::EQTLAdaptor;
$| = 1; # Autoflushes all print statements

Bio::EnsEMBL::HDF5::set_log(2); # Small talk from the C layer

main();

sub main {
  my $options = get_options();

  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(
	  -host => $options->{host},
	  -user => $options->{user},
	  -pass => $options->{pass},
	  -port => $options->{port},
	  -db_version => 73,
  );

  my $eqtl_adaptor = Bio::EnsEMBL::HDF5::EQTLAdaptor->new(
            -filename => $options->{hdf5},
            -core_db_adaptor => $registry->get_DBAdaptor('human', 'core'),
            -var_db_adaptor => $registry->get_DBAdaptor('human', 'variation'),
  );

  my $results = $eqtl_adaptor->fetch( {
      gene        => $options->{gene},
      tissue      => $options->{tissue},
      snp         => $options->{snp},
      statistic   => $options->{statistic},
      chromosome  => $options->{chromosome},
      position    => $options->{position},
      });
  print scalar(@$results)."\n";

  print join("\t", qw/tissue snp gene statistic value chromosome position/)."\n";
  foreach my $result (@$results) {
    foreach my $column (qw/tissue snp gene statistic value chromosome position/) {
      if (defined $options->{$column}) {
        print "*$options->{$column}\t";
      } else {
        print "$result->{$column}\t";
      }
    }
    print "\n";
  }

  $eqtl_adaptor->close;
}

sub get_options {
  my %options = ();
  GetOptions(\%options, "help=s", "host|h=s", "port|p=s", "user|u=s", "pass|p=s", "tissue=s", "gene=s", "snp=s", "statistic=s", "hdf5=s", "sqlite3|d=s");
  if (defined $options{tissues}
      && defined $options{files}
      && (scalar @{$options{tissues}} != scalar @{$options{files}})) {
    die("You should provide as many tissue names as filenames!");
  }
  return \%options;
}
