=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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
use feature qw(say);
$| = 1; # Autoflushes all print statements

Bio::EnsEMBL::HDF5::hdf5_set_log(1); # Small talk from the C layer

main();

sub main {
  my %options = %{get_options()};
  my $fill = 0;

  if (-e $options{hdf5}) {
    $fill = 1;
  }

  my $eqtl_adaptor = build_eqtl_table(\%options);

  ## Stash the content of the files
  if ($fill) {
    foreach my $file (@{$options{files}}) {
      my $tissue = shift @{$options{tissues}};
      say "Loading file $file for tissue $tissue";
      extract_file_data($eqtl_adaptor, $file, $tissue);
    }
  }

  $eqtl_adaptor->close;
}

sub get_options {
  my %options = ();
  GetOptions(\%options, "help=s", "host|h=s", "port|p=s", "species|s=s", "user|u=s", "pass|p=s", "tissues|t=s@", "files|f=s@","hdf5=s", "sqlite3|d=s", "snp_info=s", "gene_info=s");
  if (defined $options{tissues} 
      && defined $options{files} 
      && (scalar @{$options{tissues}} != scalar @{$options{files}})) {
    die("You must provide as many tissue names as filenames!");
  }
  if (!defined $options{tissues} || scalar @{$options{tissues}} < 1) {
	  die("No tissues!");
  }
  $options{species} ||= 'homo_sapiens';
  return \%options;
}

sub build_eqtl_table {
  my ($options) = @_;

  say "Extracting SNP ids from input files";
  my $snp_id_file = $options->{snp_info};
  
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(
	  -host => $options->{host}, 
	  -user => $options->{user}, 
	  -pass => $options->{pass}, 
	  -port => $options->{port}, 
  );

  print "Building initial file\n";
  return Bio::EnsEMBL::HDF5::EQTLAdaptor->new(
          -filename         => $options->{hdf5},
          -core_db_adaptor  => $registry->get_DBAdaptor($options->{species}, 'core'),
          -var_db_adaptor   => $registry->get_DBAdaptor($options->{species}, 'variation'),
          -tissues          => $options->{tissues},
          -statistics       => ['beta','p-value'],
          -dbfile           => $options->{sqlite3},
          -snp_ids          => $options->{snp_info},
          -gene_ids          => $options->{gene_info},
  )
}

sub extract_file_data {
  my ($eqtl_adaptor, $filename, $tissue) = @_;
  my @result = ();
  open my $fh, "gunzip -c $filename |";

  # Go through all the lines in the file
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split("\t", $line);
    if ($items[0] !~ /^rs[0-9]*/) {
      next;
    }
    $items[1] =~ s/\.[0-9]+//;
    push @result, {gene => $items[1], snp => $items[0], tissue => $tissue, statistic => 'beta', value => $items[2]};
    push @result, {gene => $items[1], snp => $items[0], tissue => $tissue, statistic => 'p-value', value => $items[4]};
    if (scalar @result > 1000000) {
      my $start = time;
      say "STORING";
      $eqtl_adaptor->store(\@result);
      say "STORED\t". (time - $start);
      @result = ();
    }
  }

  if (scalar @result > 0) {
    $eqtl_adaptor->store(\@result);
    @result = ();
  }

  close $fh;
  return \@result;
}
