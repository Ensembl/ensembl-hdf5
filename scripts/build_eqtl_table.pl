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

Bio::EnsEMBL::HDF5::set_log(1); # Small talk from the C layer

main();

sub main {
  my %options = %{get_options()};
  my $fill = 0;

  if (-e $options{hdf5}) {
    $fill = 1;
  }

  $options{sqlite3} ||= $options{hdf5} . ".sqlite3";

  print "Opening file $options{sqlite3}\n";
  my $eqtl_adaptor = build_eqtl_table(\%options);

  ## Stash the content of the files
  if ($fill) {
    foreach my $file (@{$options{files}}) {
      my $tissue = shift @{$options{tissues}};
      print "Loading file $file for tissue $tissue\n";
      extract_file_data($eqtl_adaptor, $file, $tissue);
    }
  }

  $eqtl_adaptor->close;
}

sub get_options {
  my %options = ();
  GetOptions(\%options, "help=s", "host|h=s", "port|p=s", "species|s=s", "user|u=s", "pass|p=s", "tissues|t=s@", "files|f=s@","hdf5=s", "sqlite3|d=s");
  if (defined $options{tissues} 
      && defined $options{files} 
      && (scalar @{$options{tissues}} != scalar @{$options{files}})) {
    die("You should provide as many tissue names as filenames!");
  }
  if (!defined $options{tissues} || scalar @{$options{tissues}} < 1) {
	  die("No tissues!");
  }
  return \%options;
}

sub extract_snp_ids_from_file {
  my ($file) = @_;
  defined $file || die;
  my ($fh, $temp) = tempfile();
  run("gzip -dc $file | cut -f1 | grep -v ^# | grep ^rs | sort | uniq > $temp");
  return $temp;
}

sub extract_snp_ids {
  my ($files, $filename) = @_;
  my $out = $filename.".gtex.snp.ids";
  if (-e $out && !-z $out) {
    return $out;
  }
  scalar @$files || die;
  my @temps = map({extract_snp_ids_from_file($_)} @$files);
  my ($fh, $temp) = tempfile();
  run('sort -m ' .join(" ", @temps) ." | uniq > $out");
  return $out;
}

sub build_eqtl_table {
  my ($options) = @_;

  print "Extracting SNP ids from input files\n";
  my @files = @{$options->{files}};
  my $snp_id_file = extract_snp_ids(\@files, $options->{hdf5});
  
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_registry_from_db(
	  -host => $options->{host}, 
	  -user => $options->{user}, 
	  -pass => $options->{pass}, 
	  -port => $options->{port}, 
	  -db_version => 73,
  );

  print "Building initial file\n";
  return Bio::EnsEMBL::HDF5::EQTLAdaptor->new(
            -filename => $options->{hdf5},
            -core_db_adaptor => $registry->get_DBAdaptor('human', 'core'),
            -var_db_adaptor => $registry->get_DBAdaptor('human', 'variation'),
            -tissues => $options->{tissues},
            -statistics => ['beta','p-value'],
            -dbfile => $options->{sqlite3},
	    -snp_ids => $snp_id_file,
  )
}

sub extract_file_data {
  my ($eqtl_adaptor, $filename, $tissue) = @_;
  my @result = ();

  open my $fh, "gunzip -c $filename |";
  <$fh>; # Skip first header line

  # Go through all the lines in the file
  while (my $line = <$fh>) {
    chomp $line;
    my @items = split("\t", $line);
    if (!($items[0] =~ /^rs[0-9]*/)) {
      next;
    }
    $items[1] =~ s/\.[0-9]+//;
    push @result, {gene => $items[1], snp => $items[0], tissue => $tissue, statistic => 'beta', value => $items[2]};
    push @result, {gene => $items[1], snp => $items[0], tissue => $tissue, statistic => 'p-value', value => $items[4]};
    if (scalar @result > 1000000) {
      my $start = time;
      print "STORING\n";
      $eqtl_adaptor->store(\@result);
      print "STORED\t". (time - $start) ."\n";
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

=head2 run

  Description: Wrapper function for system calls
  Arg1: Command line
  Returntype: undef
  Side effects: Runs command, prints out error in case of failure

=cut

sub run {
  my ($cmd) = @_;
  print "Running $cmd\n";
  my $exit_code = system($cmd);
  if ($exit_code != 0) {
    die("Failure when running command\n$cmd\n")
  }
}
