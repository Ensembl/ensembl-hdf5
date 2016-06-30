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
use Data::Dumper;
use feature qw(say);
$| = 1; # Autoflushes all print statements

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
	-user => 'anonymous'
);

print "Building initial file\n";
my $eqtl_adaptor = Bio::EnsEMBL::HDF5::EQTLAdaptor->new(
	-filename         => 'test.hd5',
	-core_db_adaptor  => $registry->get_DBAdaptor('human', 'core'),
	-var_db_adaptor   => $registry->get_DBAdaptor('human', 'variation'),
	-tissues          => ['arm'],
	-statistics       => ['beta','p-value'],
	-dbfile           => 'test.sqlite3',
	-snp_ids          => 'Whole_Blood.snps.txt',
	-gene_ids         => 'Whole_Blood.genes.txt',
);

## Stash the content of the files
my @data_points = ();

# Go through all the lines in the file
open my $fh, "gunzip -c Whole_Blood.cis.eqtl.10.gz |";
<$fh>; # Skip first header line
while (my $line = <$fh>) {
  chomp $line;
  my @items = split("\t", $line);
  $items[1] =~ s/\.[0-9]+//;
  push @data_points, {gene => $items[1], snp => $items[0], tissue => 'arm', statistic => 'beta', value => $items[2]};
  push @data_points, {gene => $items[1], snp => $items[0], tissue => 'arm', statistic => 'p-value', value => $items[4]};
}
close $fh;

$eqtl_adaptor->store(\@data_points);

print Dumper($eqtl_adaptor->fetch({gene => "ENSG00000223972"}));
print Dumper($eqtl_adaptor->fetch({gene => "ENSG00000223972"}, 1));
print Dumper($eqtl_adaptor->fetch({gene => "ENSG00000223972", tissue => 'arm'}, 1));
eval {print Dumper($eqtl_adaptor->fetch({gene => "ENSG00000223972", tissue => 'leg'}, 1))};

$eqtl_adaptor->close;
