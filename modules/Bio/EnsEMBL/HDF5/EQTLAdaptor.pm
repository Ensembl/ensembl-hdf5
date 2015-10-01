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

EQTLAdaptor - An array adaptor for HDF5 files specialised in eQTL
data. In particular, it populates gene and variation IDs from the 
Ensembl DBs, and silently converts user inputs into Ensembl IDs

=cut

package Bio::EnsEMBL::HDF5::EQTLAdaptor;

use strict;
use warnings;

use base qw( Bio::EnsEMBL::HDF5::ArrayAdaptor );

use POSIX qw/ strftime /;
use File::Copy qw/ copy /;
use List::Util qw/reduce max/;
use File::Temp qw/ tempfile /;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;

=head2 new

    Constructor
    Argument [1] : Hash
      -FILENAME        : path to HDF5 file (required: to be opened or created)
      -CORE_DB_ADAPTOR : Bio::EnsEMBL::DBSQL::DBAdaptor (required if creating new HDF5)
      -VAR_DB_ADAPTOR  : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor (required if creating new HDF5)
      -TISSUES         : Array ref of string names (required if creating new HDF5)
      -DBFILE          : path to SQLite3 file, if non standard
    Returntype   : Bio::EnsEMBL::HDF5::EQTLAdaptor

=cut

sub new {
  my $class = shift;
  my ($filename, $core_db, $variation_db, $tissues, $statistics, $db_file, $snp_id_file) = rearrange(['FILENAME','CORE_DB_ADAPTOR','VAR_DB_ADAPTOR','TISSUES','STATISTICS','DBFILE','SNP_IDS'], @_);

  if (! defined $filename) {
    die("Cannot create HDF5 adaptor around undef filename!\n");
  }

  $db_file ||= $filename . ".sqlite3";
  my ($fh, $temp) = tempfile;

  if (-e $db_file) {
    copy($db_file, $temp);
  }

  my $self;
  ## If creating a new database
  if (! -e $filename || -z $filename) {
    my $curated_snp_id_file = _curate_variant_names($variation_db, $snp_id_file);
    `sort -k1,1 -k2,2n $curated_snp_id_file`;
  
    my $snp_count = `wc -l $curated_snp_id_file | sed -e 's/ .*//'`;
    chomp $snp_count;
    my $snp_max_length = `awk 'length(\$3) > max {max = length(\$3)} END {print max}' $curated_snp_id_file`;
    my $gene_stats = _fetch_gene_stats($core_db);
    $self = $class->SUPER::new(
	    -FILENAME => $filename, 
	    -SIZES => {
		    gene => $gene_stats->{count}, 
		    snp => $snp_count,
		    tissue => scalar @$tissues,
		    statistic => scalar @$statistics,
	    }, 
	    -LABEL_LENGTHS => {
		    gene => $gene_stats->{max_length},
		    snp => $snp_max_length, 
		    tissue => max(map(length, @$tissues)),
		    statistic => max(map(length, @$statistics)),
	    }, 
	    -DBNAME => $temp,
    );
    $self->_store_gene_labels($core_db);
    $self->_store_variation_labels($variation_db, $curated_snp_id_file);
    $self->store_dim_labels('tissue', $tissues);
    $self->store_dim_labels('statistic', $statistics);
    $self->index_tables;
    copy($temp, $db_file);
  } else {
    $self = $class->SUPER::new(-FILENAME => $filename, -DBNAME => $db_file);
  }

  $self->{tissue_ids} = $self->dim_indices('tissue');
  $self->{gene_ids} = $self->dim_indices('gene');
  $self->{statistic_ids} = $self->dim_indices('statistic');
  $self->{variation_adaptor} = $variation_db->get_adaptor("variation");
  $self->{gene_adaptor} = $core_db->get_adaptor("gene");
  bless $self, $class;
  return $self;
}

=head2 _fetch_gene_stats

  Counts all the genes and their max ID length throughout the database
  Argument [1] : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub _fetch_gene_stats {
  my ($core_db) = @_;
  my $count = 0;
  my $max_length = -1;
  my $gene_adaptor = $core_db->get_adaptor("gene");
  my $slice_adaptor = $core_db->get_adaptor("slice");

  print "Going through gene IDs\n";
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    foreach my $gene (@{$gene_adaptor->fetch_all_by_Slice($slice)}) {
      $count++;
      if ($max_length < length $gene->stable_id) {
        $max_length = length $gene->stable_id;
      }
    }
  }

  $count > 0 || die("No genes found?");

  return {count => $count, max_length => $max_length};
}

=head2 _store_gene_labels

  Counts all the genes and their max ID length throughout the database
  Argument [1] : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub _store_gene_labels {
  my ($self, $core_db) = @_;
  my $gene_adaptor = $core_db->get_adaptor("gene");
  my $slice_adaptor = $core_db->get_adaptor("slice");

  print "Storing gene IDs\n";
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel')}) {
    my @labels = map { $_->stable_id } sort _by_position @{$gene_adaptor->fetch_all_by_Slice($slice)};
    $self->store_dim_labels('gene', \@labels);
  }
}

=head2 _sort_and_curate_variant_names

=cut

sub _curate_variant_names {
  my ($variation_db, $snp_id_file) = @_;

  print "Storing dataset variation IDs in temporary table\n";
  my $va = $variation_db->get_adaptor("variation");

  ## We start by storing all the variation labels
  ## unsorted into a temporary file
  my ($fh, $temp) = tempfile;
  open my $file, "<", $snp_id_file;
  while (my $line = <$file>) {
    chomp $line;
    my $feature = $va->fetch_by_name($line)->get_all_VariationFeatures->[0];
    if (! defined $feature) {
      next;
    }
    print $fh join("\t", ($feature->seq_region_name, $feature->seq_region_start, $feature->variation_name))."\n";
  }
  close $file;

  return $temp;
}

=head2 _store_variation_labels

  Counts all the variations and their max ID length throughout the database
  Argument [1] : Bio::EnsEMBL::DBSQL::DBAdaptor
  Argument [2] : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor

=cut

sub _store_variation_labels {
  my ($self, $variation_db, $snp_id_file) = @_;

  ## We then read the sorted file 
  print "Streaming variation IDs from database\n";
  open my $file2, "<", $snp_id_file;
  my @labels = ();
  while (my $line = <$file2>) {
    chomp $line;
    my @items = split("\t", $line);
    my $name = $items[2];
    push @labels, $name;

    # If buffer full, push into SQLite and HDF5 storage
    if (scalar @labels > 10000) {
       $self->store_dim_labels('snp', \@labels);
       @labels = ();
    }
  }

  # Flush out remaining buffer
  if (scalar @labels) {
    $self->store_dim_labels('snp', \@labels);
  }
}

=head2 _by_position

  Compares two Bio::EnsEMBL::Feature objects ($a and $b) based on coordinates
  Returntype : int 

=cut 

sub _by_position {
  my $a_pos = $a->strand > 0? $a->seq_region_start: $a->seq_region_end;
  my $b_pos = $b->strand > 0? $b->seq_region_start: $b->seq_region_end;
  my $chrom_cmp = $a->seq_region_name cmp $b->seq_region_name;
  return $chrom_cmp ? $chrom_cmp : $a_pos <=> $b_pos;
}

=head2 _convert_to_Ensembl

  Converts external IDs (gene names or variation rsIDs) into Ensembl IDs, possibly uncovering redundancies
  Argument [1] : Hashrefs, where each key/value pair is a dimension/value
    Possible dimensions are:
    * gene   : Ensembl gene ID or common name
    * snp    : rsID or Ensembl ID
    * tissue : string
    * value  : floating point scalar
  Returntype : Arrayref of hashrefs

=cut 

sub _convert_coords {
  my ($self, $coords) = @_;
  my ($gene, $snp, $tissue, $statistic);

  if (defined $coords->{'snp'}) {
    $snp = $self->_get_numerical_value('snp', $coords->{snp});
    if (! defined $snp) {
      my $id = $self->{variation_adaptor}->fetch_by_name($coords->{'snp'})->name;
      $snp = $self->_get_numerical_value('snp', $id);
    }

    if (! defined $snp) {
      die("Did not recognize ".$coords->{'snp'}."\n");
    }
  }

  if (defined $coords->{'gene'}) {
    ## Hacky regexp of Ensembl IDs
    if ($coords->{'gene'} =~ /ENS[A-Z]+[0-9]{11}/) {
      $gene = $coords->{'gene'};
    } else {
      $gene = $self->{gene_adaptor}->fetch_all_by_external_name($coords->{'gene'})->display_id->[0];
    }
    my $gene_id = $self->{gene_ids}{$gene}; 

    if (!defined $gene_id) {
      my $EnsemblGene = $self->{gene_adaptor}->fetch_by_stable_id($gene);
      if (!defined $EnsemblGene) {
        die("No gene with name $gene\n");
      }
      $gene_id = $self->{gene_ids}{$EnsemblGene->stable_id}; 
    }

    if (!defined $gene_id) {
      die("Did not recognize ".$coords->{'gene'}."\n".(scalar keys $self->{gene_ids}));
    } else {
      $gene = $gene_id;
    }
  }

  if (defined $coords->{tissue}) {
    $tissue = $self->{tissue_ids}{$coords->{tissue}};
    if (! defined $tissue) {
      die("Did not recognise tissue $coords->{tissue}\n");
    }
  }

  if (defined $coords->{statistic}) {
    $statistic = $self->{statistic_ids}->{$coords->{statistic}};
    if (! defined $statistic) {
      die("Did not recognise statistic $coords->{statistic}\n");
    }
  }

  my $res = {};
  if (defined $gene) {
    $res->{gene} = $gene;
  } 
  if (defined $snp) {
    $res->{snp} = $snp;
  }
  if (defined $tissue) {
    $res->{tissue} = $tissue;
  } 
  if (defined $statistic) {
    $res->{statistic} = $statistic;
  }
  if (exists $coords->{value}) {
    $res->{value} = $coords->{value};
  }
  return $res;
}

1;
