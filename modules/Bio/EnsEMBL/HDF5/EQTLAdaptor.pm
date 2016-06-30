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
use feature qw/say/;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::HDF5_sqlite qw (
      hdf5_store_dim_labels
    );

# use Cache::FileCache;
# DBI->trace(1);
use Data::Dumper;
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
  my ($hdf5_file, $core_db, $variation_db,
    $tissues, $statistics, $db_file, $snp_id_file, $gene_ids) =
  rearrange(['FILENAME','CORE_DB_ADAPTOR','VAR_DB_ADAPTOR',
    'TISSUES','STATISTICS','DBFILE','SNP_IDS', 'GENE_IDS'], @_);

  if (! defined $hdf5_file) {
    die("Cannot create HDF5 adaptor around undef filename!\n");
  }
  # If not hdf5 sqlite3 file has provided, create one using base of hdf5 file
  $db_file ||= $hdf5_file . ".sqlite3";
  my ($fh, $temp) = tempfile;

  if (-e $db_file) {
    copy($db_file, $temp);
  }

  my $self;
  ## If creating a new database
  if (! -e $hdf5_file || -z $hdf5_file) {
    say 'Creating a new DB (no hdf5 file passed or size 0)';
    my $curated_snp_id_file = _curate_variant_names($variation_db, $snp_id_file, $hdf5_file.".snp.ids");

    my $snp_count = `wc -l $curated_snp_id_file | sed -e 's/ .*//'`;
    chomp $snp_count;

    my $snp_max_length = `awk 'length(\$0) > max {max = length(\$0)} END {print max}' $curated_snp_id_file`;
    chomp($snp_max_length);

    my $gene_stats;
    if (defined $gene_ids) {
      my $gene_count = `wc -l $gene_ids | sed -e 's/ .*//'`;
      chomp $gene_count;
      my $gene_max_length = `awk 'length(\$1) > max {max = length(\$1)} END {print max}' $gene_ids`;
      chomp($gene_max_length);
      $gene_stats->{count} = $gene_count;
      $gene_stats->{max_length} = $gene_max_length;
    } else {
      $gene_stats = _fetch_gene_stats($core_db);
    }

    $self = $class->SUPER::new(
      -FILENAME   => $hdf5_file,
      -DBNAME     => $temp,
      -SIZES      => {
        gene      => $gene_stats->{count},
        snp       => $snp_count,
        tissue    => scalar @$tissues,
        statistic => scalar @$statistics,
	      },
      -LABEL_LENGTHS => {
        gene      => $gene_stats->{max_length},
        snp       => $snp_max_length,
        tissue    => max(map(length, @$tissues)),
        statistic => max(map(length, @$statistics)),
	    },
    );
    if (defined $gene_ids) {
      $self->_store_my_gene_labels($gene_ids);
    } else {
      $self->_store_gene_labels($core_db);
    }
    $self->_store_variation_labels($curated_snp_id_file);
    $self->store_dim_labels('tissue', $tissues);
    $self->store_dim_labels('statistic', $statistics);
    $self->index_tables;
    copy($temp, $db_file);
    $self->{variation_adaptor}  = $variation_db->get_adaptor("variation");
    $self->{gene_adaptor}       = $core_db->get_adaptor("gene");
  } else {
    say "$hdf5_file";
    $self = $class->SUPER::new(-FILENAME => $hdf5_file, -DBNAME => $temp);
    $self->{hdf5_file} = $hdf5_file;
  }

  $self->{tissue_ids}         = $self->dim_indices('tissue');
  $self->{gene_ids}           = $self->dim_indices('gene');
  $self->{statistic_ids}      = $self->dim_indices('statistic');

  bless $self, $class;
  return $self;
}

=head2 _fetch_gene_stats

  Counts all the genes and their max ID length throughout the database
  Argument [1] : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub _fetch_gene_stats {
  my ($core_db) = @_;

  my $count         = 0;
  my $max_length    = -1;
  my $gene_adaptor  = $core_db->get_adaptor("gene");
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

=head2 _store_my_gene_labels

  Counts all the genes and their max ID length throughout the database
  Argument [1] : File with Gene Ids.

=cut

sub _store_my_gene_labels {
  my ($self, $filename) = @_;

  open my $file, "<", $filename;
  my @labels= ();
  while (my $line = <$file>) {
    chomp $line;
    push @labels, $line;
  }
  close $file;
  $self->store_dim_labels('gene', \@labels);
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

=head2 _curate_variant_names
  $snps_id_file = list of unique IDs (rs or GTEX) parsed from GTEX file
  $file_hdf5_snps = contains chr, pos, rsID, ID from $snps_id_file.
                    sorted first by chr, then by pos

  Output:
  $seq_region_name	$seq_region_start	$seq_region_end	$rs_id	$old_rs_id	$display_consequence
=cut

sub _curate_variant_names {
  my ($variation_db, $snps_id_file, $file_hdf5_snps) = @_;

  if (-e $file_hdf5_snps && ! -z $file_hdf5_snps) {
    return $file_hdf5_snps;
  }

  print "Storing dataset variation IDs in temporary table\n";
  my $va  = $variation_db->get_adaptor("variation");

  ## We start by storing all the variation labels
  ## unsorted into a temporary file
  my ($out, $temp) = tempfile;
#  die "Empty..." if(-e $file_gtex_snps and -z $file_gtex_snps);
  open my $in, "<", $file_hdf5_snps;
  while (my $line = <$in>) {
    chomp $line;
    $line =~ /^(rs\d+)/;
    my $rsid = $1;
    if (! defined $rsid) {
      die "Expect rsID (/^rs\\d+/), not $rsid";
    }

    my $variant = $va->fetch_by_name($rsid);

    my $features = $variant->get_all_VariationFeatures;
    if (! scalar @$features) {
      next;
    }
    my $feature = $features->[0];

    print $out join("\t", ($feature->seq_region_name, $feature->seq_region_start, $feature->seq_region_end, $feature->variation_name, $rsid, $feature->display_consequence))."\n";
  }

  `sort -k1,1 -k2,2n $temp > $file_hdf5_snps`;

  close $in;
  close $out;
  return $file_hdf5_snps;
}

=head2 _store_variation_labels

  Counts all the variations and their max ID length throughout the database
  Argument [1] : File location

=cut

sub _store_variation_labels {
  my ($self, $snp_id_file) = @_;

  $self->{snp_ids} = {};

  ## We then read the sorted file
  print "Streaming variation IDs from database\n";
  open my $file2, "<", $snp_id_file;
  my @labels = ();
  while (my $line = <$file2>) {
    chomp $line;
    my ($chrom, $start, $end, $name, $given_name, $consequence) = split("\t", $line);
    if (!($given_name eq $name)) {
      $self->{snp_ids}{$given_name} = $name;
    }
    push @labels, join("\t", ($name, $chrom, $start, $end, $consequence));

    # If buffer full, push into SQLite and HDF5 storage
    if (scalar @labels > 10000) {
       hdf5_store_dim_labels($self->{hdf5}, 'snp', \@labels);
       my @rsIds = map {my @array = split("\t", $_); $array[0] } @labels;
       $self->_insert_into_sqlite3_table('snp', \@rsIds);
       @labels = ();
    }
  }

  # Flush out remaining buffer
  if (scalar @labels) {
    hdf5_store_dim_labels($self->{hdf5}, 'snp', \@labels);
    my @rsIds = map {my @array = split("\t", $_); $array[0] } @labels;
    $self->_insert_into_sqlite3_table('snp', \@rsIds);
  }
}

=head2 _load_snp_aliases

  Reads off list of SNP id replacements
  Argument [1] : Bio::EnsEMBL::DBSQL::DBAdaptor
  Argument [2] : Bio::EnsEMBL::Variation::DBSQL::DBAdaptor

=cut

sub _load_snp_aliases {
  my ($self) = @_;

  $self->{snp_ids} = {};

  ## We then read the sorted file
  open my $file, "<", $self->{hdf5_file}.".snp.ids";
  while (my $line = <$file>) {
    chomp $line;
    my @items = split("\t", $line);
    my $name = $items[2];
    my $given_name = $items[3];
    if (!($given_name eq $name)) {
      $self->{snp_ids}{$given_name} = $name;
    }
  }
  close $file;
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

=head2 _convert_coords

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

  if (defined $coords->{snp}) {
    if (! defined $self->{snp_ids}) {
	    $self->_load_snp_aliases;
    }
    if (defined $self->{snp_ids} && exists $self->{snp_ids}{$coords->{snp}}) {
       $snp = $self->_get_numerical_value('snp', $self->{snp_ids}{$coords->{snp}});
    } else {
       $snp = $self->_get_numerical_value('snp', $coords->{snp});
    }
    if (! defined $snp) {
      printf("CONNECTING TO VARIANT SERVER $coords->{snp}\n");
      my $EnsemblSnp = $self->{variation_adaptor}->fetch_by_name($coords->{snp});
      if (!defined $EnsemblSnp) {
        die("No SNP or variant with name $coords->{snp}\n");
      }
      $snp = $self->_get_numerical_value('snp', $EnsemblSnp->name);
    }

    if (! defined $snp) {
      die("Did not recognize ".$coords->{snp}."\n");
    }
  }

  if (defined $coords->{gene}) {
    ## Hacky regexp of Ensembl IDs
    if ($coords->{gene} =~ /ENS[A-Z]+[0-9]{11}/) {
      $gene = $coords->{gene};
    } else {
      printf("CONNECTING TO CORE SERVER $coords->{gene}\n");
      $gene = $self->{gene_adaptor}->fetch_all_by_external_name($coords->{gene})->[0]->stable_id;
    }
    my $gene_id = $self->{gene_ids}{$gene};

    if (!defined $gene_id) {
      printf("CONNECTING TO CORE SERVER2 $gene\n");
      my $EnsemblGene = $self->{gene_adaptor}->fetch_by_stable_id($gene);
      if (!defined $EnsemblGene) {
        die("No gene with name $gene\n");
      }
      $gene_id = $self->{gene_ids}{$EnsemblGene->stable_id};
    }

    if (!defined $gene_id) {
      die("Did not recognize ".$coords->{gene}."\n".(scalar keys $self->{gene_ids}));
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

=head2 fetch_all_tissues

  Returns all known tissue identifiers
  Returntype : List ref of strings

=cut

sub fetch_all_tissues {
  my ($self) = @_;
  return $self->get_dim_labels("tissue");
}

=head2 fetch

  Returns all data subject to constraints
  Arg[1]: hash ref of { $dim => $value } constraints
  Arg[2]: Bool, true if return value augmented with data for web display (snp data and snp consequence)
  Returntype : List ref of hashrefs of {$dim => $value} data points

=cut

sub fetch{
  my ($self, $constraints, $web) = @_;
  my $res = $self->SUPER::fetch($constraints);
  if (! exists $constraints->{snp}) {
    foreach my $correlation (@$res) {
      my ($rs_id, $seq_region_name, $seq_region_start, $seq_region_end, $display_consequence) = split("\t", $correlation->{snp});
      $correlation->{snp} = $rs_id;
      if ($web) {
        $correlation->{seq_region_name} = $seq_region_name;
        $correlation->{seq_region_start} = $seq_region_start;
        $correlation->{seq_region_end} = $seq_region_end;
        $correlation->{display_consequence} = $display_consequence;
      }
    }
  }
  return $res;
}

1;
