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

package Bio::EnsEMBL::HDF5::ArrayAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::HDF5;

=head2 new

    Constructor
    Argument [1] : HDF5 Filename 
    Argument [2] : Optional: Hash ref of dimension name => array ref of allowed values 
    Returntype   : Bio::EnsEMBL::HDF5::ArrayAdaptor

=cut

sub new {
  my $class = shift;
  my ($filename, $dim_sizes, $dim_label_lengths, $dbname) = rearrange(['FILENAME','SIZES', 'LABEL_LENGTHS','DBNAME'], @_);

  defined $filename || die ("Must specify HDF5 filename!");

  $dbname ||= $filename . ".sqlite3";
  
  my $self = {
    hdf5 => undef,
    sqlite3 => Bio::EnsEMBL::DBSQL::DBConnection->new(-DBNAME => $dbname, -DRIVER => 'SQLite'),
    dim_labels => undef,
    st_handles => {}
  };

  bless $self, $class;

  if (! (-e $filename) || -z $filename) {
    if (! defined $dim_sizes || !defined $dim_label_lengths) {
      die("Cannot create a new HDF5 store without dimensional label info");
    }

    if (-e $filename && -z $filename) {
      unlink $filename;
    }

    Bio::EnsEMBL::HDF5::create($filename, $dim_sizes, $dim_label_lengths);
    my @dim_names = keys %$dim_sizes;
    $self->_create_sqlite3_file($filename, \@dim_names);
  }
  
  $self->{hdf5} = Bio::EnsEMBL::HDF5::open($filename); 

  return $self;
}

=head2 _create_sqllite3_file

    Argument [1] : HDF5 Filename 
    Argument [2] : Hash ref of dimension name => array ref of allowed values 

=cut

sub _create_sqlite3_file {
  my ($self, $filename, $dim_names) = @_;

  foreach my $key (@$dim_names) {
    $self->_create_sqlite3_table($key);
  }
}

=head2 _create_sqlite3_table

  Argument [1]: Dimension name

=cut

sub _create_sqlite3_table {
  my ($self, $dim_name) = @_;

  my $sql = "
  CREATE TABLE $dim_name (
    hdf5_index	INTEGER PRIMARY KEY AUTOINCREMENT,
    external_id	INTEGER
  )
  ";
  $self->{sqlite3}->do($sql);
}


=head2 store_dim_labels

  Argument [1]: dim name
  Argument [2]: array of labels

=cut

sub store_dim_labels {
  my ($self, $dim_name, $dim_labels) = @_;
  Bio::EnsEMBL::HDF5::store_dim_labels($self->{hdf5}, $dim_name, $dim_labels);
  $self->_insert_into_sqlite3_table($dim_name, $dim_labels);
}

=head2 _insert_into_sqlite3_table

  Argument [1]: Dimension name
  Argument [2]: Array of ids

=cut

sub _insert_into_sqlite3_table {
  my ($self, $dim_name, $dim_labels) = @_;
  my $sql2 = "INSERT INTO $dim_name (external_id) VALUES (?)";
  my $sth = $self->{sqlite3}->prepare($sql2);
  $sth->execute_array({}, $dim_labels);
}

=head2 _get_numerical_value

  Arguments [1]: Name of dimension / key
  Arguments [2]: Name of value
  Return type: integer index of value within HDF5 matrix

=cut

sub _get_numerical_value {
  my ($self, $dim, $label) = @_;
  if ($dim eq 'value') {
    return $label;
  }
  if (! exists $self->{st_handles}{$dim}) {
    $self->{st_handles}{$dim} = $self->{sqlite3}->prepare("SELECT hdf5_index FROM $dim WHERE external_id=?")
  }
  $self->{st_handles}{$dim}->execute($label);
  my @row =  $self->{st_handles}{$dim}->fetchrow_array;

  scalar @row > 0 || die("Label $label unknown in dimension $dim");
  ## Note that SQLite3 autoincrement starts at 1, instead of 0 for HDF5 indices
  return $row[0] - 1;
}

=head2 _convert_coords

  Arguments [1]: Hash ref of coord string => string label
  Returntype: Hash ref of coord string => integer index

=cut

sub _convert_coords {
  my ($self, $coords) = @_;
  my $numerical_coords = {};
  foreach my $key (keys %$coords) {
    $numerical_coords->{$key} = $self->_get_numerical_value($key, $coords->{$key});
  }
  return $numerical_coords;
}

=head2 store 

  Arguments [1]: Arrayref of hashrefs: {dim_name1 => label1, ...,  "value" => scalar}

=cut

sub store {
  my ($self, $data_points) = @_;
  my @converted_points = map {$self->_convert_coords($_)} @$data_points;
  Bio::EnsEMBL::HDF5::store($self->{hdf5}, \@converted_points);
}

=head2 fetch 

  Arguments [1]: Hashref of dimension name => label
  Returntype: Arrayref of hashrefs: dimension name => label

=cut

sub fetch {
  my ($self, $constraints) = @_;
  $constraints ||= {};
  my $res = Bio::EnsEMBL::HDF5::fetch($self->{hdf5}, $self->_convert_coords($constraints));
  return $res;
}

=head2 close

=cut

sub close {
  my ($self) = @_;
  Bio::EnsEMBL::HDF5::close($self->{hdf5});
  foreach my $key (keys %{$self->{st_handles}}) {
    $self->{st_handles}{$key}->finish;
  }
  $self->{sqlite3}->db_handle->disconnect;
}

1;
