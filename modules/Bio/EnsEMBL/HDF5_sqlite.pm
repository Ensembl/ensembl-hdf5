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

HDF5 Mockup: A replacement class for the HDF5 module using SQLite for testing purposes.

=cut

package Bio::EnsEMBL::HDF5_sqlite;

use strict;
use warnings;
use Data::Dumper;

use List::Util qw/ max /;

use Bio::EnsEMBL::DBSQL::DBConnection;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration use Bio::EnsEMBL::HDF5 ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
  hdf5_close
  hdf5_create
  hdf5_fetch
  hdf5_get_dim_labels
  hdf5_open
  hdf5_store
  hdf5_store_dim_labels
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
  
);

our $VERSION = '0.01';


=head2 create

    Constructor
    Argument [1] : Path to HDF5 file
    Argument [2] : Hash ref of { dimension name => length }
    Argument [3] : Hash ref of { dimension name => longest length }

=cut

sub hdf5_create {
  my ($filename, $dim_sizes, $dim_label_lengths) = @_;

  if (! defined $filename) {
    die("Cannot create HDF5 adaptor around undef filename!\n");
  }

  if (!-z $filename) {
    unlink $filename;
  }

  my $sqlite = Bio::EnsEMBL::DBSQL::DBConnection->new(-DBNAME => $filename, -DRIVER => 'SQLite'),

  my @dim_names = keys %$dim_sizes;

  # Create table for dim names
  my $max_dim_name_length = max(map(length, @dim_names));
  my $sql_command = "
  CREATE TABLE IF NOT EXISTS dim_names (
    name VARCHAR($max_dim_name_length)
  )
  ";
  $sqlite->do($sql_command);

  # Store dim names
  $sql_command = "INSERT INTO dim_names (name) VALUES (?)";
  $sqlite->db_handle->begin_work;
  my $sth = $sqlite->prepare($sql_command);
  $sth->bind_param_array(1, \@dim_names);
  $sth->execute_array({});
  $sqlite->db_handle->commit;

  # Create table for dim labels
  foreach my $dim_name (keys %$dim_label_lengths) {
    $sql_command = "
    CREATE TABLE IF NOT EXISTS $dim_name (
      label VARCHAR($dim_label_lengths->{$dim_name})
    )
    ";
    $sqlite->do($sql_command);
  }

  # Create table for main matrix
  my $column_descriptions = join(", ", map { $_." VARCHAR(".$dim_label_lengths->{$_}.")" } @dim_names);
  $sql_command = "
  CREATE TABLE IF NOT EXISTS matrix (
    $column_descriptions,
    value FLOAT
  )
  ";
  $sqlite->do($sql_command);
  $sqlite->db_handle->disconnect;
}

=head2 DESTROY
  Destructor
  Description: Restores default state of the STDOUT filehandle as it is a copy
               and may not flush correctly.
=cut

sub DESTROY {
  my $self = shift;
  if ($self->{'stdout'}) {
    close $self->{'filehandle'};
  }
}

=head2 open

  Opens existing file
  Argument [1] : File location
  Returntype   : Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub hdf5_open {
  my ($filename) = @_;
  return Bio::EnsEMBL::DBSQL::DBConnection->new(-DBNAME => $filename, -DRIVER => 'SQLite'),
}

=head2 store_dim_labels

  Store dimension labels
  Argument [1] : Bio::EnsEMBL::DBSQL::DBConnection
  Argument [2] : Dimension name
  Argument [2] : Listref of labels for that dimension

=cut

sub hdf5_store_dim_labels {
  my ($sqlite, $dim_name, $dim_labels) = @_;
  my $sql_command = "INSERT INTO $dim_name (label) VALUES (?)";
  my $sth = $sqlite->prepare($sql_command);
  $sth->execute_array({}, $dim_labels);
}

=head2 _get_all_dim_names

  Get all dimension names
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection
  Returntype: Listref(dimension names)

=cut

sub _get_all_dim_names {
  my ($sqlite) = @_;
  my $sth = $sqlite->prepare("SELECT name FROM dim_names");
  $sth->execute;
  my @names = map { $_->[0] } @{$sth->fetchall_arrayref};
  return \@names;
}

=head2 get_all_dim_labels

  Get all labels associated to all dimensions
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection
  Returntype: Hashref { dimension name => listref(dimension labels) }

=cut

sub get_all_dim_labels {
  my ($sqlite) = @_;
  my %hash = map { $_ => hdf5_get_dim_labels($sqlite, $_) } @{_get_all_dim_names($sqlite)};
  return \%hash;
}

=head2 get_dim_labels

  Get all labels associated to particular dimension
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection
  Argument [2]: Dimension name
  Returntype: Listref of dimension labels

=cut

sub hdf5_get_dim_labels {
  my ($sqlite, $dim_name) = @_;
  my $sql = "SELECT label FROM $dim_name";
  my $sth = $sqlite->prepare($sql);
  $sth->execute;
  my @labels = map { $_->[0] } @{$sth->fetchall_arrayref};
  return \@labels;
}

=head2 store

  Stores a bunch of datapoints into the matrix
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection
  Argument [2]: Listref of hashrefs { dimension name => dimension label, value => scalar }

=cut

sub hdf5_store {
  my ($sqlite, $points) = @_;
  my $dim_names = _get_all_dim_names($sqlite);
  push @$dim_names, "value";

  my $sql_command = "INSERT INTO matrix (".join(",", @$dim_names).") VALUES (".join(",", map {"?"} @$dim_names).")";
  my $sth = $sqlite->prepare($sql_command);
  for (my $i = 1; $i <= scalar @$dim_names; $i++) {
    my @column = map { $_->{$dim_names->[$i-1]} } @$points;
    $sth->bind_param_array($i, \@column);
  }

  $sth->execute_array({});

}

=head2 fetch

  Fetches all values that fit a given pattern
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection
  Argument [2]: Hashref { dimension name => required dimension_value }
  Returntype: Listref of hashrefs { dimension name => dimension label, value => scalar }

=cut

sub hdf5_fetch {
  my ($sqlite, $constraints) = @_;

  # Remove null constraints
  foreach my $key (keys %$constraints) {
    if (!defined $constraints->{$key}) {
      delete $constraints->{$key};
    }
  } 
  my $dim_names = _get_all_dim_names($sqlite);
  my $dim_labels = get_all_dim_labels($sqlite);
  push @$dim_names, "value";

  my @constrained_dims = keys %$constraints;
  my %hash = map { $_ => 1 } @constrained_dims;
  my @free_dims = grep { ! exists $hash{$_} } @$dim_names;
  my $free_dims_string = join(", ", @free_dims);
  my $constraints_string = '';
  if (scalar @constrained_dims) {
    defined $constraints->{$_} or delete $constraints->{$_} for keys %{$constraints};
    if( scalar(keys %{$constraints}) > 0 ) {
      $constraints_string = "WHERE ".join(" AND ", map { "$_ = $constraints->{$_}" } @constrained_dims);
    }
  }

  my $sql_command = "SELECT $free_dims_string FROM matrix $constraints_string";
  my $sth = $sqlite->prepare($sql_command);
  $sth->execute;
  my @array = ();
  while (my $row = $sth->fetchrow_hashref) {
    my %hash = map { $_ => $dim_labels->{$_}->[$row->{$_}]} keys $row;
    $hash{value} = $row->{value};
    push @array, \%hash;
  }
  return \@array;
}

=head2 close

  Closes connection
  Argument [1]: Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub hdf5_close {
  my ($sqlite) = @_;
  $sqlite->db_handle->disconnect;
}

1;
