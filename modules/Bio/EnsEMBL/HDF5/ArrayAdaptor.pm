=pod

=head1 LICENSE

  Copyright (c) 1999-2015 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

  http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

ArrayAdaptor - A generic array adaptor for HDF5 files 

=cut

package Bio::EnsEMBL::HDF5::ArrayAdaptor;

use strict;
use warnings;

use PDL;
use PDL::Char;
use PDL::IO::HDF5;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

=head2 new

    Constructor
    Argument     : HDF5 Filename 
    Returntype   : Bio::EnsEMBL::HDF5::ArrayAdaptor

=cut

sub new {
  my ($class, $filename) = @_;
  
  # The HDF5 library gets confused if it finds an empty stub 
  if (-z $filename) {
    unlink $filename;
  }
  my $exists = -e $filename;

  my $self = {
    hdf5 => new PDL::IO::HDF5($filename),
    dim_names => undef,
    dim_indices => undef,
    dim_labels => undef,
    label_indices => undef
  };

  bless $self, $class;

  # Might as well preload labels...
  if ($exists) {
    $self->_load_dim_data;
  }
  
  return $self;
}

=head2 _load_dim_data

    Loads dimension names and labels into memory

=cut

sub _load_dim_data {
  my $self = shift;
  $self->_load_dim_names;
  $self->_load_dim_labels;
}

=head2 _load_dim_names

    Loads dimension names into a PDL 

=cut

sub _load_dim_names {
  my $self = shift;
  $self->{dim_names} = $self->{hdf5}->group('/dims')->dataset("names")->get;
  $self->{dim_indices} = {};
  foreach my $dim_index (0..($self->{dim_names}->getdim(1) - 1)) {
    $self->{dim_indices}{$self->{dim_names}->atstr($dim_index)} = $dim_index;
  }
}

=head2 _load_dim_labels

    Loads string labels into memory

=cut

sub _load_dim_labels {
  my $self = shift;
  $self->{dim_labels} = {};
  $self->{label_indices} = {};
  foreach my $dim_index (0..($self->{dim_names}->getdim(1) - 1)) {
    my $dim_name = $self->{dim_names}->atstr($dim_index);
    my @res = @{_load_dim_labels2($self->{hdf5}->group('/labels')->dataset($dim_name)->get)};
    $self->{dim_labels}{$dim_name} = $res[0];
    $self->{label_indices}{$dim_name} = $res[1];
  } 
}

sub _load_dim_labels2 {
  my ($dim_labels) = @_;
  my %label_indices = ();
  my $index = 0;

  foreach my $index (0..($dim_labels->getdim(1) - 1)) {
    $label_indices{$dim_labels->atstr($index)} = $index;
  }

  return [$dim_labels, \%label_indices];
}

=head2 set_data

  Arguments [1]: Arrayref of hashrefs: {dim_name1 => label1, ...,  "value" => scalar}

=cut

sub set_data {
  my ($self, $data_points) = @_;
  $self->_set_dim_labels(_extract_dim_labels($data_points));
  $self->_load_dim_data;
  $self->{hdf5}->dataset('data')->set($self->_create_pdl($data_points));
}

=head2 _extract_dim_label

  Arguments [1]: Arrayref of hashrefs: {dim_name1 => label1, ...,  "value" => scalar}
  Returntype: Hashref of label => arrayref of values

=cut

sub _extract_dim_labels {
  my ($data_points) = @_;
  if (length @$data_points == 0) {
    throw("Empty dataset");
  }

  my %data_point0 = %{$data_points->[0]};
  my $dim = scalar keys %data_point0;
  my %dim_labels = ();

  foreach my $data_point (@$data_points) {
    if (scalar keys %$data_point != $dim) {
      throw("Error, data point does not have the right keys. Expected:\n".join(", ", keys %data_point0)."\nFound:\n".join(", ", keys %$data_point));
    }
    foreach my $key (keys %$data_point) {
      if ($key eq 'value') {
        next;
      }
      if (! exists $data_point0{$key}) {
        throw("Error, data point does not have the right keys. Expected:\n".join(", ", keys %data_point0)."\nFound:\n".join(", ", keys %$data_point));
      }
      $dim_labels{$key} ||= {};
      $dim_labels{$key}{$data_point->{$key}} = 1;
    }
  }

  my %final = ();
  foreach my $key (keys %dim_labels) {
    my @labels = keys %{$dim_labels{$key}};
    $final{$key} = \@labels;
  }

  return \%final;
}

=head2 _set_dim_labels

  Arguments [1]: Hashref of dimension label -> array ref of values

=cut

sub _set_dim_labels {
  my ($self, $dim_labels) = @_;
  foreach my $dim_name (keys %$dim_labels) {
    $self->{hdf5}->group('/labels')->dataset($dim_name)->set(PDL::Char->new($dim_labels->{$dim_name}));
  }
  $self->{hdf5}->group('/dims')->dataset('names')->set(PDL::Char->new(sort {scalar(@{$dim_labels->{$a}}) <=> scalar(@{$dim_labels->{$b}})} keys %$dim_labels));
}

=head2 _create_pdl

  Arguments [1]: Arrayref of hashrefs: {dim_name1 => label1, ...,  "value" => scalar} 
  Returntype: PDL

=cut

sub _create_pdl {
  my ($self, $data_points) = @_;

  my $pdl = $self->_empty_pdl;
  foreach my $data_point (@$data_points) {
    $self->_populate_pdl($pdl, $data_point);
  }
  return $pdl;
}

=head2 _empty_pdl

  Returntype: Null PDL with the dimensions expected by the number of dimensions and their labels

=cut


sub _empty_pdl {
  my ($self) = @_;
  my @dim_sizes = ();
  foreach my $dim_index (0..($self->{dim_names}->getdim(1) - 1)) {
    my $dim_name = $self->{dim_names}->atstr($dim_index);
    $dim_sizes[$dim_index] = $self->{dim_labels}{$dim_name}->getdim(1);
  }
  return zeroes(@dim_sizes);
}

=head2 _populate_pdl

  Arguments [1]: Hashref of dimension name => label

=cut

sub _populate_pdl {
  my ($self, $pdl, $data_point) = @_;
  my @coords = ();
  foreach my $dim_name (keys %$data_point) {
    if ($dim_name eq "value") {
      next;
    }
    $coords[$self->{dim_indices}{$dim_name}] = $self->{label_indices}{$dim_name}{$data_point->{$dim_name}};
  }
  set($pdl, @coords, $data_point->{value});
}

=head2 get

  Arguments [1]: Hashref of dimension name => label

=cut

sub get {
  my ($self, $constraints) = @_;
  $constraints ||= {};
  my ($start, $finish) = @{$self->_constraint_corners($constraints)};
  my $pdl = $self->{hdf5}->dataset('data')->get(pdl($start), pdl($finish));
  return $self->_extract_sparse_pdl_data($constraints, $pdl);
}

=head2 _constraint_corners

  Defines the corners of a multidimensional slice of an HDF5 dataset to be extracted
  Arguments [1]: Key-value hashref defining the constraints of the search
  Returntype: Arrayref with two PDL vectors, representing the start and end corner

=cut

sub _constraint_corners {
  my ($self, $constraints) = @_;
  my $dataset = $self->{hdf5}->dataset('data');
  my @start = ();
  my @finish = ();
  foreach my $dim_index (0..($self->{dim_names}->getdim(1) - 1)) {
    my $dim_name = $self->{dim_names}->atstr($dim_index);
    if (exists $constraints->{$dim_name}) {
       my $selected_index = $self->{label_indices}{$dim_name}{$constraints->{$dim_name}};
       if (!defined $selected_index) {
         throw("Unknown $dim_name in request: $constraints->{$dim_name}");
       }
       $start[$dim_index] = $selected_index;
       $finish[$dim_index] = $selected_index;
    } else {
       $start[$dim_index] = 0;
       $finish[$dim_index] = $self->{dim_labels}{$dim_name}->getdim(1) - 1;
    }
  }
  return [\@start, \@finish];
}

=head2 _extract_sparse_pdl_data

  Recursively breaks down a PDL extracted from the HDF5 file into an arrayref of hashrefs
  Arguments [1]: Key-value hashref defining the constraints of the search
  Arguments [2]: PDL
  Returntype: Arrayref of hashrefs

=cut

sub _extract_sparse_pdl_data {
  my ($self, $constraint, $pdl) = @_;
  my @unconstrained_dims = grep {!exists $constraint->{$_}} keys $self->{dim_indices};
  my @data_points = ();
  my $coords = whichND($pdl);
  foreach my $row (dog $coords) {
    my $data_point = {value => $pdl->at(@{unpdl $row})};
    foreach my $dim (@unconstrained_dims) {
      $data_point->{$dim} = $self->{dim_labels}{$dim}->atstr($row->at($self->{dim_indices}{$dim}));
    }
    push @data_points, $data_point;
  }
  return \@data_points;
}

=head2 _extract_dense_pdl_data

  Recursively breaks down a PDL extracted from the HDF5 file into an arrayref of hashrefs
  Arguments [1]: Key-value hashref defining the constraints of the search
  Arguments [2]: PDL
  Returntype: Arrayref of hashrefs

=cut

sub _extract_dense_pdl_data {
  my ($self, $constraint, $pdl) = @_;
  my $dimensionality = scalar($pdl->getndims);
  if ($dimensionality == 0) {
    return [{value => $pdl->at}];
  } else {
    my $dim_name = $self->{dim_names}->atstr($dimensionality - 1);
    my @data_points = ();
    if (exists $constraint->{$dim_name}) {
      foreach my $sub_pdl (dog($pdl)) {
        push @data_points, @{$self->_extract_pdl_data($constraint, $sub_pdl)};
      }
    } else {
      my $index = 0;
      foreach my $sub_pdl (dog($pdl)) {
        push @data_points, @{_add_dim($dim_name, $self->{dim_labels}{$dim_name}->atstr($index++), $self->_extract_pdl_data($constraint, $sub_pdl))};
      }
    }
    return \@data_points;
  }
}

=head2 _add_dim

  Adds key/value pairs to all the hashrefs in an arrayref
  Arguments [1]: Key scalar
  Arguments [2]: Value scalar
  Arguments [3]: Arrayref of hashrefs
  Returntype: Arrayref of hashrefs

=cut

sub _add_dim {
  my ($dim_name, $label, $data_points) = @_;
  foreach my $data_point (@$data_points) {
    $data_point->{$dim_name} = $label;
  }
  return $data_points;
}

1;
