use 5.014002;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Bio::EnsEMBL::HDF5',
    VERSION_FROM      => '../modules/Bio/EnsEMBL/HDF5.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => '../modules/Bio/EnsEMBL/HDF5.pm', # retrieve abstract from module
       AUTHOR         => 'Daniel Zerbino <zerbino@ebi.ac.uk>') : ()),
    LIBS              => ['-L../c -lhdf5_wrapper -L/nfs/users/nfs_d/dz1/.local/lib -lhdf5'], # e.g., '-lm'
    DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
    INC               => '-I../c -I/nfs/users/nfs_d/dz1/.local/include', # e.g., '-I. -I/usr/include/other'
	# Un-comment this if you add C files to link with later:
    # OBJECT            => '$(O_FILES)', # link all the C files too
);
