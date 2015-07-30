default: xs

xs: xs/blib/arch/auto/Bio/EnsEMBL/HDF5

xs/blib/arch/auto/Bio/EnsEMBL/HDF5: c/libhdf5_wrapper.a
	cd xs; perl Makefile.PL; make; make test; make install

c/libhdf5_wrapper.a:
	cd c; make
