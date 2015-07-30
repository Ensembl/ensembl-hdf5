CFLAGS=-g -Wall -fPIC
INC=
LIB_PATHS=-L./
LIBS=-lhdf5_wrapper -lhdf5

default: lib test

lib: hdf5_wrapper.o
	ar rcs libhdf5_wrapper.a hdf5_wrapper.o

test: hdf5_test.o lib
	${CC} ${CFLAGS} ${LIB_PATHS} hdf5_test.o ${LIBS} -o test
	./test
	rm TEST.hd5

%.o: %.c; ${CC} ${CFLAGS} ${INC} ${OPTS} -c $< -o $@

clean:
	rm -Rf *.o *.a

