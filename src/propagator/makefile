include ../var.mk

all: libpropagator.a

libpropagator.a: propagator-on-box-1d.o
	$(AR) rcs -o $@ $^

propagator-on-box-1d.o: propagator-on-box-1d.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

clean:
	$(RM) *.o *.a

