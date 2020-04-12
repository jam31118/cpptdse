include ../var.mk

all: libpropagator.a

OBJS = propagator-on-box-1d.o \
	   propagator-on-box-1d-under-field.o \
	   propagator-for-spharm-zero-m.o

libpropagator.a: $(OBJS)
	$(AR) rcs -o $@ $^

propagator-on-box-1d.o: propagator-on-box-1d.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

propagator-on-box-1d-under-field.o: propagator-on-box-1d-under-field.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

propagator-for-spharm-zero-m.o: propagator-for-spharm-zero-m.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

clean:
	$(RM) *.o *.a

