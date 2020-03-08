include ../var.mk

all: libpropagator.a

libpropagator.a: 1d.o
	$(AR) rcs -o $@ $^

1d.o: 1d.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

clean:
	$(RM) *.o *.a

