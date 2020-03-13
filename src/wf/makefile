include ../var.mk

all: libwf.a

libwf.a: wavefunction-on-box-1d.o
	$(AR) rcs -o $@ $^

wavefunction-on-box-1d.o: wavefunction-on-box-1d.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS)

clean:
	$(RM) *.o *.a

