include ../var.mk

all: libwf.a

libwf.a: wavefunction-on-box-1d.o wavefunction-for-spharm-zero-m.o
	$(AR) rcs -o $@ $^

wavefunction-on-box-1d.o: wavefunction-on-box-1d.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS)

wavefunction-for-spharm-zero-m.o: wavefunction-for-spharm-zero-m.cc
	$(CXX) -fPIC -c -o $@ $< $(OPT_FLAGS)

clean:
	$(RM) *.o *.a

