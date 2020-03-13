include var.mk

LIB_DIR=../lib

all: install

install: libtdse.so
	if [ ! -d "$$LIB_DIR" ]; then mkdir -p $(LIB_DIR); mv $^ $(LIB_DIR); fi

libtdse.so: tridiag.o array.o ./propagator/libpropagator.a libwf.a
	$(CXX) -shared *.o -o $@ \
		-Wl,--whole-archive ./propagator/libpropagator.a ./wf/libwf.a \
		-Wl,--no-whole-archive $(MATRIX_LIB_FLAGS)

tridiag.o: tridiag.cc
	$(CXX) -c -o $@ $< -fPIC $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

array.o: array.cc
	$(CXX) -c -o $@ $< -fPIC $(OPT_FLAGS)

./propagator/libpropagator.a:
	$(MAKE) -C ./propagator

libwf.a:
	$(MAKE) -C ./wf

clean:
	$(RM) *.o *.so $(LIB_DIR)/*.so
	if [ -d "$$LIB_DIR" ]; then rmdir $(LIB_DIR); fi
	$(MAKE) clean -C ./propagator
	$(MAKE) clean -C ./wf

