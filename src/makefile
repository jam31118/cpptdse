LIB_DIR=../lib

OPT_FLAGS = -O3

all: install

install: libtdse.so
	if [ ! -d "$$LIB_DIR" ]; then mkdir -p $(LIB_DIR); mv $^ $(LIB_DIR); fi

libtdse.so: tridiag.o
	$(CXX) -shared $^ -o $@ -fPIC

tridiag.o: tridiag.cc
	$(CXX) -c -o $@ $< -fPIC $(OPT_FLAGS)

clean:
	$(RM) *.o *.so $(LIB_DIR)/*.so
	if [ -d "$$LIB_DIR" ]; then rmdir $(LIB_DIR); fi

