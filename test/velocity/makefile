LIB_DIR=../../lib
LIB_FLAGS = -L$(LIB_DIR) -lcbohm

all: main

main: main.cc libcbohm.so
	$(CXX) -o $@ $< $(LIB_FLAGS)

libcbohm.so:
	+$(MAKE) -C ../../src

clean:
	+$(MAKE) clean -C ../../src
	$(RM) -f main

