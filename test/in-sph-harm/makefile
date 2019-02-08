include ../../var.mk

LIB_DIR=../../lib
LIB_FLAGS = -L$(LIB_DIR) -lcbohm

all: main

main: main.cc $(LIB_DIR)/libcbohm.so
	$(CXX) -o $@ $< $(LIB_FLAGS)

$(LIB_DIR)/libcbohm.so:
	+$(MAKE) -C ../../src

clean:
	+$(MAKE) clean -C ../../src
	$(RM) -f main

