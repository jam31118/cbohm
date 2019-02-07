include ../var.mk

LIB_DIR=../lib/

OPTIMIZE_FLAGS = -O3

all: install

install: libcbohm.so
	mkdir -p $(LIB_DIR) && mv $^ $(LIB_DIR)

libcbohm.so: velocity.o
	$(CXX) -shared $^ -o $@ -fPIC $(LAPACK_LIB_FLAGS)
#	$(AR) rcs $@ $^

velocity.o: velocity.cc
	$(CXX) -c $< -o $@ $(OPTIMIZE_FLAGS) $(LAPACK_LIB_FLAGS) -fPIC #-DDEBUG

clean:
	$(RM) -f *.a *.o *.so $(LIB_DIR)/*.so $(LIB_DIR)/*.a


