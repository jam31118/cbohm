SRC_DIR = ./src
TEST_DIR = ./test
TEST_SUBDIRS := $(wildcard $(TEST_DIR)/*/.)

.PHONY: $(TEST_SUBDIRS)
all: install

install:
	+$(MAKE) -C $(SRC_DIR)

test: $(TEST_SUBDIRS)

$(TEST_SUBDIRS):
	+$(MAKE) -C $@

clean:
	+$(MAKE) clean -C $(SRC_DIR)

clean-test:
	for dir in $(TEST_SUBDIRS); do \
		$(MAKE) clean -C $$dir; \
	done

