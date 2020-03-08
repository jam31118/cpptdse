SRC_DIR = ./src

all: install

install:
	$(MAKE) -C $(SRC_DIR)

clean:
	$(MAKE) clean -C $(SRC_DIR)

