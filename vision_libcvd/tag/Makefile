all: installed

SOURCE_DIR = build/tag-20090724
TARBALL = build/tag-20090724.tar.gz
TARBALL_URL = http://www.edwardrosten.com/cvd/tag-20090724.tar.gz
UNPACK_CMD = tar xzf
include $(shell rospack find mk)/download_unpack_build.mk

INSTALL_DIR = $(shell rospack find tag)

installed: $(SOURCE_DIR)/unpacked
	mkdir -p $(INSTALL_DIR)
	cd $(SOURCE_DIR) && ./configure --prefix=$(INSTALL_DIR) --with-TooN=$(shell rospack find TooN)/include && make install
	touch installed
	
clean:
	-cd $(SOURCE_DIR) && make clean
	rm -rf $(SOURCE_DIR) $(INSTALL_DIR)/lib $(INSTALL_DIR)/include installed
	
wipe: clean
	rm -rf build
	
.PHONY : clean download