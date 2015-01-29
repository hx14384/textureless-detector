all: installed

SOURCE_DIR = build/TooN-2.0
TARBALL = build/TooN-2.0.tar.gz
TARBALL_URL = http://www.edwardrosten.com/cvd/TooN-2.0.tar.gz
UNPACK_CMD = tar xzf
include $(shell rospack find mk)/download_unpack_build.mk

INSTALL_DIR = $(shell rospack find TooN)

export PKG_CONFIG = $(INSTALL_DIR)

installed: $(SOURCE_DIR)/unpacked
	mkdir -p $(INSTALL_DIR)
	mkdir $(INSTALL_DIR)/lib
	mkdir $(INSTALL_DIR)/lib/pkgconfig
	cd $(SOURCE_DIR) && ./configure --prefix=$(INSTALL_DIR) && make install
	touch installed
	
clean:
	-cd $(SOURCE_DIR) && make clean
	rm -rf $(SOURCE_DIR) $(INSTALL_DIR)/lib $(INSTALL_DIR)/include installed
	
wipe: clean
	rm -rf build
	
.PHONY : clean download