all: installed

SOURCE_DIR = build/gvars-3.0
TARBALL = build/gvars-3.0.tar.gz
TARBALL_URL = http://www.edwardrosten.com/cvd/gvars-3.0.tar.gz
UNPACK_CMD = tar xzf
#TARBALL_PATCH = cstdio.patch
include $(shell rospack find mk)/download_unpack_build.mk

INSTALL_DIR = $(shell rospack find gvars3)

installed: $(SOURCE_DIR)/unpacked
	mkdir -p $(INSTALL_DIR)
	cd $(SOURCE_DIR) && ./configure --prefix=$(INSTALL_DIR) CXXFLAGS="-I$(shell rospack find TooN)/include" && make install
	touch installed
	
clean:
	-cd $(SOURCE_DIR) && make clean
	rm -rf $(SOURCE_DIR) $(INSTALL_DIR)/lib $(INSTALL_DIR)/include installed
	
wipe: clean
	rm -rf build
	
.PHONY : clean download