all: installed

SOURCE_DIR = build/libcvd-20121025
TARBALL = build/libcvd-20121025.tar.gz
TARBALL_URL = http://www.edwardrosten.com/cvd/libcvd-20121025.tar.gz
UNPACK_CMD = tar xzf
TARBALL_PATCH = aligned_mem.patch
include $(shell rospack find mk)/download_unpack_build.mk

INSTALL_DIR = $(shell rospack find cvd)

installed: $(SOURCE_DIR)/unpacked
	mkdir -p $(INSTALL_DIR)
	cd $(SOURCE_DIR) && ./configure --prefix=$(INSTALL_DIR) --without-v4l1buffer CXXFLAGS="-I$(shell rospack find TooN)/include" && make install
	touch installed
	
clean:
	-cd $(SOURCE_DIR) && make clean
	rm -rf $(SOURCE_DIR) $(INSTALL_DIR)/bin $(INSTALL_DIR)/lib $(INSTALL_DIR)/include installed
	
wipe: clean
	rm -rf build
	
.PHONY : clean download