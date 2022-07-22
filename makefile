CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
BIN_DIR = ./bin

cutNapAdapter: cutNapAdapter.o cutNapAdapterMain.o bioUtils.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o $(BIN_DIR)/cutNapAdapter cutNapAdapter.o cutNapAdapterMain.o bioUtils.o $(LIBS) 

cutNapAdapter.o: cutNapAdapter.cpp cutNapAdapter.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c cutNapAdapter.cpp

cutNapAdapter_main.o: cutNapAdapterMain.cpp
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c cutNapAdapterMain.cpp

bioUtils.o: bioUtils.cpp bioUtils.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c bioUtils.cpp
	
clean:
	rm -f *.o
