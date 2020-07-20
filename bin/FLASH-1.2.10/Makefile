#
# Makefile for FLASH.  Process this using GNU make by simply running
#
# $ make
#
# You may set CPPFLAGS, LDFLAGS, or CC in the environment if you would like to
# use a non-default compiler, extra preprocessor flags (e.g. to find headers),
# or extra linker flags (e.g. to find libraries).

CFLAGS += -O2 -Wall -std=c99 -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64
LDLIBS := -lz -lpthread
OBJ    := combine_reads.o flash.o iostream.o read_io.o read_queue.o read_util.o util.o
EXE    := flash

$(EXE):$(OBJ)

clean:
	rm -f $(OBJ) $(EXE)

.PHONY: clean
