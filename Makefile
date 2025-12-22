# Adjust these paths if your folder structure changes
LIB_DIR = ./build/lib
INC_DIR = ./build/include
JULIA_LIB_DIR = ./build/lib/julia

# Compiler flags
CFLAGS = -I$(INC_DIR)
LDFLAGS = -L$(LIB_DIR) -lGradusXSPEC \
          -Wl,-rpath,$(LIB_DIR) \
          -Wl,-rpath,$(JULIA_LIB_DIR)

all: main

main: main.c
	gcc main.c -o main $(CFLAGS) $(LDFLAGS)

clean:
	rm -f main
