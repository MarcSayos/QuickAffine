# Compiler and flags
CC = gcc
CFLAGS = -Wall -O3 -g
LDFLAGS =

# Target executable
TARGET = quickaffine

# Source files
SRC = main.c GapAffine_SWG.c GapAffine_Windowed.c GapAffine_Banded.c
OBJ = $(SRC:.c=.o)

# Header files
HEADERS = main.h GapAffine_SWG.h GapAffine_Windowed.h  GapAffine_Banded.h

# Default rule to build the executable
all: $(TARGET)

# Linking the object files to create the executable
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

# Compiling each source file into object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule to remove object files and the executable
clean:
	rm -f $(OBJ) $(TARGET)

# Phony targets to prevent issues with files named 'all' or 'clean'
.PHONY: all clean
