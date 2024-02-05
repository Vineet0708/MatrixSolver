# Makefile for a program with two source files (main.c and functions.c)

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -g -O3
LDFLAGS = -lm -lrt -pg

# Source files and objects
SOURCES = main.c functions.c
OBJECTS = $(SOURCES:.c=.o)

# Executable name
EXECUTABLE = myprogram

# Default target
all: $(EXECUTABLE)

# Linking step
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE) $(LDFLAGS)

# Compilation step for each source file
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target to remove generated files
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
