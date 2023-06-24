# Makefile for the Sphere Packing project

# Author Info
EMAIL=arstgri@gmail.com

# GLPK Directory
GLPK_DIR = /home/amirreza/apps/glpk

# Compiler and Linker Options
CC = gcc
CFLAGS = -I$(SRC_DIR) -I$(GLPK_DIR)/include 
LDFLAGS = -L$(GLPK_DIR)/lib -lglpk -lm

# Source Files
SRC_DIR = src
SRCS = $(wildcard $(SRC_DIR)/*.c)

# Build Directory
BUILD_DIR = build

# Object Files
OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))

# Output Executable
TARGET = packing

# Default Target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
