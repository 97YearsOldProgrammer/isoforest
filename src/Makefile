# ==============================================================================
# RFHMM Makefile
# ==============================================================================

CC = gcc
TARGET = hmm

# Directories
INCDIR = include
PARSEDIR = parser
ENCDIR = encoder
DECDIR = decoder

# Source files
CORE_SRC = main.c output.c cJSON.c batch.c
PARSER_SRC = $(wildcard $(PARSEDIR)/*.c)
ENCODER_SRC = $(wildcard $(ENCDIR)/*.c)
DECODER_SRC = $(wildcard $(DECDIR)/*.c)

# Object files
OBJS = $(CORE_SRC:.c=.o) \
       $(PARSER_SRC:.c=.o) \
       $(ENCODER_SRC:.c=.o) \
       $(DECODER_SRC:.c=.o)

# Header dependencies
HEADERS = $(wildcard $(INCDIR)/*.h) model.h cJSON.h $(DECDIR)/randomf.h batch.h

# ==============================================================================
# Build Configurations
# ==============================================================================

# Debug (default) - with address sanitizer
CFLAGS_DEBUG = -fsanitize=address -g -I$(INCDIR) -I.

# Release - optimized, no debug
CFLAGS_RELEASE = -O3 -I$(INCDIR) -I.

# OpenMP - parallel processing
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIBOMP := $(shell test -d /opt/homebrew/opt/libomp/include && echo "yes")
    ifeq ($(LIBOMP),yes)
        CFLAGS_OPENMP = -O3 -Xpreprocessor -fopenmp -I$(INCDIR) -I. -I/opt/homebrew/opt/libomp/include
        LDFLAGS_OPENMP = -L/opt/homebrew/opt/libomp/lib -lomp
    endif
else
    CFLAGS_OPENMP = -O3 -fopenmp -I$(INCDIR) -I.
    LDFLAGS_OPENMP = -fopenmp
    LIBOMP = yes
endif

# Default flags
CFLAGS = $(CFLAGS_DEBUG)
LDFLAGS =

# ==============================================================================
# Build Targets
# ==============================================================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm $(LDFLAGS)

release: clean
	@echo "Building release..."
	$(MAKE) CFLAGS="$(CFLAGS_RELEASE)" $(TARGET)

openmp: clean
ifeq ($(LIBOMP),yes)
	@echo "Building with OpenMP ($(shell nproc 2>/dev/null || sysctl -n hw.ncpu) cores available)..."
	$(MAKE) CFLAGS="$(CFLAGS_OPENMP)" LDFLAGS="$(LDFLAGS_OPENMP)" $(TARGET)
else
	@echo "Error: OpenMP not available. Install with: brew install libomp"
	@exit 1
endif

# ==============================================================================
# Pattern Rules
# ==============================================================================

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

$(PARSEDIR)/%.o: $(PARSEDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

$(ENCDIR)/%.o: $(ENCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

$(DECDIR)/%.o: $(DECDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# ==============================================================================
# Test Targets
# ==============================================================================

test: $(TARGET)
	./$(TARGET) -s seq -M ../models/worm.splicemodel -v

rf: $(TARGET)
	./$(TARGET) -s seq -M ../models/worm.splicemodel -S -N 100

clean:
	rm -f $(TARGET) *.o $(PARSEDIR)/*.o $(ENCDIR)/*.o $(DECDIR)/*.o

help: $(TARGET)
	./$(TARGET) --help

.PHONY: all release openmp test rf clean help
