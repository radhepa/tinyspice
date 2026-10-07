# tinyspice - build with GNU make (on Windows: mingw32-make)
CC      := gcc
CFLAGS  := -std=c99 -O2 -Wall -Wextra -pedantic
LDLIBS  := -lm

ifeq ($(OS),Windows_NT)
  EXE := .exe
  RM  := del /Q
  PY  := python
else
  EXE :=
  RM  := rm -f
  PY  := python3
endif

CORE := src/analysis.c src/devices.c src/netlist.c src/solver.c
HDRS := $(wildcard src/*.h)

all: tinyspice$(EXE)

tinyspice$(EXE): $(CORE) src/main.c $(HDRS)
	$(CC) $(CFLAGS) -o $@ $(CORE) src/main.c $(LDLIBS)

tests/test_core$(EXE): tests/test_core.c $(CORE) $(HDRS)
	$(CC) $(CFLAGS) -Isrc -o $@ tests/test_core.c $(CORE) $(LDLIBS)

test: tests/test_core$(EXE)
	./tests/test_core$(EXE)

validate: tinyspice$(EXE)
	$(PY) tests/validate.py

clean:
	-$(RM) tinyspice$(EXE) tests/test_core$(EXE)

.PHONY: all test validate clean
