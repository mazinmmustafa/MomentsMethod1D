# Executable
EXE = program.exe

# DLL
DLL = myLib.dll

# Compilers
CC = gcc

# Compiler Library
CLIB = -lgfortran -lquadmath

# Compiler Flags
CFLG = -Wall -Wextra -Ofast

# Directories
BDIR = bin
HDIR = include
ODIR = obj
SDIR = src
LDIR = lib

# Variables
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))

# Targets
all: $(BDIR)/$(EXE)

$(BDIR)/$(EXE): $(COBJ)
	$(CC) -o $@ $^ -I $(HDIR) -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLG) -c -o $@ $^

dll: $(COBJ)
	$(CC) -shared -o $(LDIR)/$(DLL) $^ -I $(HDIR) -L./ $(CLIB)

# Clean
.PHONY: clean

clean:
	del $(BDIR)\$(EXE) $(ODIR)\*.o $(LDIR)\$(DLL)
	cls

cleanMesh:
	del Mesh\*.txt Mesh\*.dat
	cls
