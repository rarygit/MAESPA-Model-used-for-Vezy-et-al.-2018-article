##
## Makefile for building a Maespa, could be nicer...
## Martin De Kauwe, 05/12/2013
## Sligthly modified by R Vezy, 07/11/2019

SRCDIR := src
BUILDDIR := build

# OS specific part
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	PROG =	maespa.out
	RM = rm -rf
	RMDIR = rm -rf $(BUILDDIR)/*
endif

ifeq ($(OS),Windows_NT)
	PROG =	maespa.exe
	RM = del /F /Q
	RMDIR = del /F /Q $(BUILDDIR)\*.*
endif

SRCS =	switches.f90 default_conditions.f90 maestcom.f90  metcom.f90 getmet.f90 \
				maindeclarations.f90 inout.f90 maespa.f90 physiol.f90 \
        radn.f90 unstor.f90 utils.f90 watbal.f90

SRCFILES := $(addprefix $(SRCDIR)/,$(SRCS))

OBJS =	switches.o default_conditions.o maestcom.o metcom.o getmet.o maindeclarations.o inout.o \
				maespa.o physiol.o radn.o unstor.o utils.o watbal.o
#OBJSFILES= $(SRCFILES:.f90=.o)
OBJSFILES := $(addprefix $(BUILDDIR)/,$(OBJS))
MODFILES := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.d,$(SRCFILES))

LIBS =
INCLS =

F90 = gfortran

FFLAGS = -g -fbounds-check -finit-local-zero -Wuninitialized -ftrapv -ffree-form -ffree-line-length-none -O3 -mcmodel=large -ffast-math

build: $(BUILDDIR)/$(PROG)

all: build cleanmod
	@echo "Build complete, MAESPA executable is in the build directory"

$(BUILDDIR)/$(PROG): $(OBJSFILES)
	$(F90) $(FFLAGS) -o $@ $(OBJSFILES) $(LIBS) $(INCLS)
	@echo 'Finished building target: $@'
	@echo ' '

#$(SRCDIR)/%.o: $(SRCDIR)/%.f90
#		$(F90) $(FFLAGS) -c $<
$(BUILDDIR)/%.o : $(SRCDIR)/%.f90
		$(F90) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(BUILDDIR)/%.mod

clean:
	$(RM) *.mod
	$(RMDIR)
	@echo 'Clean complete'

cleanmod:
	$(RM) *.mod
# If something doesn't work right, have a 'make debug' to
# show what each variable contains.
debug:
	@echo "SRCFILES= $(SRCFILES)"
	@echo "PROG = $(PROG)"
	@echo "OBJSFILES = $(OBJSFILES)"
	@echo "F90 = $(F90)"
	@echo "FFLAGS = $(FFLAGS)"
	@echo "command: $(F90) $(FFLAGS) -o $@ $(OBJFILES) $(LIBS) $(INCLS)"
	@echo "RMDIR: $(RMDIR)"

.SUFFIXES: $(SUFFIXES) .f90

default_conditions.o: switches.o
getmet.o: maestcom.o metcom.o switches.o
maindeclarations.o: maestcom.o
inout.o: maestcom.o switches.o
maespa.o: maestcom.o metcom.o switches.o maindeclarations.o
physiol.o: maestcom.o metcom.o
radn2.o: maestcom.o
unstor.o: maestcom.o
utils.o: maestcom.o
watbal.o: maestcom.o metcom.o
