#/*
#This file is part of 3Dec, an accurate base-calling software for sequences.
#
#Copyright (c) 2015, Bo Wang, Academy of Mathematics and Systems Science,
#Chinese Academy of Sciences, Beijing 100190, China

#This Source Code Form is subject to the terms of the Creative Commons 
#Attribution-NonCommercial-ShareAlike 4.0 International Public License.
#If a copy of the licence was not distributed with this file, you can obtain
#one at http://creativecommons.org/licenses/by-nc-sa/4.0/
#*/

SRCDIR = ./src
BINDIR = ./bin
MANDIR = ./man
DESTDIR = /usr/local/bin
BUILDDIR = ./build
CXX = g++
FC = gfortran
CPPFLAGS += -I./include
CFLAGS = -Wall -O3 -funroll-loops -DNDEBUG -std=c++98
ifneq ($(openmp),disabled)
CFLAGS += -fopenmp
endif

LDFLAGS =  -lm
ifneq ($(dependency),included)
LDFLAGS += -llinear
endif

INCFLAGS = 
DEFINES =
objects =  basecall.o cacc_funs.o correct_acc.o io.o operation.o
header = cacc_funs.h classtype.h purity.h
target = oacc cacc bsca 3Dec 3Dec-train
VPATH = $(SRCDIR):$(LIBDIR):$(BUILDDIR):$(BINDIR)

.PHONY: default
default: 3Dec 3Dec-train

ifeq ($(dependency),included)
obj_linear = tron.o linear.o blas.a
CPPFLAGS += -I./include/liblinear  
DEFINES += -D_linear_2

$(obj_linear) : ./include/liblinear/train
	cp ./include/liblinear/$@ $(BUILDDIR)/

./include/liblinear/train :
	$(MAKE) -C ./include/liblinear
	cp ./include/liblinear/blas/blas.a ./include/liblinear
endif

obj_all= $(objects) $(obj_linear)


.PHONY: all
all: $(target)

.PHONY: install
install: install3dec installtrain
	@echo "All work has been done."

.PHONY:install3dec
install3dec: 3Dec
	cp $(BINDIR)/3Dec $(DESTDIR);chmod +rx $(DESTDIR)/3Dec

.PHONY:installtrain
installtrain:3Dec-train 
	cp $(BINDIR)/3Dec-train $(DESTDIR);chmod +rx $(DESTDIR)/3Dec-train
	

.PHONY: uninstall
uninstall:
	rm $(DESTDIR)/3Dec $(DESTDIR)/3Dec-train

$(target): % : %.o $(objects) $(obj_linear)
	$(CXX) $(DEFINES) $(CFLAGS) $(CPPFLAGS) $(INCFLAGS) -o $(BINDIR)/$@   $(obj_all:%=$(BUILDDIR)/%)  $(BUILDDIR)/$(notdir $<) $(LDFLAGS)

$(objects) $(target:=.o): %.o: %.cpp $(header) Makefile
	$(CXX) $(DEFINES) $(CFLAGS) $(CPPFLAGS) $(INCFLAGS) -o $(BUILDDIR)/$@ -c $<

.f.o:
	$(FC) -m64 -O3 -o $@ -c $<

.PHONY: clean
clean:
	rm -f $(target:%=$(BINDIR)/%) *~ *.o  *.obj  *.exe *.lib ./*.lis $(MANDIR)/*.1 $(MANDIR)/*.1.html $(BUILDDIR)/*
	$(MAKE) -C include/liblinear  clean
	find . -name "*.log" -exec rm {} \;

