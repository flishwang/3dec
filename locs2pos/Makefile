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

CXX = g++
FC = gfortran
CFLAGS = -Wall -O3 -funroll-loops -DNDEBUG -std=c++98
INCFLAGS = 
DEFINES =
objects =  locs2txt locs2txt_all.sh
header =
CPPFLAGS =

.PHONY: default
default: $(objects)

locs2txt: locs2txt.cpp
	$(CXX) $(DEFINES) $(CFLAGS) $(CPPFLAGS) $(INCFLAGS) -o $@ $<;chmod +x $@

.PHONY: clean
clean:
	rm -f $(objects)

