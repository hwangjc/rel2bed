#makefile for rel_to_bed

debug_or_optimize = -O3

CXX = g++
CXXFLAGS = -Wall -Werror -pedantic --std=c++17 $(debug_or_optimize)

main: main.cpp utr.h
	$(CXX) $(CXXFLAGS) main.cpp utr.h -o rel_to_bed


