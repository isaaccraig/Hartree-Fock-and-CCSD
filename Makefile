
OBJS = HartreeFock.o Read.o QuantumUtils.o main.o

CC = g++ -std=c++11
DEBUG = -ggdb
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

prog : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o HF

main.o : main.cpp HartreeFock.hpp
	$(CC) $(CFLAGS) main.cpp

HartreeFock.o : HartreeFock.cpp HartreeFock.hpp Read.hpp QuantumUtils.hpp
	$(CC) $(CFLAGS) HartreeFock.cpp

Read.o : Read.cpp Read.hpp QuantumUtils.hpp
	$(CC) $(CFLAGS) Read.cpp

QuantumUtils.o : QuantumUtils.cpp QuantumUtils.hpp
	$(CC) $(CFLAGS) QuantumUtils.cpp
