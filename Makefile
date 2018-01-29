
OBJS = HartreeFock.o Read.o main.o

CC = g++ -std=c++11
DEBUG = -ggdb
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

prog : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o HF

main.o : main.cpp HartreeFock.hpp 
	$(CC) $(CFLAGS) main.cpp

HartreeFock.o : HartreeFock.hpp HartreeFock.cpp Read.hpp
	$(CC) $(CFLAGS) HartreeFock.cpp

Read.o : Read.hpp Read.cpp
	$(CC) $(CFLAGS) Read.cpp
