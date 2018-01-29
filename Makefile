
OBJS = HartreeFock.o Utils.o main.o

CC = g++ -std=c++11
DEBUG = -ggdb
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

prog : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o HF

main.o : main.cpp HartreeFock.hpp
	$(CC) $(CFLAGS) main.cpp

HartreeFock.o : HartreeFock.hpp HartreeFock.cpp Utils.hpp
	$(CC) $(CFLAGS) HartreeFock.cpp

Utils.o : Utils.hpp Utils.cpp
	$(CC) $(CFLAGS) Utils.cpp
