#   Makefile

CC = g++ 

CFLAGS = -std=c++11 -O3

all: MMult0 sor

sor: sor.cpp
	$(CC) $(CFLAGS) -o sor sor.cpp

MMult0: MMult0.cpp
	$(CC) $(CFLAGS) -o MMult0 MMult0.cpp

clean :
	rm MMult0 sor
