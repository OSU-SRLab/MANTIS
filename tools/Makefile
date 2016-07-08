CC=g++
CFLAGS=-O3
LDFLAGS=-O3

all: RepeatFinder

RepeatFinder: main.o circularLinkedList.o
	$(CC) $(LDFLAGS) main.o circularLinkedList.o -o RepeatFinder
	strip RepeatFinder

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

circularLinkedList.o: circularLinkedList.cpp circularLinkedList.hpp
	$(CC) $(CFLAGS) -c circularLinkedList.cpp -o circularLinkedList.o

clean:
	rm -f *.o RepeatFinder

