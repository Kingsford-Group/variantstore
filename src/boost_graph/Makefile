CC=g++
CFLAGS=-std=c++11 -Wall -lssl -lcrypto
LIBPATH=-I/../../include
FILE=boost_graph
all: $(FILE) clean
$(FILE): $(FILE).o
	$(CC) $(CFLAGS) $(LIBPATH) -o $(FILE) $(FILE).o 
$(FILE).o: $(FILE).cc
	$(CC) $(CFLAGS) $(LIBPATH) -c $(FILE).cc
clean:
	rm -f *.o
