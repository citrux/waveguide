CC=g++
CFLAGS=-std=c++11 -lmgl
SRC=dispersion_relation.cpp
NAME=dispersion_relation

happy:
	$(CC) $(CFLAGS) $(SRC) -o $(NAME)

clean:
	rm $(NAME)

