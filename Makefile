CC=g++
CFLAGS=-std=c++11 -lmgl
SRC=prog.cpp
NAME=prog

happy:
	$(CC) $(CFLAGS) $(SRC) -o $(NAME)

clean:
	rm $(NAME)

