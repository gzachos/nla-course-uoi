
CC = gcc
CFLAGS = -g -O2 -Wall -Wundef
LDLIBS = -lm 
OBJECTS =

.PHONY: cholesky cholesky-optimal clean
all: cholesky cholesky-optimal

cholesky: cholesky.c
	$(CC) $(CFLAGS) -UOPTIMIZED $^ -o cholesky $(LDLIBS)

cholesky-optimal: cholesky.c
	$(CC) $(CFLAGS) -DOPTIMIZED $^ -o cholesky-optimal $(LDLIBS)

clean:
	rm -rf cholesky cholesky-optimal


