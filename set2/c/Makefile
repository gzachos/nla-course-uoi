CC = gcc
CFLAGS = -g -O2 -Wall -Wundef
LDLIBS = -lm 
OBJECTS =

.PHONY: set2 set2-optimal clean
all: set2 set2-optimal

set2: set2.c
	$(CC) $(CFLAGS) -UOPTIMIZED $^ -o set2 $(LDLIBS)

set2-optimal: set2.c
	$(CC) $(CFLAGS) -DOPTIMIZED $^ -o set2-optimal $(LDLIBS)

clean:
	rm -rf set2 set2-optimal


