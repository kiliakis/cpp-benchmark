CC=g++
CFLAGS=-std=c++11 -Ofast -Wall
LDFLAGS=
HOME=/afs/cern.ch/user/k/kiliakis/opt/lib
all: sin-war sin-LUT

sin-war: sin-war.cpp
	$(CC) $(CFLAGS) -L${HOME} -lvdt -o $@ $<

sin-LUT: sin-LUT.cpp
	$(CC) $(CFLAGS) -o $@ $<

%.o : %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<


clean:
	rm sin-war sin-LUT *.o
