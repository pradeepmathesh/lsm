CXXFLAGS = -I/home/pradeep/sc1/pauls_lib
LDFLAGS = -lm
all: r
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c -o $@ $<
lsm: lsm.o matrix.o bitmaplib.o paulslib.o
#	tput reset;
	$(CXX) $(LDFLAGS) -o $@ $^
bitmaplib.o: /home/pradeep/sc1/pauls_lib/bitmaplib.c
	$(CXX) $(CFLAGS) -c -o $@ $<
paulslib.o: /home/pradeep/sc1/pauls_lib/paulslib.c
	$(CXX) $(CFLAGS) -c -o $@ $<
c:
	rm *.o
r: lsm
	rm *.o;time ./lsm