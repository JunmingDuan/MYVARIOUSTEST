source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
LDFLAGS +=-L/  -O3 -lpthread

all : main

%.o : %.cpp
	g++ -c -o $@ $< $(CXXFLAGS)

main : $(object)
	g++ -o $@ $(object) $(LDFLAGS) $(LIBS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.txt
	-rm -f *~

.PHONY : default clean
