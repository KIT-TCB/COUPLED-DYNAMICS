#CXX= g++
#CXX=icpc
CXXOPTS= -O2 -Wall -pedantic -g
OBJ= do_rksuite.o rksuite.o

librksuite: $(OBJ)
	ar rvs librksuite.a $(OBJ)
%.o: %.cpp
	$(CXX) $(CXXOPTS) -c $<
clean:
	rm -f librksuite.a $(OBJ)
