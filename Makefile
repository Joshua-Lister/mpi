CXX = mpic++
CXXFLAGS = -Wall -O3
DEPENDICIES = $(wildcard include/*.h)

all: run_parallel

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -Iinclude -c $< -o $@
run_parallel: src/interface.o src/algo.o src/main.o $(DEPENDICIES)
	$(CXX) $(CXXFLAGS) -Iinclude src/interface.o src/algo.o src/main.o -o run_parallel

clean:
	rm run_parallel
