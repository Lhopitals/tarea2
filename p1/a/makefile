CXX = g++

CXXFLAGS = -Wall -std=c++17 -Og -g

target = miprograma.exe

SRCS = $(wildcard *.cxx)

OBJS = $(patsubst %.cxx, %.o, $(SRCS))

FLGS = -Wall

$(target): $(OBJS)
		$(CXX) -o $@ $^
$(OBJS):%.o:%.cxx
		$(CXX) -c $< -o $@

clean:
	@rm -f *.o *~
borra:
	@rm -f *.o *.exe *~
run:
	@./$(target)
clear:
	@rm -f *.o *~
