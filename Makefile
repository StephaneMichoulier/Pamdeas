
CC = g++
CXXFLAGS = -Wall -std=c++11

$(info $(shell mkdir -p build))

EXEC = build/eden
SRCS = source/main.cpp source/general_functions.cpp source/porosity.cpp source/evol.cpp source/readwrite.cpp source/disruption.cpp

DEPEND_FILE = .depend
OBJS = $(SRCS:%.cpp=%.o)

all: $(DEPEND_FILE) $(EXEC)

$(EXEC): $(OBJS)
	g++ $(OBJS) -o $(EXEC)

clean:
	$(RM) $(EXEC)
	$(RM) $(OBJS)

$(DEPEND_FILE): $(SRCS)
	g++ -MM $(SRCS) > $(DEPEND_FILE)

-include $(DEPEND_FILE)
