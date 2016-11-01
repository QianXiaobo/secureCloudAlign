# Makefile
#

OBJECTS = 

CXX = g++
#CXX = clang++

CXXFLAGS = -std=c++11

ifeq ($(dbg), 1)
	CXXFLAGS += -g3 -ggdb -W -Wall -Wextra -pedantic -Wno-variadic-macros -Wno-long-long 
else
	CXXFLAGS += -O3
endif

ifeq ($(pg), 1)
	CXXFLAGS += -pg
endif

INC = -I /home/yongzhao/include

LIB = -L /home/yongzhao/lib/

all: build encrypt decrypt exdare

build: build.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o build build.cpp $(OBJECTS) $(INC) $(LIB) -lsdsl -ldivsufsort -ldivsufsort64 -lboost_serialization -lboost_system -lboost_thread

encrypt: encrypt.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o encrypt encrypt.cpp $(OBJECTS) $(INC) $(LIB) -lpthread -lcryptopp

decrypt: decrypt.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o decrypt decrypt.cpp $(OBJECTS) $(INC) $(LIB) -lpthread -lcryptopp

exdare: exdare.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o exdare exdare.cpp $(OBJECTS) $(INC) $(LIB) -lpthread -lcryptopp -lboost_serialization -lboost_system -lboost_thread

index: fm_sdsl.cpp $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o index fm_sdsl.cpp $(OBJECTS) $(INC) $(LIB) -lsdsl -ldivsufsort -ldivsufsort64 -lboost_system -lboost_thread -lboost_serialization

$(OBJECTS): %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(INC)

clean:
	rm *.o *~ build encrypt decrypt exdare
