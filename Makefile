CXX      ?= g++
CXXFLAGS ?= -std=c++17
CPPFLAGS ?= -O3 -Wall -pedantic -I. -I${mkEigenInc} -I/home/simoripa96/Scrivania/pacs-examples/Examples/include -isystem$(mkBoostInc) #Cambiare i percorsi
LDFLAGS  ?= -L/home/simoripa96/Scrivania/pacs-examples/Examples/lib -Wl,-rpath=/home/simoripa96/Scrivania/pacs-examples/Examples/lib -L$(mkBoostLib) #Cambiare i percorsi

LDLIBS   ?= -lmuparser -lboost_iostreams -lboost_system -lboost_filesystem
LINK.o := $(LINK.cc) # Use C++ linker.


EXEC = main
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)


.PHONY = all $(EXEC) $(OBJS) clean distclean $(DEPEND)

all: $(DEPEND) $(EXEC)

$(EXEC): $(OBJS)



$(OBJS): %.o: %.cpp

clean:
	$(RM) $(DEPEND)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~


