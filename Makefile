CXX      ?= g++
CXXFLAGS ?= -std=c++17
CPPFLAGS ?= -O3 -Wall -pedantic -I. -I${mkEigenInc} -I/home/simoripa96/Scrivania/pacs-examples/Examples/include -isystem$(mkBoostInc)
LDFLAGS  ?= -L/home/simoripa96/Scrivania/pacs-examples/Examples/lib -Wl,-rpath=/home/simoripa96/Scrivania/pacs-examples/Examples/lib -L$(mkBoostLib)

LDLIBS   ?= -lmuparser -lboost_iostreams -lboost_system -lboost_filesystem
LINK.o := $(LINK.cc) # Use C++ linker.

DEPEND = make.dep

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

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)
