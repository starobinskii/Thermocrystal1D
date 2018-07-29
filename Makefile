#
# 'make depend' uses makedepend to automatically generate dependencies 
#			   (dependencies are added to end of Makefile)
# 'make'		build executable file './Executions/task'
# 'make clean'  removes all .o and executable files
#

#*# ************************************************************************ #*#

CC = g++

CFLAGS = -w -std=c++11 -O3
CFLAGS += $(shell mpicxx -showme:compile)

INCLUDES = 

LFLAGS = $(shell mpicxx -showme:link)

LIBS = 

SRCS = ./Sources/main.cc ./Sources/thermocrystal1D.cc

MAIN = ./Executions/task

#*# ************************************************************************ #*#

OBJS = $(SRCS:.cc=.o)

.PHONY: depend clean

all:	$(MAIN)
		@echo  $(MAIN) has been compiled

$(MAIN): $(OBJS) 
		$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.cc.o:
		$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
		$(RM) ./Sources/*.o ./Sources/*~ $(MAIN)

depend: $(SRCS)
		makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE - make depend needs it