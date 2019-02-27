TARGETS= main

ifdef D
	DEBUG=-g
	OPT=
else
	DEBUG=
	OPT=-Ofast
endif

ifdef NH
	ARCH=
else
	ARCH=-msse4.2 -D__SSE4_2_
endif

ifdef P
	PROFILE=-pg -no-pie # for bug in gprof.
endif

CXX = g++ -std=c++17
CC = gcc -std=gnu11
LD= g++ -std=c++17

LOC_INCLUDE=include
LOC_LIB=lib
LOC_SRC=src
OBJDIR=obj

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE) 
						

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -L$(LOC_LIB) -lm -lvcflib -lhts -lz \
					 -lbz2 -llzma -lrt -lpthread -lssl -lcrypto

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
main:										$(OBJDIR)/main.o $(OBJDIR)/gqf.o $(OBJDIR)/gqf_file.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o

# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/main.o: 			$(LOC_SRC)/main.cc $(LOC_INCLUDE)/gqf_cpp.h \
												$(LOC_INCLUDE)/graph.h

$(OBJDIR)/gqf.o: 				$(LOC_SRC)/gqf/gqf.c $(LOC_INCLUDE)/gqf/gqf.h
$(OBJDIR)/gqf_file.o: 	$(LOC_SRC)/gqf/gqf_file.c $(LOC_INCLUDE)/gqf/gqf_file.h
$(OBJDIR)/hashutil.o: 	$(LOC_INCLUDE)/gqf/hashutil.h

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/gqf/%.c | $(OBJDIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) core $(TARGETS)

