TARGETS= main generate_seed

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

CXX = g++ -std=c++11
CC = gcc -std=gnu11
LD= g++ -std=c++11

LOC_INCLUDE=include
LOC_SRC=src
OBJDIR=obj

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -lpthread -lboost_system \
-lboost_thread -lm -lbz2 -lz -lrt

#
# declaration of dependencies
#

all: $(TARGETS)

# dependencies between programs and .o files
main:							$(OBJDIR)/main.o $(OBJDIR)/sketch.o \
									$(OBJDIR)/marsaglia.o $(OBJDIR)/misc.o \
									$(OBJDIR)/fsutil.o

generate_seed: $(OBJDIR)/generate_seed.o

# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/sketch.o: 		$(LOC_SRC)/sketch.cc  $(LOC_INCLUDE)/sketch.h \
												$(LOC_INCLUDE)/marsaglia.hpp \
												$(LOC_INCLUDE)/seeded_prg.hpp $(LOC_INCLUDE)/misc.hpp
$(OBJDIR)/misc.o: 			$(LOC_SRC)/misc.cc $(LOC_INCLUDE)/misc.hpp
$(OBJDIR)/marsaglia.o: 	$(LOC_SRC)/marsaglia.cc $(LOC_INCLUDE)/marsaglia.hpp
$(OBJDIR)/fsutil.o: 		$(LOC_SRC)/fsutil.cc $(LOC_INCLUDE)/fsutil.h
$(OBJDIR)/main.o: 			$(LOC_SRC)/main.cc \
												$(LOC_INCLUDE)/marsaglia.hpp \
												$(LOC_INCLUDE)/seeded_prg.hpp $(LOC_INCLUDE)/misc.hpp
$(OBJDIR)/generate_seed.o: $(LOC_SRC)/generate_seed.cc $(LOC_INCLUDE)/seeded_prg.hpp

#
# generic build rules
#

$(TARGETS):
	$(LD) $^ $(LDFLAGS) -o $@

$(OBJDIR)/%.o: $(LOC_SRC)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR)/%.o: $(LOC_SRC)/%.c | $(OBJDIR)
	$(CXX) $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) core $(TARGETS)

