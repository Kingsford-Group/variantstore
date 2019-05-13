TARGETS=variantdb test_graphcontainer test_variantgraph test_index

ifdef D
	DEBUG=-g -DDEBUG_MODE
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
dot= dot -Grankdir=LR -Tpng 

LOC_INCLUDE=include
LOC_LIB=lib
LOC_SRC=src
OBJDIR=obj
SER=ser

CXXFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -fopenmp -m64 -I. -I$(LOC_INCLUDE)

CFLAGS += -Wall $(DEBUG) $(PROFILE) $(OPT) $(ARCH) -m64 -I. -I$(LOC_INCLUDE)

LDFLAGS += $(DEBUG) $(PROFILE) $(OPT) -L$(LOC_LIB) -lm -lvcflib -lhts -lz \
					 -lbz2 -llzma -lrt -lpthread -lssl -lcrypto -lboost_system -lsdsl \
						`pkg-config --cflags --libs protobuf` -ltcmalloc 

#
# declaration of dependencies
#

all: variantdb
tests:	test_graphcontainer test_variantgraph test_index

graph:	$(SER)/graph.png

# dependencies between programs and .o files
test_graphcontainer:		$(OBJDIR)/test_graphcontainer.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_variantgraph:			$(OBJDIR)/test_variantgraph.o \
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_index:							$(OBJDIR)/test_index.o \
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
variantdb:							$(OBJDIR)/variantdb.o \
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/test_graphcontainer.o: 	$(LOC_SRC)/test_graphcontainer.cc \
																	$(LOC_INCLUDE)/gqf_cpp.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/test_variantgraph.o: 		$(LOC_SRC)/test_variantgraph.cc \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/graph.h \
																	$(LOC_INCLUDE)/dot_graph.h \
																	$(LOC_INCLUDE)/stream.hpp
$(OBJDIR)/test_index.o: 					$(LOC_SRC)/test_index.cc \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/variantdb.o: 						$(LOC_SRC)/variantdb.cc \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/graph.h \
																	$(LOC_INCLUDE)/stream.hpp

$(OBJDIR)/variantgraphvertex.pb.o: 	$(LOC_SRC)/variantgraphvertex.pb.cc \
																		$(LOC_INCLUDE)/variantgraphvertex.pb.h
$(OBJDIR)/gqf.o: 				$(LOC_SRC)/gqf/gqf.c $(LOC_INCLUDE)/gqf/gqf.h
$(OBJDIR)/gqf_file.o: 	$(LOC_SRC)/gqf/gqf_file.c $(LOC_INCLUDE)/gqf/gqf_file.h
$(OBJDIR)/hashutil.o: 	$(LOC_INCLUDE)/gqf/hashutil.h


# create dot graph
$(SER)/graph.png: $(SER)/graph.dot

# generate proto files
proto:
	protoc --cpp_out=./ ./include/variantgraphvertex.proto && \
		mv include/variantgraphvertex.pb.cc src

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

$(SER)/%.png: $(SER)/%.dot | $(SER)
	$(dot) -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(SER):
	@mkdir -p $(SER)

clean:
	rm -rf $(OBJDIR) $(SER)/graph.png core $(TARGETS)

