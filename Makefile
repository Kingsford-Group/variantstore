TARGETS=variantstore test_graphcontainer test_variantgraph test_index test_query bm_query

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
					 -lbz2 -llzma -lpthread -lsdsl \
						`pkg-config --cflags --libs protobuf` -ltcmalloc

#
# declaration of dependencies
#

all: variantstore
tests:	test_graphcontainer test_variantgraph test_index test_query
bm: bm_query

graph:	$(SER)/graph.png

# dependencies between programs and .o files
variantstore:						$(OBJDIR)/variantstore.o $(OBJDIR)/variantstore_fs.o \
												$(OBJDIR)/commands.o \
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_graphcontainer:		$(OBJDIR)/test_graphcontainer.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_variantgraph:			$(OBJDIR)/test_variantgraph.o  $(OBJDIR)/variantstore_fs.o\
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_index:							$(OBJDIR)/test_index.o  $(OBJDIR)/variantstore_fs.o\
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
test_query:							$(OBJDIR)/test_query.o  $(OBJDIR)/variantstore_fs.o\
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o
bm_query:								$(OBJDIR)/bm_query.o  $(OBJDIR)/variantstore_fs.o\
												$(OBJDIR)/variantgraphvertex.pb.o $(OBJDIR)/gqf.o \
												$(OBJDIR)/gqf_file.o $(OBJDIR)/rank_select.o \
												$(OBJDIR)/hashutil.o $(OBJDIR)/util.o

# dependencies between .o files and .cc (or .c) files
$(OBJDIR)/variantstore.o: 				$(LOC_SRC)/variantstore.cc 
$(OBJDIR)/commands.o: 						$(LOC_SRC)/commands.cc \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/query.h \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/variantstore_fs.h \
																	$(LOC_INCLUDE)/graph.h \
																	$(LOC_INCLUDE)/stream.hpp
$(OBJDIR)/test_graphcontainer.o: 	$(LOC_SRC)/test_graphcontainer.cc \
																	$(LOC_INCLUDE)/gqf_cpp.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/test_variantgraph.o: 		$(LOC_SRC)/test_variantgraph.cc \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/variantstore_fs.h \
																	$(LOC_INCLUDE)/graph.h \
																	$(LOC_INCLUDE)/dot_graph.h \
																	$(LOC_INCLUDE)/stream.hpp
$(OBJDIR)/test_index.o: 					$(LOC_SRC)/test_index.cc \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/query.h \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/variantstore_fs.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/test_query.o: 					$(LOC_SRC)/test_query.cc \
																	$(LOC_INCLUDE)/query.h \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/variantstore_fs.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/bm_query.o: 						$(LOC_SRC)/bm_query.cc \
																	$(LOC_INCLUDE)/index.h \
																	$(LOC_INCLUDE)/variant_graph.h \
																	$(LOC_INCLUDE)/variantgraphvertex.pb.h \
																	$(LOC_INCLUDE)/variantstore_fs.h \
																	$(LOC_INCLUDE)/graph.h
$(OBJDIR)/variantgraphvertex.pb.o: 	$(LOC_SRC)/variantgraphvertex.pb.cc \
																		$(LOC_INCLUDE)/variantgraphvertex.pb.h
$(OBJDIR)/gqf.o: 				$(LOC_SRC)/gqf/gqf.c $(LOC_INCLUDE)/gqf/gqf.h
$(OBJDIR)/gqf_file.o: 	$(LOC_SRC)/gqf/gqf_file.c $(LOC_INCLUDE)/gqf/gqf_file.h
$(OBJDIR)/hashutil.o: 	$(LOC_INCLUDE)/gqf/hashutil.h
$(OBJDIR)/rank_select.o: 	$(LOC_INCLUDE)/gqf/rank_select.h


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
