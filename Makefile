CXX ?= g++
BOOST_INSTALL_PATH ?= $(BOOST_ROOT)
OBABEL_INSTALL_PATH ?= $(OPEN_BABEL_ROOT)
CXXFLAGS = -std=c++11
ifeq ($(Debug),Y)
CXXFLAGS += -g -O0
else
CXXFLAGS += -O2
endif

ifeq ($(Symbols),Y)
CXXFLAGS += -g
endif

ifeq ($(OpenMP),Y)
OMPFLAG = -fopenmp
else
OMPFLAG =
endif

CXXFLAGS += $(OMPFLAG)

ifeq ($(BOOST_INSTALL_PATH),)
BOOSTIP =
BOOSTLP =
else
BOOSTIP = -I$(BOOST_INSTALL_PATH)/include
BOOSTLP = -L$(BOOST_INSTALL_PATH)/lib
endif
OBABELIP = -I$(OBABEL_INSTALL_PATH)/include/openbabel-2.0
OBABELLP = -L$(OBABEL_INSTALL_PATH)/lib

.SUFFIXES: .cc .o

_GRID_OBJS = grid_main.o main_utils.o utils.o infile_reader.o Molecule.o Fragment.o Vector3d.o Atom.o AtomInterEnergyGrid.o InterEnergyGrid.o EnergyCalculator.o log_writer_stream.o OBMol.o
_FRAGQUERY_OBJS = fragquery_main.o main_utils.o Vector3d.o InterEnergyGrid.o Molecule.o Fragment.o InterEnergyGrid.o EnergyCalculator.o infile_reader.o utils.o AtomInterEnergyGrid.o FragmentInterEnergyGrid.o Atom.o log_writer_stream.o OBMol.o QueryPriority.o
_OPENDX_OBJS = fragdx_main.o main_utils.o Vector3d.o InterEnergyGrid.o Molecule.o Fragment.o InterEnergyGrid.o EnergyCalculator.o infile_reader.o OpenDx.o utils.o AtomInterEnergyGrid.o FragmentInterEnergyGrid.o Atom.o log_writer_stream.o OBMol.o

# ALL = atomgrid-gen fragment-query opendx-gen
ALL = atomgrid-gen fragment-query

GRID_OBJS = $(patsubst %,objs/%,$(_GRID_OBJS))
FRAGQUERY_OBJS = $(patsubst %,objs/%,$(_FRAGQUERY_OBJS))
OPENDX_OBJS = $(patsubst %,objs/%,$(_OPENDX_OBJS))


all: objs $(ALL)
	rm -r objs

query: objs fragment-query
	rm -r objs

objs:
	mkdir -p objs
	mkdir -p objs/test

atomgrid-gen: $(GRID_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

fragment-query: $(FRAGQUERY_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel

opendx-gen: $(OPENDX_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(BOOSTLP) $(OBABELLP) -lboost_regex -lboost_program_options -lopenbabel


objs/%.o: src/%.cc src/%.hpp src/common.hpp
	$(CXX) $(CXXFLAGS) $(BOOSTIP) $(OBABELIP) -o $@ -c $<

objs/%.o: src/%.cc src/common.hpp
	$(CXX) $(CXXFLAGS) $(BOOSTIP) $(OBABELIP) -o $@ -c $<

clean:
	rm -rf objs $(ALL)

