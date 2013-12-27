#------------------------------------------------------------------------------#
# This makefile was generated by 'cbp2make' tool rev.147                       #
#------------------------------------------------------------------------------#


WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -Wall
RESINC = 
LIBDIR = -L/usr/lib
LIB = 
LDFLAGS = -lgsl -lgslcblas -lboost_program_options

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -g
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/Simul_Prog

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O2
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS) -s
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/Simul_Prog

OBJ_DEBUG = $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/Random.o $(OBJDIR_DEBUG)/Population.o $(OBJDIR_DEBUG)/Parameters.o $(OBJDIR_DEBUG)/OutputFormat.o $(OBJDIR_DEBUG)/Individual.o $(OBJDIR_DEBUG)/Haplotype.o $(OBJDIR_DEBUG)/Allele.o $(OBJDIR_DEBUG)/Genotype.o $(OBJDIR_DEBUG)/GeneticMap.o $(OBJDIR_DEBUG)/Fitness.o $(OBJDIR_DEBUG)/Environment.o $(OBJDIR_DEBUG)/Architecture.o $(OBJDIR_DEBUG)/ArchiMultilinear.o $(OBJDIR_DEBUG)/ArchiAdditive.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/Random.o $(OBJDIR_RELEASE)/Population.o $(OBJDIR_RELEASE)/Parameters.o $(OBJDIR_RELEASE)/OutputFormat.o $(OBJDIR_RELEASE)/Individual.o $(OBJDIR_RELEASE)/Haplotype.o $(OBJDIR_RELEASE)/Allele.o $(OBJDIR_RELEASE)/Genotype.o $(OBJDIR_RELEASE)/GeneticMap.o $(OBJDIR_RELEASE)/Fitness.o $(OBJDIR_RELEASE)/Environment.o $(OBJDIR_RELEASE)/Architecture.o $(OBJDIR_RELEASE)/ArchiMultilinear.o $(OBJDIR_RELEASE)/ArchiAdditive.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

$(OBJDIR_DEBUG)/Random.o: Random.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Random.cpp -o $(OBJDIR_DEBUG)/Random.o

$(OBJDIR_DEBUG)/Population.o: Population.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Population.cpp -o $(OBJDIR_DEBUG)/Population.o

$(OBJDIR_DEBUG)/Parameters.o: Parameters.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Parameters.cpp -o $(OBJDIR_DEBUG)/Parameters.o

$(OBJDIR_DEBUG)/OutputFormat.o: OutputFormat.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c OutputFormat.cpp -o $(OBJDIR_DEBUG)/OutputFormat.o

$(OBJDIR_DEBUG)/Individual.o: Individual.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Individual.cpp -o $(OBJDIR_DEBUG)/Individual.o

$(OBJDIR_DEBUG)/Haplotype.o: Haplotype.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Haplotype.cpp -o $(OBJDIR_DEBUG)/Haplotype.o

$(OBJDIR_DEBUG)/Allele.o: Allele.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Allele.cpp -o $(OBJDIR_DEBUG)/Allele.o

$(OBJDIR_DEBUG)/Genotype.o: Genotype.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Genotype.cpp -o $(OBJDIR_DEBUG)/Genotype.o

$(OBJDIR_DEBUG)/GeneticMap.o: GeneticMap.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c GeneticMap.cpp -o $(OBJDIR_DEBUG)/GeneticMap.o

$(OBJDIR_DEBUG)/Fitness.o: Fitness.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Fitness.cpp -o $(OBJDIR_DEBUG)/Fitness.o

$(OBJDIR_DEBUG)/Environment.o: Environment.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Environment.cpp -o $(OBJDIR_DEBUG)/Environment.o

$(OBJDIR_DEBUG)/Architecture.o: Architecture.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c Architecture.cpp -o $(OBJDIR_DEBUG)/Architecture.o

$(OBJDIR_DEBUG)/ArchiMultilinear.o: ArchiMultilinear.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c ArchiMultilinear.cpp -o $(OBJDIR_DEBUG)/ArchiMultilinear.o

$(OBJDIR_DEBUG)/ArchiAdditive.o: ArchiAdditive.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c ArchiAdditive.cpp -o $(OBJDIR_DEBUG)/ArchiAdditive.o

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/Random.o: Random.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Random.cpp -o $(OBJDIR_RELEASE)/Random.o

$(OBJDIR_RELEASE)/Population.o: Population.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Population.cpp -o $(OBJDIR_RELEASE)/Population.o

$(OBJDIR_RELEASE)/Parameters.o: Parameters.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Parameters.cpp -o $(OBJDIR_RELEASE)/Parameters.o

$(OBJDIR_RELEASE)/OutputFormat.o: OutputFormat.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c OutputFormat.cpp -o $(OBJDIR_RELEASE)/OutputFormat.o

$(OBJDIR_RELEASE)/Individual.o: Individual.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Individual.cpp -o $(OBJDIR_RELEASE)/Individual.o

$(OBJDIR_RELEASE)/Haplotype.o: Haplotype.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Haplotype.cpp -o $(OBJDIR_RELEASE)/Haplotype.o

$(OBJDIR_RELEASE)/Allele.o: Allele.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Allele.cpp -o $(OBJDIR_RELEASE)/Allele.o

$(OBJDIR_RELEASE)/Genotype.o: Genotype.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Genotype.cpp -o $(OBJDIR_RELEASE)/Genotype.o

$(OBJDIR_RELEASE)/GeneticMap.o: GeneticMap.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c GeneticMap.cpp -o $(OBJDIR_RELEASE)/GeneticMap.o

$(OBJDIR_RELEASE)/Fitness.o: Fitness.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Fitness.cpp -o $(OBJDIR_RELEASE)/Fitness.o

$(OBJDIR_RELEASE)/Environment.o: Environment.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Environment.cpp -o $(OBJDIR_RELEASE)/Environment.o

$(OBJDIR_RELEASE)/Architecture.o: Architecture.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c Architecture.cpp -o $(OBJDIR_RELEASE)/Architecture.o

$(OBJDIR_RELEASE)/ArchiMultilinear.o: ArchiMultilinear.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ArchiMultilinear.cpp -o $(OBJDIR_RELEASE)/ArchiMultilinear.o

$(OBJDIR_RELEASE)/ArchiAdditive.o: ArchiAdditive.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ArchiAdditive.cpp -o $(OBJDIR_RELEASE)/ArchiAdditive.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

