CC     = g++
CCFLAG = -march=native -O3 -pipe -Wall -Wextra -Wunused-parameter \
                        -Wunused-variable
INC    = -I$(METIS_DIR)/include
LDFLAGS = -lm -lmetis -L$(METIS_DIR)/lib
LDFLAGS += -Wl,-rpath=$(METIS_DIR)/lib


HDR = $(wildcard *.h)
SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,%.o,$(SRC))

PROGRAM  = graphpart

all:     $(PROGRAM)

.cc.o: $(HDR)
	$(CC) $(CCFLAG) $(INC) -c $*.cc

$(PROGRAM): $(OBJ)
	$(CC) -o $(PROGRAM) $(OBJ) $(LDFLAGS)


plot:
	python3 plot.py

clean:
	rm -rf $(OBJ) $(PROGRAM) point
