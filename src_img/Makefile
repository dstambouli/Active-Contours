CC=mpicc 

DEBUG=-g -Wall -lm -pedantic -std=c99


OPTIM=-O3



CFLAGS= $(DEBUG) $(OPTIM) 


EXE = cv
SRC = main.c 

.SUFFIXES: .c .o

OBJS =  $(SRC:.c=.o)

all: $(EXE)
 
.c.o :
	$(CC)  -c $< $(CFLAGS)
 
$(EXE): $(OBJS)
	$(CC) -o $(EXE) $(OBJS) $(CFLAGS) $(LDFLAGS) 

clean: 
	rm -f $(OBJS) $(EXE)

mrproper: clean
	 rm -f $(EXE)

