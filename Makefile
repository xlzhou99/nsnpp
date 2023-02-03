
FLIB =-L/export/home/xlzhou/sav_nsnppz/nsnpp/lib

FFLAGS = -O3 -fp-model strict -assume protect_parens
OBJS = main.o   input_var.o  	nltermc_ns.o   dlepsw.o 	leinitz.o  error.o 

%.o: %.f90
	ifort $(FFLAGS)  -c $<

a.out: $(OBJS)
	ifort  $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(FLIB)  -llegen  -llapack  -lblas

clean:
	/bin/rm -f *% *~ $(OBJS)

virgin: clean
	/bin/rm -f a.out
###
