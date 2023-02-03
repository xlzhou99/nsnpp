.SUFFIXES: .o .f 

####### FLAGS ############
DEST	      	= .
F77FLAGS	=-O2
LIBRARY	      	=liblegen.a
MAKEFILE	= Makefile
IFORTFLAGS	=-O3 -fp-model strict -assume protect_parens
####### COMPILE #########
COMPILE.ifort	=ifort $(IFORTFLAGS) -c -o

.f.o :
	$(COMPILE.ifort) $<

####### SOURCE #########

OBJS	      = dlel2no.o dleh1no.o dln.o dln2.o legcheb.o\
		dlestopx.o \
		dlestopxy.o \
		dnftos.o \
		ftos2.o \
		ftos2a.o \
		ftos2b.o \
		jacobl.o \
		ledeig.o \
		ledhmz1.o \
		ledhmz1b.o \
		ledhmz2.o \
		ledhmz2a.o \
		ledhmz2b.o \
		lednhmz.o \
		lednstof.o \
		ledx.o \
		ledxx.o \
		ledy.o \
		ledyy.o \
		leinit.o  leinitc.o\
		lendhmz.o \
		lendstof.o \
		leneig.o leneigw.o\
		lenhmz0.o lenhmz1w.o\
		lenhmz1.o \
		lenhmz1b.o \
		lenhmz2.o \
		lenstof2.o lenhmz2w.o lenhmz2wb.o lenhmz2wc.o\
		lepdhmz2b.o ledxc.o ledyc.o\
		leptos.o \
		leptos2.o \
		lestof2.o \
		lestof2b.o \
		lestop.o \
		lestop2.o \
		letouni1.o \
		letouni2.o letrfm1.o letrfm2.o \
		levalue2.o \
		ndftos.o \
		nftos2.o\
                lebih1db.o  lebih1d.o  leftosbih.o   \
                leinitbih.o leinitbihc.o lestofbih.o sympenta.o  \


########## MAKE LIBRARY#################
all:		$(LIBRARY)

$(LIBRARY): 	 $(OBJS)
		@echo -n "Loading $(LIBRARY) ... "
		@ar cr $(LIBRARY)  $(OBJS)
		@echo "done"
###
