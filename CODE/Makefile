ARCH		= $(DOW_MACH)

SHELL		= /bin/sh
#.SILENT:

DOWSER_DIR	= .
BINDIR = ../bin/$(ARCH)

COMMANDS	= $(DOWSER_DIR)/Makearch.$(ARCH)
FILELIST	= $(DOWSER_DIR)/Makedata.files
DEPENDFILE	= $(DOWSER_DIR)/Makedata.depend

include	$(COMMANDS)
include $(FILELIST)

CFLAGS  = $(ARCH_COPTS)
FFLAGS  = $(ARCH_FOPTS)

LFLAGS	= $(ARCH_LOPTS)
FLFLAGS	= $(ARCH_FLOPTS)

.c.o:
	$(CC) $(CFLAGS) -c $(DOWSER_DIR)/$<
.f.o:
	$(F77) $(FFLAGS) -c $(DOWSER_DIR)/$<

all: $(ARCH) bin/$(ARCH) bin/dowserx bin/dowser-repeat
	cd $(ARCH) ; $(MAKE) -f ../Makefile.arch
	$(MAKE) $(DOWSER_EXEC) 
	@$(ECHO) have made all dowser parts

$(BINDIR)/reformatPDB: $(ARCH)/reformatPDB.o $(ARCH)/readPDB10.o $(ARCH)/AddAtoms.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/reformatPDB.o $(ARCH)/readPDB10.o $(ARCH)/AddAtoms.o $(LFLAGS)

$(BINDIR)/chooser: $(ARCH)/chooser.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/chooser.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(LFLAGS)

$(BINDIR)/drain: $(ARCH)/drain.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/drain.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(LFLAGS)

$(BINDIR)/scrape: $(ARCH)/scrape.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/scrape.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(LFLAGS)

$(BINDIR)/pdb2ms: $(ARCH)/pdb2ms.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/pdb2ms.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(LFLAGS)

$(BINDIR)/xms: $(ARCH)/xms.o
	$(F77) $(FFLAGS) -o $@ $(ARCH)/xms.o $(FLFLAGS)

$(BINDIR)/qms: $(ARCH)/qms.o $(ARCH)/MakeApex.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/qms.o $(ARCH)/MakeApex.o $(LFLAGS)

$(BINDIR)/placeWat: $(ARCH)/placeWat.o $(ARCH)/readPDB10.o $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/placeWat.o $(ARCH)/Align.o $(ARCH)/MinimizeFRPR.o $(ARCH)/readPDB10.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/RefineWat: $(ARCH)/RefineWat.o $(ARCH)/readPDB10.o $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/RefineWat.o $(ARCH)/Align.o  $(ARCH)/MinimizeFRPR.o $(ARCH)/readPDB10.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/CompareWat: $(ARCH)/CompareWat.o $(ARCH)/readPDB10.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/CompareWat.o $(ARCH)/readPDB10.o $(LFLAGS)

$(BINDIR)/SetOH: $(ARCH)/SetOH.o $(ARCH)/readPDB10.o $(ARCH)/Align.o $(ARCH)/WriteStr.o $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/SetOH.o $(ARCH)/Align.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/AmideFlip: $(ARCH)/AmideFlip.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Align.o  $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/AmideFlip.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Align.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/FindHBond: $(ARCH)/FindHBond.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/FindHBond.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/ReduceIntSurf: $(ARCH)/ReduceIntSurf.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Math.o
	$(CC) $(CFLAGS) -o $@ $(ARCH)/ReduceIntSurf.o $(ARCH)/readPDB10.o $(ARCH)/WriteStr.o $(ARCH)/Math.o $(LFLAGS)

$(BINDIR)/DictMake: $(ARCH)/DictMake.o $(ARCH)/readPDB10.o 
	$(CC) $(CFLAGS) -o $@ $(ARCH)/DictMake.o $(ARCH)/readPDB10.o $(LFLAGS)

clean:
	@if test -d $(ARCH) ; then cd $(ARCH) ; $(MAKE) -f ../Makefile.arch clean ; fi ; true
	@$(DELETE) $(DOWSER_EXEC) ../bin/dowserx ../bin/dowser-repeat
	@/bin/rm -fR $(ARCH) ../bin/$(ARCH)
	@echo "     Remove object files and executables for $(ARCH) architecture"

rmdata:
	$(DELETE) core

# include $(DEPENDFILE)

$(ARCH):
	@if test ! -d $(ARCH) ; then mkdir $(ARCH) ; fi ; true

bin/$(ARCH):
	@if test ! -d ../bin/$(ARCH) ; then mkdir ../bin/$(ARCH) ; fi ; true

bin/dowserx:
	@if test ! -r ../bin/dowserx ; then cd ../bin ; ln -s dowser dowserx; fi; true

bin/dowser-repeat:
	@if test ! -r ../bin/dowser-repeat ; then cd ../bin ; ln -s dowser dowser-repeat; fi; true
