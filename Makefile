CC=			gcc
CXX=		g++
#CFLAGS=		-Wall -g -O2 -static -funroll-loops -march=nocona -maccumulate-outgoing-args -fno-stack-protector
CFLAGS=		-Wall -g -funroll-loops -march=nocona -maccumulate-outgoing-args -fno-stack-protector

CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DHAVE_PTHREAD #-D_NO_SSE2 #-D_FILE_OFFSET_BITS=64	
OBJS =      bwtaln.o bwtgap.o BWT.o BWTConstruct.o utils.o \
            dictionary.o DNACount.o HSP.o HSPstatistic.o iniparser.o \
            inistrlib.o karlin.o MemManager.o MiscUtilities.o QSufSort.o \
            Socket.o 2BWT-Builder.o TextConverter.o Timing.o bamlite.o \
			2BWT-Interface.o bwaseqio.o r250.o cs2nt.o bwtse.o kstring.o \
			stdaln.o bwt_array.o

PROG=		bwsplice
INCLUDES=	
LIBS=		-lm -lz -lpthread
#SUBDIRS=	. bwt_gen

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

bwsplice:$(OBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) main.o $(LIBS) -o $@ 

bwtaln.o:BWT.h bwtaln.h kseq.h
bwtse.o:stdaln.h bwtaln.h utils.h kstring.h 2BWT-Interface.h
bwaseqio.o:bwaseqio.c bwtaln.h utils.h bamlite.h
bwt_array.o:bwt_array.h bwt_array.c
bamlite.o:bamlite.h
bwtgap.o:bwtgap.h bwtaln.h 2BWT-Interface.h bwt_array.h
#fastmap.o:BWT.h
BWT.o:	BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h
BWTConstruct.o:	BWTConstruct.h TypeNLimit.h BWT.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h QSufSort.h r250.h
dictionary.o:	dictionary.h inistrlib.h dictionary.h
DNACount.o:	DNACount.h TypeNLimit.h MiscUtilities.h
HSP.o:	TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h HSP.h HSPstatistic.h
HSPstatistic.o:	karlin.h HSPstatistic.h
iniparser.o:	iniparser.h dictionary.h inistrlib.h
inistrlib.o:	inistrlib.h
r250.o:	r250.h
karlin.o:	karlin.h
MemManager.o:	MiscUtilities.h TypeNLimit.h MemManager.h
MiscUtilities.o:	MiscUtilities.h TypeNLimit.h
QSufSort.o:	QSufSort.h TypeNLimit.h MiscUtilities.h
Socket.o:	Socket.h TypeNLimit.h MemManager.h MiscUtilities.h
2BWT-Builder.o:	2BWT-Builder.c
TextConverter.o:	TextConverter.h TypeNLimit.h MemManager.h MiscUtilities.h r250.h
2BWT-Interface.o:   2BWT-Interface.h
Timing.o:	Timing.h
clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a
