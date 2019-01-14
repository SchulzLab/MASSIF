#sources
SRC =  	PASTAA.cpp \
	PSCM_to_PSEM.cpp \
	TRAP.cpp \
	Matrix_new.cpp \
	TFtoMotifDomainInfo.cpp 

platform=$(shell uname -s)

ifeq ($(platform), Darwin)

CC = clang++
CFLAGS =  -std=c++11  -Xpreprocessor -fopenmp
#CFLAGS =  -std=c++11  -stdlib=libc++ 
LDFLAGS = -stdlib=libc++
all: $(SRC) Makefile
	$(CC) PASTAA.cpp $(CFLAGS) $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp $(CFLAGS) $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP.cpp $(CFLAGS) $(LDFLAGS)  -o TRAP -lomp 
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Marix
	$(CC) TFtoMotifDomainInfo.cpp $(CFLAGS) $(LDFLAGS) -o TFtoMotifDomainInfo
else

CC = g++ 
CFLAGS = -std=c++11 -fopenmp 
LDFLAGS = -static
all: $(SRC) Makefile
	$(CC) PASTAA.cpp $(CFLAGS) $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp $(CFLAGS) $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP.cpp $(CFLAGS) $(LDFLAGS)  -o TRAP 
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Marix	
	$(CC) TFtoMotifDomainInfo.cpp $(CFLAGS) $(LDFLAGS) -o TFtoMotifDomainInfo

endif

clean:
	-/bin/rm -f core *.o *~
