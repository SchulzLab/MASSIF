CC = g++ 
CFLAGS = -std=c++11 -fopenmp 
LDFLAGS = -static



#sources
SRC =  	PASTAA.cpp \
	PSCM_to_PSEM.cpp \
	TRAP.cpp \
	Matrix_new.cpp \
	TFtoMotifDomainInfo.cpp 

	

all: $(SRC) Makefile
	$(CC) PASTAA.cpp $(CFLAGS) $(LDFLAGS) -o PASTAA
	$(CC) PSCM_to_PSEM.cpp $(CFLAGS) $(LDFLAGS) -o PSCM_to_PSEM
	$(CC) TRAP.cpp $(CFLAGS) $(LDFLAGS)  -o TRAP 
	$(CC) Matrix_new.cpp $(CFLAGS) $(LDFLAGS) -o Marix
	$(CC) TFtoMotifDomainInfo.cpp $(CFLAGS) $(LDFLAGS) -o TFtoMotifDomainInfo
	

clean:
	-/bin/rm -f core *.o *~

