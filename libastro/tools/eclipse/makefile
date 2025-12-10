
APPNAME = inex_list

INCASTRO = ../../..
LIBASTRO = ../..

CC = gcc

COPTIONS = -Wall -ansi -pedantic -std=c99 -O2

INCLUDES = -I$(INCASTRO)
CFLAGS = $(COPTIONS) $(INCLUDES)

LDFLAGS = -L$(LIBASTRO)
LIBS = -lastro -lm

INEXDATAFILENAME = data/total.solar.eclipse.saros.inex.list.txt
ECLIPSEDATAFILENAME = data/total.solar.eclipse.list.txt
SHA512FILENAME = data/SHA512SUMS.eclipse


all: $(APPNAME)

clean:
	rm -f $(APPNAME)

cleantest:
	rm -f $(INEXDATAFILENAME)
	rm -f $(ECLIPSEDATAFILENAME)

$(INEXDATAFILENAME): $(APPNAME)
	./$(APPNAME) 0 180 > $(INEXDATAFILENAME)

$(ECLIPSEDATAFILENAME): $(INEXDATAFILENAME)
	grep . $(INEXDATAFILENAME) | sort -n -k 4 > $(ECLIPSEDATAFILENAME)

data: $(ECLIPSEDATAFILENAME)

hash:
	sha512sum $(INEXDATAFILENAME) $(ECLIPSEDATAFILENAME) > $(SHA512FILENAME)

test: data
	@sha512sum -c $(SHA512FILENAME)

$(APPNAME): $(APPNAME).c ../../libastro.a ../../*.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(APPNAME) $(APPNAME).c $(LIBS)
