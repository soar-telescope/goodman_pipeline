# Copyrihgt (C) 2015 by Wojtek Pych, CAMK PAN
# 2015.07.10 version 2.2
#
# the C compiler
CC = gcc
# options for the C compiler
COPTS = -Wall -O3 -pedantic
# libraries
LIBS = -lm
#
NAME = dcr
#
${NAME} : ${NAME}.o swap.o pfitshead.o pfitsin.o pfitsio.o
	$(CC) -o ${NAME} $^ $(LIBS)
#
%.o: $(STDH) %.c
	$(CC) $(COPTS) -c $^
