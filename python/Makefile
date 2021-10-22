CFLAGS=		-g -Wall -O2 -Wc++-compat #-Wextra
CPPFLAGS=	-DHAVE_KALLOC
INCLUDES=
OBJS=		mapnd.o
PROG=		mapnd.so
LIBS=		-lm -lz -lpthread
FPIC =      -fPIC -shared

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) -fPIC $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

mapnd:mapnd.o
		@+make -C ../
		$(CC) $(CFLAGS) $^ -o $@ -L../ -lminimap2 $(LIBS)

mapnd.so:mapnd.o
		@+make -C ../ CPPFLAGS+=-fPIC
		$(CC) $(CFLAGS) $(FPIC) $^ -o $@ -L../ -lminimap2 $(LIBS)

clean:
	rm -fr mapnd mapnd.o mapnd.so