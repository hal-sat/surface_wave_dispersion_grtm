targ = main 
cc = latexmk
CFLAGS = -lualatex

all = $(targ) 

main:
	$(cc) $(CFLAGS) $(targ)

view:
	evince $(targ:=.pdf) &

clean:
	$(cc) -C

	
