.NOTPARALLEL:
.PHONY: all compile mer31 mer63 mer127 install clean

all:		compile install
	
compile:	mer31 mer63 mer127

mer31:
		@cd ./src/31mer; make
mer63:
		@cd ./src/63mer; make
mer127:
		@cd ./src/127mer; make
install:
		@cp ./src/*/SOAPdenovo* ./bin/
		@printf "Installation done.\n"
clean:
		@cd ./src/31mer; make clean
		@cd ./src/63mer; make clean
		@cd ./src/127mer; make clean
		@cd ./bin; rm -f SOAPdenovo*

