CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)
	@echo "*********************************************** "
	@echo "********* BUILDING seq-gen ******************** "
	@echo "*********************************************** "	
	make seq-gen
	
seq-gen:
	echo $(CWD)	
	tar -zxvf ../inst/ext/Seq-Gen.v1.3.3.tgz
	make -C Seq-Gen.v1.3.3/source 
	mv $(CWD)/Seq-Gen.v1.3.3/source/seq-gen ../inst/ext/seq-gen
	rm -rf Seq-Gen.v1.3.3	
