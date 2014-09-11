DEST= $(R_LIBS_USER)/rPANGEAHIVsim
EXTERNAL_CODE= $(DEST)/ext
MISC_DIR= $(DEST)/misc
CWD= $$(pwd)

all:	
	echo "current directory is "$(CWD)
	echo "destination directory for external code is "$(EXTERNAL_CODE)
	mkdir -p $(EXTERNAL_CODE)
	mkdir -p $(MISC_DIR)
	@echo "*********************************************** "
	@echo "********* BUILDING VirusTreeSimulator ********* "
	@echo "*********************************************** "	
	make VirusTreeSimulator.jar	
	@echo "*********************************************** "
	@echo "********* BUILDING seq-gen ******************** "
	@echo "*********************************************** "	
	make seq-gen
	@echo "*********************************************** "
	@echo "********* COPY MISC FILES ********************* "
	@echo "*********************************************** "
	cp ../misc/* $(MISC_DIR)	
	
seq-gen:
	echo $(CWD)	
	tar -zxvf ../ext/Seq-Gen.v1.3.3.tgz
	make -C Seq-Gen.v1.3.3/source 
	mv $(CWD)/Seq-Gen.v1.3.3/source/seq-gen $(EXTERNAL_CODE)/seq-gen
	rm -rf Seq-Gen.v1.3.3
		
VirusTreeSimulator.jar:
	echo $(CWD)
	cp 	../ext/VirusTreeSimulator.jar $(EXTERNAL_CODE)/VirusTreeSimulator.jar
	
