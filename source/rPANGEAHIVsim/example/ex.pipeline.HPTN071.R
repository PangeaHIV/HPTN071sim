##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given HPTN071 epi simulation  
##--------------------------------------------------------------------------------------------------------
indir			<- system.file(package="rPANGEAHIVsim", "misc")
indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140914'
dir.create(tmpdir, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.PREV.min=0.01, s.PREV.max=0.2, epi.dt=1/48, epi.import=0.1 )	
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)

