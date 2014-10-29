##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given HPTN071 epi simulation  
##--------------------------------------------------------------------------------------------------------
\dontrun{
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp141023'
dir.create(tmpdir, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '211014_RUN123_SCENARIO_0_IND.csv'
infile.trm		<- '211014_RUN123_SCENARIO_0_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42,
												s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.11, 
												epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
												v.N0tau=1, v.r=2.851904, v.T50=-2,
												wher.mu=log(0.00447743)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=0,
												dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='normal', startseq.mode='sample')						
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)
}
