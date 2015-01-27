##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given HPTN071 epi simulation  
##--------------------------------------------------------------------------------------------------------
\dontrun{
indir			<- '/Users/Oliver/git/HPTN071sim/source/rPANGEAHIVsim/inst/misc'
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp150126'
dir.create(tmpdir, showWarnings=FALSE)
tmpdir.HPTN071	<- paste(tmpdir,'/TrChains',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '150127_RUN001_IND.csv'
infile.trm		<- '150127_RUN001_TRM.csv'
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.MODEL='Prop2DiagB4I',
						s.PREV.max=0.08, s.INTERVENTION.start=2015, s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50,
						epi.model='HPTN071', epi.dt=1/48, epi.import=0.05,
						v.N0tau=1, v.r=2.851904, v.T50=-2,
						wher.mu=log(0.003358613)-0.3^2/2, wher.sigma=0.3, bwerm.mu=log(0.002239075)-0.13^2/2, bwerm.sigma=0.13, er.gamma=NA,
						dbg.GTRparam=0, dbg.rER=0, index.starttime.mode='fix1955', startseq.mode='one', seqtime.mode='AtDiag')								
infile.args		<- paste(tmpdir,'/',substr(infile.ind, 1, nchar(infile.ind)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)						
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)
}
