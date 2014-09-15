\name{rPANGEAHIVsim.pipeline}
\alias{rPANGEAHIVsim.pipeline}
\title{rPANGEAHIV simulation pipeline}
\usage{
  rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm,
    outdir, pipeline.args = rPANGEAHIVsim.pipeline.args())
}
\arguments{
  \item{indir}{Input directory}

  \item{infile.ind}{Input file with individual
  metavariables}

  \item{infile.trm}{Input file with transmission events}

  \item{pipeline.args}{Input arguments for the simulation}
}
\value{
  file name of qsub or UNIX batch file
}
\description{
  Reads two files \code{infile.ind} and \code{infile.trm}
  from the epi simulator in directory \code{indir} and
  produces a UNIX batch file that contains all the
  simulation steps in directory \code{outdir}.
}
\examples{
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

##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given DSPS epi simulation  
##--------------------------------------------------------------------------------------------------------
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140912'
dir.create(tmpdir, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.trm		<- '140911_DSPS_RUN002_TRM.csv'
infile.ind		<- NULL
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.PREV.min=0.01, s.PREV.max=0.70, epi.model='DSPS', epi.dt=1/48, epi.import=0.1 )	
infile.args		<- paste(tmpdir,'/',substr(infile.trm, 1, nchar(infile.trm)-7), 'PipeArgs.R',sep='')
save(pipeline.args, file=infile.args)
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)
}

