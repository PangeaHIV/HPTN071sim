\name{rPANGEAHIVsim.pipeline}
\alias{rPANGEAHIVsim.pipeline}
\title{rPANGEAHIV simulation pipeline}
\usage{
  rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm,
    outdir)
}
\arguments{
  \item{indir}{Input directory}

  \item{infile.ind}{Input file with individual
  metavariables}

  \item{infile.trm}{Input file with transmission events}
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
##	example pipeline to simulate sequences for a given epi simulation  
##--------------------------------------------------------------------------------------------------------
indir			<- system.file(package="rPANGEAHIVsim", "misc")
indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140909'
dir.create(tmpdir, showWarnings=FALSE)
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir)
cat(file)

}

