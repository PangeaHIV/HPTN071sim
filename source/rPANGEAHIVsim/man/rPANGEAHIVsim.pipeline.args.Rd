\name{rPANGEAHIVsim.pipeline.args}
\alias{rPANGEAHIVsim.pipeline.args}
\title{Arguments for rPANGEAHIV simulation pipeline}
\usage{
  rPANGEAHIVsim.pipeline.args(yr.start = 1980,
    yr.end = 2020, seed = 42, s.INC.recent = 0.1,
    s.INC.recent.len = 5, s.PREV.min = 0.01,
    s.PREV.max = 0.25, epi.dt = 1/48, epi.import = 0.1,
    startseq.backdate = 40)
}
\arguments{
  \item{seed}{Random number seed}

  \item{yr.start}{Start year of the epi simulation (default
  1980)}

  \item{yr.end}{First year after the epi simulation
  (default 2020)}

  \item{s.INC.recent}{Proportion of incident cases sampled
  (not used)}

  \item{s.INC.recent.len}{Number of last year in which an
  exact proportion of incident cases is sampled (not used)}

  \item{s.PREV.min}{Proportion of infected cases sampled at
  start of the simulation}

  \item{s.PREV.max}{Proportion of infected cases sampled at
  the end of the simulation}

  \item{epi.dt}{Time increment of the epi simulation
  (default 1/48)}

  \item{epi.import}{Proportion of imported cases of all
  cases (default 0.1)}

  \item{startseq.backdate}{Hack to backdate the time of
  index cases to draw a starting sequence (default 0)}
}
\value{
  data.table
}
\description{
  Construct table with all input arguments to the
  rPANGEAHIV simulation pipeline.
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
#	input arguments for the pipeline
pipeline.args	<- rPANGEAHIVsim.pipeline.args( yr.start=1980, yr.end=2020, seed=42, s.PREV.min=0.01, s.PREV.max=0.25, epi.dt=1/48, epi.import=0.1 )	
#	
#	call simulation pipeline
#	this generates a UNIX batch file if no HPC system is detected, or
#	this generates and runs a qsub file if an HPC system is detected 
#
file			<- rPANGEAHIVsim.pipeline(indir, infile.ind, infile.trm, tmpdir, pipeline.args=pipeline.args)
cat(file)

}

