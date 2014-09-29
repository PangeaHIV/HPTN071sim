\name{rPANGEAHIVsim.pipeline.args}
\alias{rPANGEAHIVsim.pipeline.args}
\title{Arguments for rPANGEAHIV simulation pipeline}
\usage{
  rPANGEAHIVsim.pipeline.args(yr.start = 1980,
    yr.end = 2020, seed = 42, s.INC.recent = 0.1,
    s.INC.recent.len = 5, s.PREV.min = 0.01,
    s.PREV.max = 0.25, epi.model = "HPTN071",
    epi.dt = 1/48, epi.import = 0.1, v.N0tau = 35800,
    v.r = 2, v.T50 = -1, wher.mu = 0.005, wher.sigma = 0.8,
    bwerm.mu = 1.5, bwerm.sigma = 0.12,
    startseq.backdate = NA)
}
\arguments{
  \item{yr.start}{Start year of the epi simulation (default
  1980)}

  \item{yr.end}{First year after the epi simulation
  (default 2020)}

  \item{seed}{Random number seed}

  \item{s.INC.recent}{Proportion of incident cases sampled
  (not used)}

  \item{s.INC.recent.len}{Number of last year in which an
  exact proportion of incident cases is sampled (not used)}

  \item{s.PREV.min}{Proportion of infected cases sampled at
  start of the simulation}

  \item{s.PREV.max}{Proportion of infected cases sampled at
  the end of the simulation}

  \item{epi.model}{The epi model used to create the epi
  simulation (default HPTN071, alternatively DSPS)}

  \item{epi.dt}{Time increment of the epi simulation
  (default 1/48)}

  \item{epi.import}{Proportion of imported cases of all
  cases (default 0.1)}

  \item{v.N0tau}{Parameter of BEAST::LogisticGrowthN0::N0
  (default 3.58e4)}

  \item{v.r}{Parameter of BEAST::LogisticGrowthN0::r
  (default 2)}

  \item{v.T50}{Parameter of BEAST::LogisticGrowthN0::T50}

  \item{wher.mu}{Mean within host evolutionary rate of log
  normal density}

  \item{wher.sigma}{Standard deviation in within host
  evolutionary rate of log normal density}

  \item{bwerm.mu}{Mean between host evolutionary rate
  multiplier of log normal density}

  \item{bwerm.sigma}{Standard deviation in between host
  evolutionary rate multiplier of log normal density}

  \item{startseq.backdate}{Hack to backdate the time of
  index cases to draw a starting sequence (default -1)}
}
\value{
  data.table
}
\description{
  Construct table with all input arguments to the
  rPANGEAHIV simulation pipeline. The demographic within
  host coealescent model is
  N0tau*(1+exp(-r*T50))/(1+exp(-r*(T50-t))). Default
  parameters are set so that the curve asymptotes at
  Ne*tau=3e5 and reaches half its asymptotic value 1 year
  post infection. Branch lengths are multiplied by within
  host evolutionary rates, and within host branch lengths
  ending in a transmission event are multiplied in addition
  with a between-host evolutionary rate multiplier.
}
\examples{
##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given HPTN071 epi simulation  
##--------------------------------------------------------------------------------------------------------
\\dontrun{
indir			<- system.file(package="rPANGEAHIVsim", "misc")
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140915'
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
}
##--------------------------------------------------------------------------------------------------------
##	example pipeline to simulate sequences for a given DSPS epi simulation  
##--------------------------------------------------------------------------------------------------------
\\dontrun{
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
}

