\name{prog.HPTN071.input.parser.v1}
\alias{prog.HPTN071.input.parser.v1}
\title{Sequence sampler (version 1)}
\usage{
  prog.HPTN071.input.parser.v1()
}
\value{
  NULL. Saves simulations to file.
}
\description{
  Reads files from the epi simulator in directory
  \code{indir} and writes csv files in directory
  \code{outdir} for the virus tree simulator. The program
  samples sequences according to an exponentially
  increasing sampling fraction. All input arguments are
  specified via the \code{argv} string, see the Examples.
}
\examples{
##--------------------------------------------------------------------------------------------------------
##	simulate sequence sampling for epi simulation
##	sequences are sampled assuming an exponentially growing sequence sampling rate so that
##	1\% is sampled in 1980
##	25\% is sampled in 2020
##--------------------------------------------------------------------------------------------------------
indir			<- system.file(package="rPANGEAHIVsim", "misc")
indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
tmpdir.HPTN071	<- paste(tmpdir,'/HPTN071parser',sep='')
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'

#	output files
outfile.ind		<- '140716_RUN001_IND.csv'
outfile.trm		<- '140716_RUN001_TRM.csv'
dir.create(tmpdir.HPTN071, showWarnings=FALSE)

#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v1(indir, infile.trm, infile.ind, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))

#	run the parser
prog.HPTN071.input.parser.v1()
}

