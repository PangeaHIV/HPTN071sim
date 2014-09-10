\name{cmd.HPTN071.input.parser.v2}
\alias{cmd.HPTN071.input.parser.v2}
\title{Command line generator for \code{HPTN071.input.parser.v2}}
\usage{
  cmd.HPTN071.input.parser.v2(indir, infile.trm,
    infile.ind, infile.args, outdir, outfile.trm,
    outfile.ind, prog = PR.HPTN071.INPUT.PARSER2)
}
\value{
  command line string
}
\description{
  Command line generator for \code{HPTN071.input.parser.v2}
}
\examples{
##--------------------------------------------------------------------------------------------------------
##	simulate sequence sampling for epi simulation
##	sequences are sampled assuming an exponentially growing sequence sampling rate so that
##	1\% is sampled in 1980
##	25\% is sampled in 2020
##	In addition, 10\% of transmissions are broken and treated as imported from outside the simulated population.
##	The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
##	simulated population
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
cmd				<- cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))

#	run the parser
prog.HPTN071.input.parser.v2()
}

