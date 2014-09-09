\name{cmd.VirusTreeSimulator}
\alias{cmd.VirusTreeSimulator}
\title{Command line generator for the \code{VirusTreeSimulator}}
\usage{
  cmd.VirusTreeSimulator(indir, infile.trm, infile.ind,
    outdir, outfile, prog = PR.VIRUSTREESIMULATOR,
    prog.args = "-demoModel Logistic -N0 0.1 -growthRate 1.5 -t50 -4")
}
\description{
  The \code{VirusTreeSimulator} reads files from the
  sequence sampler in directory \code{indir} and writes
  detailed nexus files in directory \code{outdir} for the
  virus tree simulator. The program generates within-host
  phylogenies for sampled and unsampled individuals in a
  transmission chain along the specified within host
  coalescent model. Within-host phylogenies are then
  concatenated into a between-host phylogeny for each
  transmission chain.
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
##	Then call virus tree simulator with input args as below
##--------------------------------------------------------------------------------------------------------
indir			<- system.file(package="rPANGEAHIVsim", "misc")
indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
#	re-name the following:
tmpdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
#	simulation input files from the epi-simulator
infile.ind		<- '140716_RUN001_IND.csv'
infile.trm		<- '140716_RUN001_TRM.csv'

#	
#	step: sequence sampler
#
tmpdir.HPTN071	<- paste(tmpdir,'/HPTN071parser',sep='')
dir.create(tmpdir.HPTN071, showWarnings=FALSE)
#	output files from the sequence sampler
outfile.ind		<- '140716_RUN001_IND.csv'
outfile.trm		<- '140716_RUN001_TRM.csv'
#	get input into 'argv'. this is needed because the input parser is usually called from the command line, and 'argv' mimics the way input is provided when the parser is called from the command line
cmd				<- cmd.HPTN071.input.parser.v2(indir, infile.trm, infile.ind, tmpdir.HPTN071,  infile.trm, infile.ind)				 
argv			<<- unlist(strsplit(cmd,' '))
#	run the sequence sampler
prog.HPTN071.input.parser.v2()
#	
#	step: virus tree sampler
#
tmpdir.VTS		<- paste(tmpdir,'/VirusTreeSimulator',sep='')
dir.create(tmpdir.VTS, showWarnings=FALSE)
#	output file prefix for virus tree sampler
outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
prog.args		<- '-demoModel Logistic -N0 100000 -growthRate 0.0001 -t50 -0.04'
#	Ne=1.5e5 times 2 days generation time
#prog.args		<- '-demoModel Logistic -N0 300000 -growthRate XXX -t50 -XXX'
cmd				<- cmd.VirusTreeSimulator(tmpdir.HPTN071, infile.trm, infile.ind, tmpdir.VTS, outfile, prog.args=prog.args)	
cat(cmd)
}

