#' This file contains R functions that provide an interface to HPC scripting
#' Shell scripts are generated that can be either run directly, or submitted to an HPC system.

PR.PACKAGE					<- "rPANGEAHIVsim"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "rPANGEAHIV.startme.R")
PR.HPTN071.INPUT.PARSER		<- paste(PR.STARTME,"-exe=HPTN071.INPUT.PARSER",sep=' ')
PR.VIRUSTREESIMULATOR		<- system.file(package=PR.PACKAGE, "ext", "VirusTreeSimulator.jar")

HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi R/2.15"

##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser'
##	olli originally written 19-08-2014
##--------------------------------------------------------------------------------------------------------
cmd.HPTN071.input.parser<- function(indir, infile.trm, infile.ind, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,'\n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run the VirusTreeSimulator	
#	olli originally written 19-08-2014
#	return 		character string 
cmd.VirusTreeSimulator<- function(indir, infile.trm, infile.ind, outdir, outfile, prog=PR.VIRUSTREESIMULATOR, prog.args='-demoModel Logistic -N0 0.1 -growthRate 1.5 -t50 -4')
{
	cmd<- "#######################################################
# start: run VirusTreeSimulator
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste('java -Xms64m -Xmx400m -jar ',prog,' ',prog.args,'  ', indir,'/',infile.trm,' ',indir,'/',infile.ind,' ',outdir,'/',outfile, '\n', sep=''))
	cmd		<- paste(cmd, 'find . -name "*simple*" -delete\n')
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run VirusTreeSimulator
#######################################################\n",sep='')
	cmd
}