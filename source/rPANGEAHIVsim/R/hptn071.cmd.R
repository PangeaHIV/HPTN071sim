PR.PACKAGE					<- "rPANGEAHIVsim"
PR.STARTME					<- system.file(package=PR.PACKAGE, "misc", "rPANGEAHIV.startme.R")
PR.HPTN071.INPUT.PARSER1	<- paste(PR.STARTME,"-exe=HPTN071.INPUT.PARSER1",sep=' ')
PR.HPTN071.INPUT.PARSER2	<- paste(PR.STARTME,"-exe=HPTN071.INPUT.PARSER2",sep=' ')
PR.SEQGEN.FILECREATOR		<- paste(PR.STARTME,"-exe=PR.SEQGEN.FILECREATOR",sep=' ')
PR.SEQGEN.READER			<- paste(PR.STARTME,"-exe=PR.SEQGEN.READER",sep=' ')
PR.VIRUSTREESIMULATOR		<- system.file(package=PR.PACKAGE, "ext", "VirusTreeSimulator.jar")
PR.SEQGEN					<- system.file(package=PR.PACKAGE, "ext", "seq-gen")

HPC.MPIRUN					<- {tmp<- c("mpirun","mpiexec"); names(tmp)<- c("debug","cx1.hpc.ic.ac.uk"); tmp}
HPC.CX1.IMPERIAL			<- "cx1.hpc.ic.ac.uk"		#this is set to system('domainname',intern=T) for the hpc cluster of choice
HPC.MEM						<- "1750mb"
HPC.CX1.IMPERIAL.LOAD		<- "module load intel-suite mpi R/2.15"


######################################################################################
cmd.hpcsys<- function()
{
	tmp<- system('domainname',intern=T)
	if(!nchar(tmp))	tmp<- "debug"
	tmp
}
##--------------------------------------------------------------------------------------------------------
##	batch file wrapper
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.hpcwrapper<- function(cmd, hpcsys= cmd.hpcsys(), hpc.walltime=24, hpc.mem="1750mb", hpc.nproc=1, hpc.q=NA)
{
	wrap<- "#!/bin/sh"
	#hpcsys<- HPC.CX1.IMPERIAL
	if(hpcsys%in%c(HPC.CX1.IMPERIAL,'(none)'))
	{				
		tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
		wrap<- paste(wrap, tmp, sep='\n')		
		tmp	<- paste("#PBS -l select=1:ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
		wrap<- paste(wrap, tmp, sep='\n')
		wrap<- paste(wrap, "#PBS -j oe", sep='\n')
		if(!is.na(hpc.q))
			wrap<- paste(wrap, paste("#PBS -q",hpc.q), sep='\n\n')
		wrap<- paste(wrap, HPC.LOAD, sep='\n')
	}
	else if(hpcsys=='debug')
		cat(paste("\ndetected no HPC system and no hpcwrapper generated, domain name is",hpcsys))
	else
		stop(paste("unknown hpc system with domain name",hpcsys))
	
	cmd<- lapply(seq_along(cmd),function(i){	paste(wrap,cmd[[i]],sep='\n')	})
	if(length(cmd)==1)
		cmd<- unlist(cmd)
	cmd	
}
##--------------------------------------------------------------------------------------------------------
##	batch file caller
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
#' @export
cmd.hpccaller<- function(outdir, outfile, cmd)
{
	if( nchar( Sys.which("qsub") ) )
	{
		file	<- paste(outdir,'/',outfile,'.qsub',sep='')
		cat(paste("\nwrite HPC script to",file,"\n"))
		cat(cmd,file=file)
		cmd		<- paste("qsub",file)
		cat( cmd )
		cat( system(cmd, intern=TRUE) )
		Sys.sleep(1)
	}
	else
	{
		file	<- paste(outdir,'/',outfile,'.sh',sep='')
		cat(paste("\nwrite Shell script to\n",file,"\nStart this shell file manually\n"))
		cat(cmd,file=file)
		Sys.chmod(file, mode = "777")	
		Sys.sleep(1)
	}
	
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v1'
##	olli originally written 19-08-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v1}
#' @example example/ex.seq.sampler.v1.R
#' @export
cmd.HPTN071.input.parser.v1<- function(indir, infile.trm, infile.ind, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER1 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v1
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v1
#######################################################\n",sep='')
	cmd
}
##--------------------------------------------------------------------------------------------------------
##	command line generator for 'prog.HPTN071.input.parser.v2'
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Command line generator for \code{HPTN071.input.parser.v2}
#' @example example/ex.seq.sampler.v2.R
#' @export
cmd.HPTN071.input.parser.v2<- function(indir, infile.trm, infile.ind, outdir, outfile.trm, outfile.ind, prog=PR.HPTN071.INPUT.PARSER2 )	
{
	cmd<- "#######################################################
# start: run HPTN071.input.parser.v2
			#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir=', indir,' -infile.trm=',infile.trm,' -infile.ind=',infile.ind,' -outdir=',outdir,' -outfile.ind=',outfile.ind,' -outfile.trm=',outfile.trm,' \n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run HPTN071.input.parser.v2
#######################################################\n",sep='')
	cmd
}

######################################################################################
#	return command line to run the VirusTreeSimulator	
#	olli originally written 19-08-2014
#	return 		character string 
#' @title Command line generator for the \code{VirusTreeSimulator}
#' @description The \code{VirusTreeSimulator} reads files from the sequence sampler in directory \code{indir} and writes detailed nexus files
#' in directory \code{outdir} for the virus tree simulator. The program generates within-host phylogenies
#' for sampled and unsampled individuals in a transmission chain along the specified within host coalescent
#' model. Within-host phylogenies are then concatenated into a between-host phylogeny.
#' @example example/ex.virus.tree.simulator.R
#' @export
cmd.VirusTreeSimulator<- function(indir, infile.trm, infile.ind, outdir, outfile, prog=PR.VIRUSTREESIMULATOR, prog.args='-demoModel Logistic -N0 0.1 -growthRate 1.5 -t50 -4')
{
	cmd<- "#######################################################
# start: run VirusTreeSimulator
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste('java -Xms64m -Xmx400m -jar ',prog,' ',prog.args,'  ', indir,'/',infile.trm,' ',indir,'/',infile.ind,' ',outdir,'/',outfile, '\n', sep=''))
	cmd		<- paste(cmd, 'find ',outdir,' -name "*simple*" -delete\n', sep='')
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run VirusTreeSimulator
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run SeqGen Input File Creator	
#	olli originally written 26-08-2014
#	return 		character string 
cmd.SeqGen.createInputFiles<- function(indir.epi, infile.epi, indir.vts, infile.vts, outdir, prog=PR.SEQGEN.FILECREATOR)
{
	cmd<- "#######################################################
# start: run SeqGen.createInputFile 
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir.epi=', indir.epi,' -infile.epi=',infile.epi,' -indir.vts=', indir.vts,' -infile.vts=',infile.vts,' -outdir=',outdir,'\n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen.createInputFile
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run SeqGen Output File Reader	
#	olli originally written 26-08-2014
#	return 		character string 
cmd.SeqGen.readOutputFiles<- function(indir.epi, infile.epi, indir.sg, outdir, prog=PR.SEQGEN.READER)
{
	cmd<- "#######################################################
# start: run SeqGen.readOutputFiles 
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' -indir.epi=', indir.epi,' -infile.epi=',infile.epi,' -indir.sg=', indir.sg,' -outdir=',outdir,'\n', sep=''))
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen.readOutputFiles
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run Seq-Gen-1.3.3	
#	olli originally written 26-08-2014
#	return 		character string 
cmd.SeqGen<- function(indir, infile, outdir, outfile, prog=PR.SEQGEN, prog.args='-n1 -k1 -on -z42', alpha=1, gamma=4, invariable=0, scale=1, 
						freq.A=0.25, freq.C=0.25, freq.G=0.25, freq.T=0.25,
						rate.AC=1, rate.AG=1, rate.AT=1, rate.CG=1, rate.CT=1, rate.GT=1)
{
	cmd<- "#######################################################
# start: run SeqGen
#######################################################"
	cmd		<- paste(cmd, paste("\necho \'run ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd, paste(prog,' ',prog.args,' ',sep=''),sep='')
	#	add substitution model
	#cmd		<- paste(cmd, paste("-mHKY -t3.0 -f0.3,0.2,0.2,0.3",sep=''),sep='')
	cmd		<- paste(cmd, paste('-mGTR -a',alpha,' -g',gamma,' -i',invariable, ' -s', scale,
									' -f',freq.A,',',freq.C,',',freq.G,',',freq.T,
									' -r',rate.AC, ',', rate.AG, ',', rate.AT, ',', rate.CG, ',', rate.CT, ',', rate.GT, sep=''),sep='')
	#	add I/O
	cmd		<- paste(cmd, paste(' < ', indir,'/',infile,' > ', outdir,'/',outfile, '\n', sep=''))
	cmd		<- paste(cmd, 'rm ',indir,'/',infile,'\n', sep='')
	cmd		<- paste(cmd,paste("echo \'end ",prog,"\'\n",sep=''))
	cmd		<- paste(cmd,"#######################################################
# end: run SeqGen
#######################################################\n",sep='')
	cmd
}
######################################################################################
#	return command line to run a batch of Seq-Gen-1.3.3 jobs for a PANGEA simulation	
#	olli originally written 26-08-2014
#	return 		character string 
cmd.PANGEA.SeqGen<- function(indir.sg, infile.prefix, plot.file=NA)
{	
	s.seed		<- 42
	set.seed(s.seed)
	file		<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#	check posterior samples
	tmp			<- log.df[, list(N=length(a)), by=c('GENE','CODON_POS')]
	stopifnot( tmp[, length(unique(N))==1] )
	n.samples	<- tmp[1,N] 
			
	#	create SeqGen input files
	df.ph.out	<- df.seqgen[,	{
				#print( table( strsplit(ANCSEQ, ''), useNA='if') )
				file	<- paste(indir.sg,'/',infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.seqgen',sep='')
				cat(paste('\nwrite to file',file))
				txt		<- paste('1\t',nchar(ANCSEQ),'\n',sep='')
				txt		<- paste(txt, 'ANCSEQ\t',toupper(ANCSEQ),'\n',sep='')					
				txt		<- paste(txt, '1\n',sep='')
				txt		<- paste(txt, NEWICK, '\n',sep='')
				cat(txt, file=file)
				list(FILE= paste(infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.seqgen',sep='') )
				# ./seq-gen -mHKY -t3.0 -f0.3,0.2,0.2,0.3 -n1 -k1 -on < /Users/Oliver/git/HPTN071sim/tmp/SeqGen/140716_RUN001_50.seqgen > example.nex
			}, by=c('GENE','CODON_POS','IDCLU')]
	df.ph.out	<- df.ph.out[, list(FILE=FILE, IDCLU=IDCLU, IDX=sample(n.samples, length(FILE), replace=FALSE)), by=c('GENE','CODON_POS')]
	if(1)
	{
		#	sample GTR parameters from posterior
		setkey(log.df, GENE, CODON_POS)
		log.df[, IDX:=seq_len(n.samples)]
		df.ph.out	<- merge(df.ph.out, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
		#	standardise mu for each FILE
		tmp			<- df.ph.out[, list(FILE=FILE, mu=mu/mean(mu)), by='IDCLU']
		df.ph.out	<- merge( subset(df.ph.out, select=which(colnames(df.ph.out)!='mu')), subset(tmp, select=c(FILE, mu)), by='FILE')
	}
	if(0)
	{
		#	pick central GTR parameters
		df.ph.out	<- merge(df.ph.out, gtr.central, by= c('GENE','CODON_POS'))
	}
	if(!is.na(plot.file))
	{
		tmp			<- subset(df.nodestat, select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.ph.out, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')		
		ggsave(file=plot.file, w=6, h=6)						
	}	
	#	create SeqGen command line
	df.ph.out	<- df.ph.out[, {												
					cmd	<- cmd.SeqGen(indir.sg, FILE, indir.sg, gsub('seqgen','phy',FILE), prog=PR.SEQGEN, prog.args=paste('-n',Global.gen.n.seq,' -k1 -or -z',Global.rng.seed,sep=''), 
							alpha=alpha, gamma=4, invariable=0, scale=mu, freq.A=a, freq.C=c, freq.G=g, freq.T=t,
							rate.AC=ac, rate.AG=ag, rate.AT=at, rate.CG=cg, rate.CT=1, rate.GT=gt)
					#cat(cmd)
					list(CMD=cmd)							
				}, by='FILE']
	#print(df.ph.out)
	cmd	<- paste(df.ph.out[, CMD], collapse='\n')
	cmd
}


