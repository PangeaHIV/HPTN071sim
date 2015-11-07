######################################################################################
#	Program to simulate sequences with Seq-Gen-1.3.2 via phyclust 	
#	olli originally written 15-09-2014
######################################################################################
prog.PANGEA.SeqGen.run.WINDOWScompatible<- function()
{	
	indir.epi			<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi			<- '140716_RUN001_SAVE.R'		
	indir.sg			<- '/Users/Oliver/git/HPTN071sim/tmp140908/SeqGen'
	infile.prefix		<- '140716_RUN001_'
	infile.args			<- NA
	outdir				<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	with.plot			<- 1	
	verbose				<- 1
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'	
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.sg= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.sg<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infile.sg= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]				
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.sg, infile.prefix, sep='\n'))
	}
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed')%in%pipeline.args[, stat] ) )	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	file		<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nLoad seqgen R input file, file=',file))
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#
	log.df[, IDX:= seq_len(nrow(log.df))]
	log.df[, FILE:=NULL]	 
	#	sample GTR parameters from posterior
	tmp			<- df.seqgen[, {
									tmp	<- log.df[['IDX']][ which( log.df[['GENE']]==GENE & log.df[['CODON_POS']]==CODON_POS ) ]
									stopifnot(length(tmp)>1)
									list( IDCLU=IDCLU, IDX=sample(tmp, length(IDCLU), replace=FALSE) )	
								}, by=c('GENE','CODON_POS')]
	tmp			<- merge(tmp, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
	df.seqgen	<- merge(df.seqgen, tmp, by=c('GENE','CODON_POS','IDCLU'))
	#
	if(with.plot)
	{
		tmp			<- subset(df.nodestat, ER!=log.df$meanRate[1], select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.seqgen, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')
		file		<- paste(indir.sg,'/',infile.prefix, 'ER_by_gene.pdf',sep='')
		ggsave(file=file, w=6, h=6)						
	}			
	#
	#	run SeqGen-1.3.2 (no seed) 
	#
	df.ph.out	<- df.seqgen[,	{									
									file	<- paste(indir.sg,'/',infile.prefix, IDCLU, '_', GENE, '_', CODON_POS,'.phy',sep='')
									cat(paste('\nwrite to file',file))
									opts	<- paste('-n1 -k1 -or -mGTR -a',alpha,' -g4 -i0 -s',mu,' -f',a,',',c,',',g,',',t,' -r',ac,',',ag,',',at,',',cg,',',1,',',gt,sep='')
									input 	<- c(paste(" 1 ", nchar(ANCSEQ), sep=''), paste('ANCSEQ\t', toupper(ANCSEQ), sep = ''), 1, NEWICK)
									z		<- seqgen(opts, input=input, temp.file=file)				
								}, by=c('GENE','CODON_POS','IDCLU')]
	#
	#	process SeqGen runs
	#
	#	collect simulated sequences
	infiles		<- list.files(indir.sg)
	infiles		<- infiles[ grepl('*phy$', infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')		
	infile.df	<- data.table(FILE=infiles)
	tmp			<- infile.df[, strsplit(FILE, '_') ]
	infile.df[, CODON_POS:= sapply(tmp, function(x) rev(x)[label.idx.codonpos])]
	infile.df[, GENE:= sapply(tmp, function(x) rev(x)[label.idx.gene])]
	infile.df[, IDCLU:= sapply(tmp, function(x) rev(x)[label.idx.clu])]
	set(infile.df, NULL, 'CODON_POS', infile.df[, substr(CODON_POS,1,3)])
	cat(paste('\nFound sequences for clusters, nclu=', infile.df[, length(unique(IDCLU))]))
	#
	#	read simulated sequences
	#
	df.seq		<- infile.df[,	{
									cat(paste('\nread seq in file',FILE))									
									file	<- paste(indir.sg,'/',FILE,sep='')
									tmp		<- as.character(read.dna(file, format='sequential'))
									tmp		<- tmp[!grepl('NOEXIST',rownames(tmp)), , drop=FALSE]
									list( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), LABEL=rownames(tmp) )				
								}, by='FILE']
	df.seq		<- merge(df.seq, infile.df, by='FILE')
	#
	#	reconstruct genes from codon positions
	#
	df.seq[, STAT:=paste(GENE,CODON_POS,sep='.')]		
	df.seq		<- dcast.data.table(df.seq, IDCLU + LABEL ~ STAT, value.var="SEQ")
	#	check that seq of correct size
	stopifnot( df.seq[, length(unique(LABEL))]==nrow(df.seq) )
	stopifnot( df.seq[, all( nchar(GAG.CP1)==nchar(GAG.CP2) & nchar(GAG.CP1)==nchar(GAG.CP3) )] )
	stopifnot( df.seq[, all( nchar(POL.CP1)==nchar(POL.CP2) & nchar(POL.CP1)==nchar(POL.CP3) )] )
	stopifnot( df.seq[, all( nchar(ENV.CP1)==nchar(ENV.CP2) & nchar(ENV.CP1)==nchar(ENV.CP3) )] )
	#
	df.seq		<- df.seq[, {
								tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
								env		<- paste(as.vector(tmp), collapse='')
								tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
								gag		<- paste(as.vector(tmp), collapse='')
								tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
								pol		<- paste(as.vector(tmp), collapse='')
								list(GAG=gag, POL=pol, ENV=env, IDCLU=IDCLU)
							}, by=c('LABEL')]
	#	check that we have indeed seq for all sampled individuals
	df.seq		<- subset( df.seq, !grepl('NOEXIST',LABEL) )	
	tmp			<- df.seq[, sapply( strsplit(LABEL, treelabel.idx.sep, fixed=TRUE), '[[', treelabel.idx.idpop )]
	df.seq[, IDPOP:=as.integer(substr(tmp,7,nchar(tmp)))]	
	stopifnot( setequal( subset( df.inds, !is.na(TIME_SEQ) )[, IDPOP], df.seq[,IDPOP]) )
	#	merge simulated data
	simu.df		<- merge(subset(df.seq, select=c(IDPOP, GAG, POL, ENV)), subset( df.inds, !is.na(TIME_SEQ) ), by='IDPOP', all.x=TRUE)
	#
	#	save simulated data -- internal
	#
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, CIRCM, DOB, DOD, TIME_SEQ ) )
	if(pipeline.args['epi.model'][,v]=='DSPS')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, TIME_SEQ ) )
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(tmp, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	tmp				<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(tmp, file, format = "fasta")
	#
	#	zip simulated files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	#
	#	zip internal files
	#
	tmp				<- list.files(outdir, pattern='*pdf$', recursive=TRUE, full.names=TRUE) 
	file.copy(tmp, outdir, overwrite=TRUE)	 
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep=''), list.files(outdir, pattern='*pdf$', recursive=FALSE, full.names=TRUE) ) 	
	zip( paste(outdir, '/', infile.prefix, 'INTERNAL.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	
	return(1)
}
######################################################################################
#	Program to add gaps into sequences  	
#	olli originally written 23-06-2015
######################################################################################
prog.PANGEA.AddGaps.run.v1<- function()
{	
	indir.simu		<- '/Users/Oliver/git/HPTN071sim/treec150623/nogaps'
	indir.gap		<-	'~/git/HPTN071sim/treec150623/PANGEAreal'
	infile.simu		<- '150227_HPTN071_TRAIN1_SIMULATED'
	infile.gap		<- '150623_GlobalAlignment_cov1.fasta'
	outdir			<- '/Users/Oliver/git/HPTN071sim/treec150623/withgaps'
	outfile.cov		<- regmatches(infile.gap,regexpr('cov[0-9]+',basename(infile.gap)))	
	gap.country		<- 'ZA'
	gap.symbol		<- '?'
	gap.seed		<- 42
	outfile			<- paste(infile.simu, '_', gap.country, outfile.cov, '.fa', sep='')
	verbose			<- 1
	#
	#	read args
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.simu= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.simu<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.gap= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.gap<- tmp[1]	
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.gap= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.gap<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.simu= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.simu<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,8),
									outfile= return(substr(arg,10,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile<- tmp[1]
		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									gap.country= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) gap.country<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									gap.symbol= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) gap.symbol<- tmp[1]					
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.simu, indir.gap, infile.simu, infile.gap, outdir, outfile.cov, outfile, gap.country, gap.symbol, gap.seed	, sep='\n'))
	}
	sgp		<- PANGEA.add.gaps(indir.simu, indir.gap, infile.simu, infile.gap, gap.country, gap.symbol, gap.seed, verbose=1)
	write.dna(sgp, file=paste(outdir, outfile, sep='/'), format='fasta', colsep='', nbcol=-1)	
}
######################################################################################
#	Program to simulate sequences with Seq-Gen-1.3.3 	
#	olli originally written 09-09-2014
#	used up to version prog.HPTN071.input.parser.v3
######################################################################################
#' @title Program to simulate gene sequences
#' @description \code{prog.PANGEA.SeqGen.run} reads file \code{infile.sg} in directory \code{indir.sg} that was
#' created with the \code{SeqGen} input file creator. The simulated partial sequences are collected, coerced back
#' into Gag, Pol, Env genes, and written in fasta format to directory \code{outdir}. Patient Metavariables are 
#' stored in the same directory, and zip files are created.
#' @return NULL. Saves zip files with simulations.
#' @example example/ex.seqgen.run.R
#' @export
prog.PANGEA.SeqGen.run<- function()
{	
	indir.epi			<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi			<- '140716_RUN001_SAVE.R'		
	indir.sg			<- '/Users/Oliver/git/HPTN071sim/tmp140908/SeqGen'
	infile.prefix		<- '140716_RUN001_'
	infile.args			<- NA
	outdir				<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	with.plot			<- 1	
	verbose				<- 1
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'	
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.sg= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.sg<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infile.sg= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]				
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.sg, infile.prefix, sep='\n'))
	}
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	file			<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nLoad seqgen R input file, file=',file))
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#
	if( pipeline.args['index.starttime.mode',][,v]=='shift' )		
		root.edge.rate	<- 1e-6
	if( pipeline.args['index.starttime.mode',][,v]!='shift' )
		root.edge.rate	<- log.df[1,meanRate]		
	stopifnot( all( c('s.seed')%in%pipeline.args[, stat] ) )	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])	
	#
	#	create SeqGen input files
	#
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
	df.ph.out	<- df.ph.out[, 	{
				tmp	<- log.df[['IDX']][ which( log.df[['GENE']]==GENE & log.df[['CODON_POS']]==CODON_POS ) ]
				stopifnot(length(tmp)>1)
				list( FILE=FILE, IDCLU=IDCLU, IDX=sample(tmp, length(FILE), replace=FALSE) )
			}, by=c('GENE','CODON_POS')]
	#	sample GTR parameters from posterior
	df.ph.out	<- merge(df.ph.out, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
	#
	if(with.plot)
	{
		tmp			<- subset(df.nodestat, ER!=root.edge.rate, select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.ph.out, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')
		file		<- paste(indir.sg,'/',infile.prefix, 'ER_by_gene.pdf',sep='')
		ggsave(file=file, w=6, h=6)						
	}	
	#
	#	call SeqGen command line
	#
	cat(paste('\nUsing Gamma rate variation, gamma=',pipeline.args['er.gamma',][, as.numeric(v)]))
	tmp		<- df.ph.out[, {	
				cat(paste('\nProcess', IDCLU, GENE, CODON_POS))				
				cmd	<- cmd.SeqGen(indir.sg, FILE, indir.sg, gsub('seqgen','phy',FILE), prog=PR.SEQGEN, prog.args=paste('-n',1,' -k1 -or -z',pipeline.args['s.seed',][, as.numeric(v)],sep=''), 
						alpha=alpha, gamma=pipeline.args['er.gamma',][, as.numeric(v)], invariable=0, scale=mu, freq.A=a, freq.C=c, freq.G=g, freq.T=t,
						rate.AC=ac, rate.AG=ag, rate.AT=at, rate.CG=cg, rate.CT=1, rate.GT=gt)
				system(cmd)				
				list(CMD=cmd)							
			}, by='FILE']
	#
	#	process SeqGen runs
	#
	#	collect simulated sequences
	infiles		<- list.files(indir.sg)
	infiles		<- infiles[ grepl('*phy$', infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')		
	infile.df	<- data.table(FILE=infiles)
	tmp			<- infile.df[, strsplit(FILE, '_') ]
	infile.df[, CODON_POS:= sapply(tmp, function(x) rev(x)[label.idx.codonpos])]
	infile.df[, GENE:= sapply(tmp, function(x) rev(x)[label.idx.gene])]
	infile.df[, IDCLU:= sapply(tmp, function(x) rev(x)[label.idx.clu])]
	set(infile.df, NULL, 'CODON_POS', infile.df[, substr(CODON_POS,1,3)])
	cat(paste('\nFound sequences for clusters, nclu=', infile.df[, length(unique(IDCLU))]))
	#
	#	read simulated sequences
	#
	df.seq		<- infile.df[,	{
				cat(paste('\nread seq in file',FILE))
				file	<- paste(indir.sg,'/',FILE,sep='')
				tmp		<- as.character(read.dna(file, format='sequential'))
				list( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), LABEL=rownames(tmp) )				
			}, by='FILE']
	df.seq		<- merge(df.seq, infile.df, by='FILE')
	#
	#	reconstruct genes from codon positions
	#
	df.seq[, STAT:=paste(GENE,CODON_POS,sep='.')]		
	df.seq		<- dcast.data.table(df.seq, IDCLU + LABEL ~ STAT, value.var="SEQ")
	#	check that seq of correct size
	stopifnot( df.seq[, all( nchar(GAG.CP1)==nchar(GAG.CP2) & nchar(GAG.CP1)==nchar(GAG.CP3) )] )
	stopifnot( df.seq[, all( nchar(POL.CP1)==nchar(POL.CP2) & nchar(POL.CP1)==nchar(POL.CP3) )] )
	stopifnot( df.seq[, all( nchar(ENV.CP1)==nchar(ENV.CP2) & nchar(ENV.CP1)==nchar(ENV.CP3) )] )
	#
	df.seq		<- df.seq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, IDCLU=IDCLU)
			}, by=c('LABEL')]
	#	check that we have indeed seq for all sampled individuals
	df.seq		<- subset( df.seq, !grepl('NOEXIST',LABEL) )	
	tmp			<- df.seq[, sapply( strsplit(LABEL, treelabel.idx.sep, fixed=TRUE), '[[', treelabel.idx.idpop )]
	df.seq[, IDPOP:=as.integer(substr(tmp,7,nchar(tmp)))]	
	stopifnot( setequal( subset( df.inds, !is.na(TIME_SEQ) )[, IDPOP], df.seq[,IDPOP]) )
	#
	#	save simulated data -- internal
	#
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
	{
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, CIRCM, DOB, DOD, TIME_SEQ, CD4_SEQ, INCIDENT_SEQ ) )
		setnames(tmp, 'INCIDENT_SEQ', 'INCIDENT_WITHIN1YEAR_SEQ')
	}		
	if(pipeline.args['epi.model'][,v]=='DSPS')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, TIME_SEQ ) )
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.gag		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(df.seq.gag, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.pol		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(df.seq.pol, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.env		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(df.seq.env, file, format = "fasta")
	if(with.plot)
	{
		#
		#	create and plot NJ tree on conc seq
		#			
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		#	concatenate sequences
		tmp				<- tolower(do.call('rbind',strsplit(df.seq[, paste(GAG,POL,ENV,sep='')],'')))
		rownames(tmp)	<- df.seq[, paste(IDCLU,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		#	plot
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJconc.pdf', sep='')				
		pdf(file=file, w=10, h=80)
		plot(seq.ph, show.tip=TRUE, cex=0.5)
		dev.off()		
	}
	#
	#	zip simulated files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	#
	#	zip internal files
	#
	tmp				<- list.files(outdir, pattern='*pdf$', recursive=TRUE, full.names=TRUE) 
	file.copy(tmp, outdir, overwrite=TRUE)	 
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep=''), list.files(outdir, pattern='*pdf$', recursive=FALSE, full.names=TRUE) ) 	
	zip( paste(outdir, '/', infile.prefix, 'INTERNAL.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	
	return(1)
}
######################################################################################
#	Program to simulate sequences with Seq-Gen-1.3.3 	
#	olli originally written 27-01-2015
######################################################################################
prog.PANGEA.SeqGen.run.v4<- function()
{	
	indir.epi			<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi			<- '140716_RUN001_SAVE.R'		
	indir.sg			<- '/Users/Oliver/git/HPTN071sim/tmp140908/SeqGen'
	infile.prefix		<- '140716_RUN001_'
	infile.args			<- NA
	outdir				<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	with.plot			<- 1
	with.NJ				<- 1
	verbose				<- 1
	label.idx.codonpos	<- 1
	label.idx.gene		<- 2
	label.idx.clu		<- 3
	treelabel.idx.idpop	<- 1
	treelabel.idx.sep	<- '|'	
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,9),
									indir.sg= return(substr(arg,11,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.sg<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									infile.sg= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]				
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]				
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.sg, infile.prefix, sep='\n'))
	}
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	set.seed(pipeline.args['s.seed',][, as.numeric(v)])
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"   "df.sp"
	#
	file			<- paste(indir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nLoad seqgen R input file, file=',file))
	load(file)	#expect df.seqgen, gtr.central, log.df, df.nodestat
	#
	if( pipeline.args['index.starttime.mode',][,v]=='shift' )		
		root.edge.rate	<- 1e-6
	if( pipeline.args['index.starttime.mode',][,v]!='shift' )
		root.edge.rate	<- log.df[1,meanRate]		
	stopifnot( all( c('s.seed')%in%pipeline.args[, stat] ) )			
	#
	#	create SeqGen input files
	#
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
	df.ph.out	<- df.ph.out[, 	{
									tmp	<- log.df[['IDX']][ which( log.df[['GENE']]==GENE & log.df[['CODON_POS']]==CODON_POS ) ]
									stopifnot(length(tmp)>1)
									list( FILE=FILE, IDCLU=IDCLU, IDX=sample(tmp, length(FILE), replace=FALSE) )
								}, by=c('GENE','CODON_POS')]
	#	sample GTR parameters from posterior
	df.ph.out	<- merge(df.ph.out, log.df, by=c('GENE', 'CODON_POS', 'IDX'))
	stopifnot( nrow(df.ph.out)<65535 )	# running out of seeds?
	df.ph.out[, SEED:=sample.int(65535, nrow(df.ph.out))]
	#
	if(with.plot)
	{
		tmp			<- subset(df.nodestat, ER!=root.edge.rate, select=c(IDPOP, ER, BWM, IDCLU))
		tmp			<- merge(tmp, subset(df.ph.out, select=c(GENE, CODON_POS, IDCLU, mu)), by='IDCLU', allow.cartesian=TRUE)
		set(tmp, NULL, 'ER', tmp[, ER*mu])		
		
		ggplot(tmp, aes(x=CODON_POS, y=ER, colour=CODON_POS, group=CODON_POS)) + geom_boxplot() +				
				facet_grid(.~GENE, scales='free_y') +
				scale_colour_discrete(guide=FALSE) +
				scale_y_continuous(breaks= seq(0, 0.05, 0.002)) + labs(linetype='Gene', y='simulated within-host evolutionary rate', x='codon position')
		file		<- paste(indir.sg,'/',infile.prefix, 'ER_by_gene.pdf',sep='')
		ggsave(file=file, w=6, h=6)						
	}	
	#
	#	call SeqGen command line
	#
	cat(paste('\nUsing Gamma rate variation, gamma=',pipeline.args['er.gamma',][, as.numeric(v)]))
	tmp		<- df.ph.out[, {	
				cat(paste('\nProcess', IDCLU, GENE, CODON_POS))				
				cmd	<- cmd.SeqGen(indir.sg, FILE, indir.sg, gsub('seqgen','phy',FILE), prog=PR.SEQGEN, prog.args=paste('-n',1,' -k1 -or -z',SEED,sep=''), 
						alpha=alpha, gamma=pipeline.args['er.gamma',][, as.numeric(v)], invariable=0, scale=mu, freq.A=a, freq.C=c, freq.G=g, freq.T=t,
						rate.AC=ac, rate.AG=ag, rate.AT=at, rate.CG=cg, rate.CT=1, rate.GT=gt)
				system(cmd)				
				list(CMD=cmd)							
			}, by='FILE']
	#
	#	process SeqGen runs
	#
	#	collect simulated sequences
	infiles		<- list.files(indir.sg)
	infiles		<- infiles[ grepl('*phy$', infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')		
	infile.df	<- data.table(FILE=infiles)
	tmp			<- infile.df[, strsplit(FILE, '_') ]
	infile.df[, CODON_POS:= sapply(tmp, function(x) rev(x)[label.idx.codonpos])]
	infile.df[, GENE:= sapply(tmp, function(x) rev(x)[label.idx.gene])]
	infile.df[, IDCLU:= sapply(tmp, function(x) rev(x)[label.idx.clu])]
	set(infile.df, NULL, 'CODON_POS', infile.df[, substr(CODON_POS,1,3)])
	cat(paste('\nFound sequences for clusters, nclu=', infile.df[, length(unique(IDCLU))]))
	#
	#	read simulated sequences
	#
	df.seq		<- infile.df[,	{
				cat(paste('\nread seq in file',FILE))
				file	<- paste(indir.sg,'/',FILE,sep='')
				tmp		<- as.character(read.dna(file, format='sequential'))
				list( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), LABEL=rownames(tmp) )				
			}, by='FILE']
	df.seq		<- merge(df.seq, infile.df, by='FILE')
	#
	#	reconstruct genes from codon positions
	#
	df.seq[, STAT:=paste(GENE,CODON_POS,sep='.')]		
	df.seq		<- dcast.data.table(df.seq, IDCLU + LABEL ~ STAT, value.var="SEQ")
	#	check that seq of correct size
	stopifnot( df.seq[, all( nchar(GAG.CP1)==nchar(GAG.CP2) & nchar(GAG.CP1)==nchar(GAG.CP3) )] )
	stopifnot( df.seq[, all( nchar(POL.CP1)==nchar(POL.CP2) & nchar(POL.CP1)==nchar(POL.CP3) )] )
	stopifnot( df.seq[, all( nchar(ENV.CP1)==nchar(ENV.CP2) & nchar(ENV.CP1)==nchar(ENV.CP3) )] )
	#
	df.seq		<- df.seq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, IDCLU=IDCLU)
			}, by=c('LABEL')]
	#	check that we have indeed seq for all sampled individuals
	df.seq		<- subset( df.seq, !grepl('NOEXIST',LABEL) )	
	tmp			<- df.seq[, sapply( strsplit(LABEL, treelabel.idx.sep, fixed=TRUE), '[[', treelabel.idx.idpop )]
	df.seq[, IDPOP:=as.integer(substr(tmp,7,nchar(tmp)))]	
	stopifnot( setequal( subset( df.inds, !is.na(TIME_SEQ) )[, IDPOP], df.seq[,IDPOP]) )
	#
	#	save simulated data -- internal
	#
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep='')
	cat(paste('\nwrite to file', file))
	save(df.epi, df.trms, df.inds, df.sample, df.seq, df.sp, file=file)
	#
	#	save simulated data -- to be shared
	#	
	if(pipeline.args['epi.model'][,v]=='HPTN071')
	{
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, DOB, DOD, DIAG_T, DIAG_CD4, ART1_T, ART1_CD4, TIME_SEQ, RECENT_TR ) )
		tmp2		<- tmp[, which(!is.na(RECENT_TR))]
		cat(paste('\nSet RECENT_TR to missing for p=',1-pipeline.args['report.prop.recent',][,as.numeric(v)]))
		tmp2		<- sample(tmp2, round(length(tmp2)*(1-pipeline.args['report.prop.recent',][,as.numeric(v)])))
		set(tmp, NULL, 'RECENT_TR', tmp[, as.character(RECENT_TR)])
		set(tmp, tmp2, 'RECENT_TR', NA_character_)
		set(tmp, NULL, 'RECENT_TR', tmp[, factor(RECENT_TR)])	
		
		set(tmp, NULL, 'GENDER', tmp[,as.character(GENDER)])
		tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ<2000)]
		cat(paste('\nSet patient variables to NA for archival seq, n=',length(tmp2)))		
		set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
		set(tmp, tmp2, 'GENDER', NA_character_)
		tmp2		<- tmp[, which(is.na(DIAG_T) & TIME_SEQ>=2000)]
		cat(paste('\nSet patient variables to NA after 200, n=',length(tmp2)))
		print(tmp[tmp2,])
		set(tmp, tmp2, c('DOB','DOD'), NA_real_) 
		set(tmp, tmp2, 'GENDER', NA_character_)
		set(tmp, NULL, 'GENDER', tmp[,factor(GENDER)])
	}		
	if(pipeline.args['epi.model'][,v]=='DSPS')
		tmp			<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, GENDER, TIME_SEQ ) )
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep='')
	cat(paste('\nwrite to file', file))
	write.csv(tmp, file)		
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, GAG],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.gag		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='')
	write.dna(df.seq.gag, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, POL],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.pol		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep='')
	write.dna(df.seq.pol, file, format = "fasta")	
	tmp				<- tolower(do.call('rbind',strsplit(df.seq[, ENV],'')))
	rownames(tmp)	<- df.seq[, LABEL]
	df.seq.env		<- as.DNAbin(tmp)
	file			<- paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep='')
	write.dna(df.seq.env, file, format = "fasta")
	
	df.seq			<- merge(df.seq, df.seq[, list(IDCLU_N=length(IDPOP)), by='IDCLU'], by='IDCLU')
	if(with.NJ)
	{
		#
		#	create and plot NJ tree on conc seq
		#			
		#	load outgroup sequences
		file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
		cat(paste('\nLoading outgroup seq from file', file))
		load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
		if(nrow(df.seq)<2000)
			df.seq.nj	<- copy(df.seq)
		if(nrow(df.seq)>=2000)
		{
			cat(paste('\nToo many seqs for quick NJ tree calculation, selecting first 2000'))
			df.seq.nj	<- unique(subset(df.seq, select=c(IDCLU, IDCLU_N)))			
			df.seq.nj[, IDCLU_CN:=df.seq.nj[, cumsum(IDCLU_N)]]
			df.seq.nj	<- df.seq.nj[seq_len(max(1, which(IDCLU_CN>2000)[1]-1)), ]
			df.seq.nj	<- merge(subset(df.seq.nj, select=IDCLU), df.seq, by='IDCLU')
		}
		#	concatenate sequences	
		df.seq.nj[, TIPCOLOR:='black']
		set(df.seq.nj, df.seq.nj[,which(IDCLU_N<3)],'TIPCOLOR','red')
		tmp				<- tolower(do.call('rbind',strsplit(df.seq.nj[, paste(GAG,POL,ENV,sep='')],'')))		
		rownames(tmp)	<- df.seq.nj[, paste(IDCLU,'_',IDCLU_N,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)
		tmp				<- cbind(outgroup.seq.gag[,1:ncol(df.seq.gag)], outgroup.seq.pol, outgroup.seq.env)
		seq				<- rbind(seq,tmp)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJconc.pdf', sep='')	
		cat(paste('\nwrite to file',file))
		pdf(file=file, w=10, h=0.1*Ntip(seq.ph))
		plot(seq.ph, show.tip=TRUE, cex=0.5, tip.color=df.seq.nj[,TIPCOLOR])
		dev.off()	
		#
		#	create and plot NJ tree on pol seq
		#			
		tmp				<- tolower(do.call('rbind',strsplit(df.seq.nj[, POL],'')))		
		rownames(tmp)	<- df.seq.nj[, paste(IDCLU,'_',IDCLU_N,treelabel.idx.sep,LABEL,sep='')]	
		seq				<- as.DNAbin(tmp)		
		seq				<- rbind(seq,outgroup.seq.pol)	
		seq.ph			<- nj(dist.dna(seq, model='raw'))		
		tmp				<- which(seq.ph$tip.label=="HXB2")
		seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
		seq.ph			<- ladderize(seq.ph)
		file			<- paste(indir.sg, '/', infile.prefix, 'INFO_NJpol.pdf', sep='')	
		cat(paste('\nwrite to file',file))
		pdf(file=file, w=10, h=0.1*Ntip(seq.ph))
		plot(seq.ph, show.tip=TRUE, cex=0.5, tip.color=df.seq.nj[,TIPCOLOR])
		dev.off()			
	}
	#
	#	zip simulated sequence files
	#
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(indir.epi, '/', gsub('SAVE.R','CROSS_SECTIONAL_SURVEY.csv',infile.epi), sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_env.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_pol.fa', sep=''), paste(outdir, '/', infile.prefix, 'SIMULATED_gag.fa', sep='') )
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED_SEQ.zip', sep=''), tmp, flags = "-FSr9XTj")
	#
	#	zip simulated tree files
	#
	tmp2			<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_metadata.csv', sep=''), paste(indir.epi, '/', gsub('SAVE.R','CROSS_SECTIONAL_SURVEY.csv',infile.epi), sep=''), 
							paste(indir.epi, '/', infile.prefix, 'DATEDTREE.newick', sep=''),
							paste(indir.epi, '/', infile.prefix, 'DATEDCLUTREES.newick', sep=''),
							paste(indir.epi, '/', infile.prefix, 'SUBSTTREE.newick', sep=''))
	zip( paste(outdir, '/', infile.prefix, 'SIMULATED_TREE.zip', sep=''), tmp2, flags = "-FSr9XTj")
	#
	tmp				<- unique(c(tmp, tmp2))
	dummy			<- file.remove(tmp)	
	#
	#	zip internal files
	#
	tmp				<- list.files(outdir, pattern='*pdf$', recursive=TRUE, full.names=TRUE) 
	file.copy(tmp, outdir, overwrite=TRUE)	 
	tmp				<- c( paste(outdir, '/', infile.prefix, 'SIMULATED_INTERNAL.R', sep=''), list.files(outdir, pattern='*pdf$', recursive=FALSE, full.names=TRUE) ) 	
	zip( paste(outdir, '/', infile.prefix, 'INTERNAL.zip', sep=''), tmp, flags = "-FSr9XTj")
	dummy			<- file.remove(tmp)
	
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 17-08-2014
##--------------------------------------------------------------------------------------------------------
prog.PANGEA.SeqGen.createInputFile.v2<- function()
{
	stop()
	verbose			<- 1
	with.plot		<- 1
	label.sep		<- '|'	
	#
	#	read I/O
	#
	indir.epi		<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.vts		<- '/Users/Oliver/git/HPTN071sim/tmp140914/VirusTreeSimulator'
	infile.prefix	<- '140716_RUN001_'	
	infile.args		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_PipeArgs.R'
	outdir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp140914/SeqGen'	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.vts= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.vts<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.vts= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir.sg<- tmp[1]		
		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu')%in%pipeline.args[, stat] ) )
	#
	#	setup samplers
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol			<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	tmp				<- PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler(er.shift=pipeline.args['bwerm.mu',][, as.numeric(v)])
	rERbw			<- tmp$rERbw
	rERbw.args		<- tmp$rERbw.args
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp				<- PANGEA.RootSeq.create.sampler(root.ctime.grace= 0.5, sample.grace= 3)
	rANCSEQ			<- tmp$rANCSEQ
	rANCSEQ.args	<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df			<- PANGEA.GTR.params( )
	if( pipeline.args['dbg.rER',][, as.numeric(v)] )
	{
		cat(paste('\nSetting mus to mean per gene and codon_pos'))
		tmp		<- log.df[, list(mu= mean(mu)), by=c('GENE','CODON_POS')]
		#tmp[, ER:=mu*log.df$meanRate[1]]
		log.df	<- merge( subset(log.df, select=which(colnames(log.df)!='mu')), tmp, by=c('GENE','CODON_POS'))		
	}	
	#
	#
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#	
	infiles		<- list.files(indir.vts)
	tmp			<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles		<- infiles[ grepl(tmp, infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph			<- vector('list', length(infiles))
	df.nodestat		<- vector('list', length(infiles))	
	for(i in seq_along(infiles))
	{				
		#i<- 4
		infile			<- infiles[i]
		cat(paste('\nprocess file',i,infile))
		file			<- paste(indir.vts, '/', infile, sep='')
		#	read brl, units from annotated nexus file. attention: () may not contain two nodes
		tmp				<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
		ph				<- tmp$tree
		node.stat		<- tmp$node.stat
		node.stat		<- subset(node.stat, STAT=='Unit')
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		node.stat[, IDPOP:= as.integer(node.stat[,substr(VALUE, 4, nchar(VALUE))])]
		node.stat		<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ, IDCLU)), subset(node.stat, select=c(IDPOP, NODE_ID)), by='IDPOP')	
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw evolutionary rates for every individual in the transmission chain
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]) ), by='IDPOP')
		#	smaller ER for transmission edge
		node.stat[, BWM:= rERbw(nrow(node.stat), rERbw.args)]
		set(node.stat, NULL, 'BWM', node.stat[, ER/BWM])	#re-set to previous notation
		#	set BWM to 1 for all non-transmission edges
		tmp				<- ph$edge
		colnames(tmp)	<- c('NODE_ID','NODE_ID_TO')
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), data.table(tmp, EDGE_ID=seq_len(nrow(tmp))), by='NODE_ID') 
		setnames(tmp, c('NODE_ID','IDPOP','NODE_ID_TO'), c('NODE_ID_FROM','IDPOP_FROM','NODE_ID'))
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), tmp, by='NODE_ID') 
		tmp				<- subset( tmp, IDPOP!=IDPOP_FROM)[, NODE_ID_FROM]	#ER slows down in transmitter leading to infection. want to slow down edge leading to NODE_FROM.
		set( node.stat, node.stat[, which(!NODE_ID%in%tmp)], 'BWM', 1. )
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
		tmp				<- node.stat[, which(NODE_ID==tmp)]		
		set(node.stat, tmp, 'ER', log.df[1,meanRate] )		
		set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	set expected numbers of substitutions per branch within individual IDPOP
		setkey(node.stat, NODE_ID)
		ph$edge.length	<- ph$edge.length * node.stat[ ph$edge[, 2], ][, ER / BWM]		 
		stopifnot(all(!is.na(ph$edge.length)))		
		#	once expected number of substitutions / site are simulated, can collapse singleton nodes
		ph				<- seq.collapse.singles(ph)	
		#	set tip label so that IDPOP can be checked for consistency	
		if(pipeline.args['epi.model',][,v]=='HPTN071')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',round(DOB,d=3),label.sep,round(TIME_SEQ,d=3),sep='')]]
		if(pipeline.args['epi.model',][,v]=='DSPS')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',NA,label.sep,round(TIME_SEQ,d=3),sep='')]]		
		setkey(node.stat, NODE_ID)
		ph$tip.label	<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		#
		df.nodestat[[i]]<- node.stat
		tmp				<- ifelse(Nnode(ph), write.tree(ph, digits = 10), paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')	)
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==log.df$meanRate[1])) )
	
	if(with.plot)
	{
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM!=1.) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM==1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated within-host rate evolutionary rate\nwithout transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1] & BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.0005))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		ggsave(file, w=6, h=6)		
	}	
	#
	#	draw ancestral sequences and add to df.ph
	#
	root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
	ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
	ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	df.ph			<- cbind(df.ph, ancseq)
	#
	#	create SEQ-GEN input file
	#	
	partition.len	<- c( ncol(rANCSEQ.args$anc.seq.gag), ncol(rANCSEQ.args$anc.seq.pol), ncol(rANCSEQ.args$anc.seq.env) )
	#partition.er	<- c( 2.5, 4, 5 )
	#	split ancestral sequence into GENE and CODON_POS
	df.ph[, ANCSEQ.GAG:= substr(ANCSEQ, 1, partition.len[1])]
	df.ph[, ANCSEQ.POL:= substr(ANCSEQ, partition.len[1]+1, partition.len[1]+partition.len[2])]
	df.ph[, ANCSEQ.ENV:= substr(ANCSEQ, partition.len[1]+partition.len[2]+1, partition.len[1]+partition.len[2]+partition.len[3])]	
	stopifnot( all( strsplit( df.ph[1, ANCSEQ], '')[[1]] == strsplit( paste( df.ph[1, ANCSEQ.GAG], df.ph[1, ANCSEQ.POL], df.ph[1, ANCSEQ.ENV], sep=''), '' )[[1]] ) )
	df.ph[, ANCSEQ:=NULL]
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK'), value.name='ANCSEQ', variable.name='GENE', variable.factor=FALSE)
	set(df.ph, NULL, 'GENE', df.ph[, substr(GENE, 8, nchar(GENE))])	
	df.ph	<- df.ph[,  {
				tmp	<- strsplit(ANCSEQ, '')
				list( 	ANCSEQ.CP1= sapply(tmp, function(x) paste(x[seq.int(1,nchar(ANCSEQ[1]),3)],collapse='')  ),
						ANCSEQ.CP2= sapply(tmp, function(x) paste(x[seq.int(2,nchar(ANCSEQ[1]),3)],collapse='') ),
						ANCSEQ.CP3= sapply(tmp, function(x) paste(x[seq.int(3,nchar(ANCSEQ[1]),3)],collapse='')  ), ROOT_CALENDAR_TIME=ROOT_CALENDAR_TIME, IDCLU=IDCLU, NEWICK=NEWICK	)				
			}, by='GENE']
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK','GENE'), value.name='ANCSEQ', variable.name='CODON_POS', variable.factor=FALSE)
	set(df.ph, NULL, 'CODON_POS', df.ph[, substr(CODON_POS, 8, nchar(CODON_POS))])
	#
	#	save to file all we need to call SeqGen
	#
	file	<- paste(outdir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nsave to file=',file))
	df.seqgen	<- df.ph
	save(df.seqgen, log.df, df.nodestat, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 26-08-2014
##--------------------------------------------------------------------------------------------------------
prog.PANGEA.SeqGen.createInputFile.v1<- function()
{
	stop()
	verbose			<- 1
	with.plot		<- 1
	label.sep		<- '|'	
	#
	#	read I/O
	#
	indir.epi		<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.vts		<- '/Users/Oliver/git/HPTN071sim/tmp140914/VirusTreeSimulator'
	infile.prefix	<- '140716_RUN001_'	
	infile.args		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_PipeArgs.R'
	outdir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp140914/SeqGen'	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.vts= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.vts<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.vts= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir.sg<- tmp[1]		
		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma')%in%pipeline.args[, stat] ) )
	#
	#	setup samplers
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol			<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	rER.bwm			<- PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1(bwerm.mu=pipeline.args['bwerm.mu',][, as.numeric(v)], bwerm.sigma=pipeline.args['bwerm.sigma',][, as.numeric(v)])
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp				<- PANGEA.RootSeq.create.sampler(root.ctime.grace= 0.5, sample.grace= 3)
	rANCSEQ			<- tmp$rANCSEQ
	rANCSEQ.args	<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df			<- PANGEA.GTR.params()	
	#
	#
	#
	file		<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#	
	infiles		<- list.files(indir.vts)
	tmp			<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles		<- infiles[ grepl(tmp, infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph			<- vector('list', length(infiles))
	df.nodestat		<- vector('list', length(infiles))	
	for(i in seq_along(infiles))
	{				
		#i<- 11
		infile			<- infiles[i]
		cat(paste('\nprocess file',i,infile))
		file			<- paste(indir.vts, '/', infile, sep='')
		#	read brl, units from annotated nexus file. attention: () may not contain two nodes
		tmp				<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
		ph				<- tmp$tree
		node.stat		<- tmp$node.stat
		node.stat		<- subset(node.stat, STAT=='Unit')
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		node.stat[, IDPOP:= as.integer(node.stat[,substr(VALUE, 4, nchar(VALUE))])]
		node.stat		<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ, IDCLU)), subset(node.stat, select=c(IDPOP, NODE_ID)), by='IDPOP')	
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw evolutionary rates for every individual in the transmission chain
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]) ), by='IDPOP')
		#	smaller ER for transmission edge
		node.stat[, BWM:= rER.bwm(nrow(node.stat))]
		#	set BWM to 1 for all non-transmission edges
		tmp				<- ph$edge
		colnames(tmp)	<- c('NODE_ID','NODE_ID_TO')
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), data.table(tmp, EDGE_ID=seq_len(nrow(tmp))), by='NODE_ID') 
		setnames(tmp, c('NODE_ID','IDPOP','NODE_ID_TO'), c('NODE_ID_FROM','IDPOP_FROM','NODE_ID'))
		tmp				<- merge( subset(node.stat, select=c(NODE_ID, IDPOP)), tmp, by='NODE_ID') 
		tmp				<- subset( tmp, IDPOP!=IDPOP_FROM)[, NODE_ID_FROM]	#ER slows down in transmitter leading to infection. want to slow down edge leading to NODE_FROM.
		set( node.stat, node.stat[, which(!NODE_ID%in%tmp)], 'BWM', 1. )
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
		tmp				<- node.stat[, which(NODE_ID==tmp)]		
		set(node.stat, tmp, 'ER', log.df[1,meanRate] )		
		stopifnot( node.stat[tmp, BWM]>1)
		set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	set expected numbers of substitutions per branch within individual IDPOP
		setkey(node.stat, NODE_ID)
		ph$edge.length	<- ph$edge.length * node.stat[ ph$edge[, 2], ][, ER / BWM]
		stopifnot(all(!is.na(ph$edge.length)))
		#	once expected number of substitutions / site are simulated, can collapse singleton nodes
		ph				<- seq.collapse.singles(ph)	
		#	set tip label so that IDPOP can be checked for consistency	
		if(pipeline.args['epi.model',][,v]=='HPTN071')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',round(DOB,d=3),label.sep,round(TIME_SEQ,d=3),sep='')]]
		if(pipeline.args['epi.model',][,v]=='DSPS')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',NA,label.sep,round(TIME_SEQ,d=3),sep='')]]		
		setkey(node.stat, NODE_ID)
		ph$tip.label	<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		#
		df.nodestat[[i]]<- node.stat
		tmp				<- ifelse(Nnode(ph), write.tree(ph, digits = 10), paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')	)
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==log.df$meanRate[1])) )
	
	if(with.plot)
	{
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		ggplot(subset(df.nodestat, ER!=log.df$meanRate[1]), aes(x=ER)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, BWM!=1) , aes(x=BWM)) + geom_histogram(binwidth=0.05) + labs(x='simulated between-host rate multiplier') +
				scale_x_continuous(breaks= seq(0.8, 2.4, 0.2))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		ggsave(file, w=6, h=6)		
	}	
	#
	#	draw ancestral sequences and add to df.ph
	#
	root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
	ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
	ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	df.ph			<- cbind(df.ph, ancseq)
	#
	#	create SEQ-GEN input file
	#	
	partition.len	<- c( ncol(rANCSEQ.args$anc.seq.gag), ncol(rANCSEQ.args$anc.seq.pol), ncol(rANCSEQ.args$anc.seq.env) )
	#partition.er	<- c( 2.5, 4, 5 )
	#	split ancestral sequence into GENE and CODON_POS
	df.ph[, ANCSEQ.GAG:= substr(ANCSEQ, 1, partition.len[1])]
	df.ph[, ANCSEQ.POL:= substr(ANCSEQ, partition.len[1]+1, partition.len[1]+partition.len[2])]
	df.ph[, ANCSEQ.ENV:= substr(ANCSEQ, partition.len[1]+partition.len[2]+1, partition.len[1]+partition.len[2]+partition.len[3])]	
	stopifnot( all( strsplit( df.ph[1, ANCSEQ], '')[[1]] == strsplit( paste( df.ph[1, ANCSEQ.GAG], df.ph[1, ANCSEQ.POL], df.ph[1, ANCSEQ.ENV], sep=''), '' )[[1]] ) )
	df.ph[, ANCSEQ:=NULL]
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK'), value.name='ANCSEQ', variable.name='GENE', variable.factor=FALSE)
	set(df.ph, NULL, 'GENE', df.ph[, substr(GENE, 8, nchar(GENE))])	
	df.ph	<- df.ph[,  {
				tmp	<- strsplit(ANCSEQ, '')
				list( 	ANCSEQ.CP1= sapply(tmp, function(x) paste(x[seq.int(1,nchar(ANCSEQ[1]),3)],collapse='')  ),
						ANCSEQ.CP2= sapply(tmp, function(x) paste(x[seq.int(2,nchar(ANCSEQ[1]),3)],collapse='') ),
						ANCSEQ.CP3= sapply(tmp, function(x) paste(x[seq.int(3,nchar(ANCSEQ[1]),3)],collapse='')  ), ROOT_CALENDAR_TIME=ROOT_CALENDAR_TIME, IDCLU=IDCLU, NEWICK=NEWICK	)				
			}, by='GENE']
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK','GENE'), value.name='ANCSEQ', variable.name='CODON_POS', variable.factor=FALSE)
	set(df.ph, NULL, 'CODON_POS', df.ph[, substr(CODON_POS, 8, nchar(CODON_POS))])
	#
	#	save to file all we need to call SeqGen
	#
	file	<- paste(outdir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nsave to file=',file))
	df.seqgen	<- df.ph
	save(df.seqgen, log.df, df.nodestat, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	program to generate files for Seq Gen from output of Matt s VirusTreeSimulator
##	olli originally written 18-09-2014
##	modified 17-01-2015
##--------------------------------------------------------------------------------------------------------
#' @title Program to generate \code{SeqGen} input files
#' @description The \code{prog.PANGEA.SeqGen.createInputFile} reads files from the virus tree simulator in directory \code{indir.vts} and writes input files for \code{SeqGen}
#' to directory \code{outdir}. The program reads simulated transmission chain phylogenies with branches in units of calendar time
#' for sampled and unsampled individuals in a transmission chain. Within host evolutionary rates are drawn from a distribution, and
#' within host branch lengths are converted into the expected number of substitutions along the branch. Transmission branches are
#' multiplied with a multiplier to allow for slower evolution between hosts. The multiplier is drawn from a distribution. Starting sequences
#' are drawn from a pool of precomputed sequences. GTR parameters are drawn from a distribution. This is all that s needed to specify 
#' the SeqGen input files for each transmission chain.
#' @return NULL. Input to call SeqGen is stored in an R file.
#' @example example/ex.seqgen.inputfilecreator.R
#' @export
prog.PANGEA.SeqGen.createInputFile<- function()
{
	verbose			<- 1
	with.plot		<- 1
	label.sep		<- '|'	
	#
	#	read I/O
	#
	indir.epi		<- '/Users/Oliver/git/HPTN071sim/tmp140914/TrChains'
	infile.epi		<- '140716_RUN001_SAVE.R'	
	indir.vts		<- '/Users/Oliver/git/HPTN071sim/tmp140914/VirusTreeSimulator'
	infile.prefix	<- '140716_RUN001_'	
	infile.args		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_PipeArgs.R'
	outdir.sg		<- '/Users/Oliver/git/HPTN071sim/tmp140914/SeqGen'	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.epi= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.epi= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.epi<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,10),
									indir.vts= return(substr(arg,12,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir.vts<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.vts= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.prefix<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir.sg<- tmp[1]		
		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir.epi, infile.epi, indir.vts, infile.prefix, outdir.sg, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}	
	stopifnot( all( c('s.seed','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma')%in%pipeline.args[, stat] ) )
	#
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )	
	#
	#	setup samplers
	#
	file				<- paste(indir.epi, '/', infile.epi, sep='')
	load(file)	#expect "df.epi"    "df.trms"   "df.inds"   "df.sample"
	#
	cat(paste('\ncreate sampler of evolutionary rates'))
	#	create sampler of within host evolutionary rates
	rER.pol				<- PANGEA.WithinHostEvolutionaryRate.create.sampler.v1(wher.mu=pipeline.args['wher.mu',][, as.numeric(v)], wher.sigma=pipeline.args['wher.sigma',][, as.numeric(v)])
	#	create sampler of between host evolutionary rate dampener
	tmp					<- PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler(er.meanlog=pipeline.args['bwerm.mu',][, as.numeric(v)], er.sdlog=pipeline.args['bwerm.sigma',][, as.numeric(v)])
	rERbw				<- tmp$rERbw
	rERbw.args			<- tmp$rERbw.args
	#	create sampler of ancestral sequences
	cat(paste('\ncreate sampler of ancestral sequences'))
	tmp					<- subset( df.inds, !is.na(TIME_SEQ) )[, length(unique(IDCLU))]
	cat(paste('\nnumber of clusters for which a sequence is required, n=', tmp))
	tmp					<- ifelse(tmp<400, 0.5, ifelse(tmp<500, 1, 2))	
	cat(paste('\nUsing root.ctime.grace=', tmp))
	tmp					<- PANGEA.RootSeq.create.sampler(root.ctime.grace= tmp, sample.grace= 3 )	
	rANCSEQ				<- tmp$rANCSEQ
	rANCSEQ.args		<- tmp$rANCSEQ.args 	
	#	read GTR parameters
	log.df				<- PANGEA.GTR.params()
	if(pipeline.args['dbg.rER',][, as.numeric(v)]==1 )
	{
		cat(paste('\nSetting mus to mean per gene and codon_pos'))
		tmp				<- log.df[, list(mu= mean(mu)), by=c('GENE','CODON_POS')]
		#tmp[, ER:=mu*log.df$meanRate[1]]
		log.df			<- merge( subset(log.df, select=which(colnames(log.df)!='mu')), tmp, by=c('GENE','CODON_POS'))		
	}		
	if( pipeline.args['dbg.GTRparam',][, as.numeric(v)]==1 )
	{
		cat(paste('\nSetting GTR parameters to MEAN (except mu)'))
		tmp				<- mean	
		log.df			<- log.df[, list(state=state, mu=mu, alpha=tmp(alpha), at=tmp(at), ac=tmp(ac), cg=tmp(cg), ag=tmp(ag), gt=tmp(gt), meanRate=tmp(meanRate), a=tmp(a), c=tmp(c),  g=tmp(g), t=tmp(t) ), by=c('GENE','CODON_POS')]		
	}		
	log.df[, IDX:= seq_len(nrow(log.df))]
	log.df[, FILE:=NULL]	
	#
	#
	#	
	infiles				<- list.files(indir.vts)
	tmp					<- paste('^',infile.prefix,'.*nex$',sep='')
	infiles				<- sort(infiles[ grepl(tmp, infiles)  ])	
	if(!length(infiles))	stop('cannot find files matching criteria')
	#
	#	read from VirusTreeSimulator and convert branch lengths in time to branch lengths in subst/site
	#
	df.ph				<- vector('list', length(infiles))
	phd					<- vector('list', length(infiles))
	phs					<- vector('list', length(infiles))
	phd.plot			<- vector('list', length(infiles))
	df.nodestat			<- vector('list', length(infiles))
	cat(paste('\nUsing StartTimeMode',pipeline.args['index.starttime.mode',][,v]))
	root.edge.rate		<- NA
	if( pipeline.args['index.starttime.mode',][,v]=='shift' )
	{		
		root.edge.rate	<- 1e-6
		cat(paste('\nFix root edge rate to =',root.edge.rate))
	}		
	if( as.numeric(pipeline.args['root.edge.fixed',][,v] )) 
	{
		root.edge.rate	<- log.df[1,meanRate]
		cat(paste('\nFix root edge rate to =',root.edge.rate))		
	}	
	for(i in seq_along(infiles))
	{				
		#i	<- 40
		#i<- 47; i<- 9		
		infile			<- infiles[i]
		cat(paste('\nprocess file',i,infile))
		file			<- paste(indir.vts, '/', infile, sep='')
		#	read brl, units from annotated nexus file. attention: () may not contain two nodes
		tmp				<- hivc.beast2out.read.nexus.and.stats(file, method.node.stat='any.node')
		ph				<- tmp$tree
		node.stat		<- tmp$node.stat
		node.stat		<- subset(node.stat, STAT=='Unit')
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		node.stat[, IDPOP:= as.integer(node.stat[,substr(VALUE, 4, nchar(VALUE))])]
		node.stat		<- merge(subset(df.inds, select=c(IDPOP, GENDER, DOB, TIME_SEQ, IDCLU)), subset(node.stat, select=c(IDPOP, NODE_ID)), by='IDPOP')
		#	produce collapsed tree with branch length in units of calendar time
		phd[[i]]		<- seq.collapse.singles(ph)
		#
		#	create collapsed Newick tree with expected substitutions / site for each branch 
		#
		#	draw within host evolutionary rates for every individual in the transmission chain, and	smaller ERs along the transmission lineages
		node.stat		<- merge(node.stat, data.table( IDPOP=node.stat[, unique(IDPOP)], ER= rER.pol(node.stat[, length(unique(IDPOP))]), BWM= rERbw(node.stat[, length(unique(IDPOP))], rERbw.args) ), by='IDPOP')
		#	re-set to previous notation in terms of BWM ( between host multiplier to ER, ie ER= within-host ER / BWM )
		set(node.stat, NULL, 'BWM', node.stat[, ER/BWM])	
		#	set BWM to 1 for all edges that are NOT leading to a transmission. 
		#	Because only one seq is sampled per patient, these are only edges that end in a tip.
		set(node.stat, node.stat[, which( NODE_ID%in%seq.int(1,Ntip(ph)) )], 'BWM', 1.)
		#	no ER possible for root node - there s no edge leading to it
		set(node.stat, node.stat[, which(NODE_ID==Ntip(ph)+1)], c('ER','BWM'), NA_real_)		
		#	set root edge evolutionary rate to overall mean between-host rate
		#	get NODE_ID of edge from root
		if(!is.na(root.edge.rate))
		{
			tmp				<- ph$edge[match(Ntip(ph)+1, ph$edge[1, ]), 2]
			tmp				<- node.stat[, which(NODE_ID==tmp)]		
			set(node.stat, tmp, 'ER', root.edge.rate )		
			set(node.stat, tmp, 'BWM', 1. )		# no need to further slow down root edge			
		}
		#	check root edge length
		if( pipeline.args['index.starttime.mode',][,v]=='fix45' )
		{
			stopifnot( ph$edge.length[ which( ph$edge[, 1] == Ntip(ph)+1 ) ]>=29.5 )			
		}
		#	check calendar time of root in simulated phylogeny for consistency
		tmp				<- seq.collapse.singles(ph)
		tmp2			<- regmatches(tmp$tip.label[1], regexpr('ID_[0-9]+',tmp$tip.label[1]))
		tmp2			<- as.numeric(substr(tmp2, 4, nchar(tmp2)))
		tmp2			<- subset(node.stat, IDPOP==tmp2)[1, TIME_SEQ]
		root.ctime		<- ifelse(Nnode(tmp), tmp2 - (node.depth.edgelength(tmp)[1] + tmp$root.edge), tmp2-tmp$root.edge)		
		tmp				<- subset(node.stat, IDPOP<0)[, unique(IDPOP)]
		stopifnot(length(tmp)==1)
		stopifnot(subset(df.trms, IDTR==tmp)[, round(IDTR_TIME_INFECTED, d=1)]==round(root.ctime, d=1))
		#	check if all sampling times are consistent with node height
		tmp				<- subset( node.stat, NODE_ID<=Ntip(ph) )
		setkey(tmp, NODE_ID)
		tmp2			<- seq.collapse.singles(ph) 
		if( Nnode(tmp2) )
			tmp[, NODE_DEPTH:=root.ctime + tmp2$root.edge + node.depth.edgelength(tmp2)[ seq_len(Ntip(tmp2)) ] ]
		if( Nnode(tmp2)==0 )
			tmp[, NODE_DEPTH:=root.ctime + tmp2$root.edge ]
		stopifnot( tmp[, max(abs(NODE_DEPTH-TIME_SEQ))<=1e-6 ] )
		#	set expected numbers of substitutions per branch within individual IDPOP
		setkey(node.stat, NODE_ID)
		ph$edge.length	<- ph$edge.length * node.stat[ ph$edge[, 2], ][, ER / BWM]
		stopifnot(all(!is.na(ph$edge.length)))		
		#	once expected number of substitutions / site are simulated, can collapse singleton nodes
		ph				<- seq.collapse.singles(ph)			
		#	set tip label so that IDPOP can be checked for consistency	
		if(pipeline.args['epi.model',][,v]=='HPTN071')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',round(DOB,d=3),label.sep,round(TIME_SEQ,d=3),sep='')]]
		if(pipeline.args['epi.model',][,v]=='DSPS')
			node.stat[, LABEL:= node.stat[, paste('IDPOP_',IDPOP,label.sep,GENDER,label.sep,'DOB_',NA,label.sep,round(TIME_SEQ,d=3),sep='')]]		
		setkey(node.stat, NODE_ID)
		ph$tip.label		<- node.stat[seq_len(Ntip(ph)), ][, LABEL]
		phd[[i]]$tip.label	<- node.stat[seq_len(Ntip(phd[[i]])), ][, LABEL]
		phd.plot[[i]]		<- seq.singleton2bifurcatingtree( phd[[i]] )
		phs[[i]]			<- seq.singleton2bifurcatingtree( ph )
		#phd[[i]]			<- seq.addrootnode( phd[[i]], dummy.label=paste('NOEXIST_IDCLU',node.stat[, unique(IDCLU)],'|NA|DOB_NA|',root.ctime,sep='') )
		#
		df.nodestat[[i]]	<- node.stat
		if(Nnode(ph))
		{
			tmp			<- write.tree(ph, digits = 10)			
			tmp			<- paste( '(',substr(tmp,1,nchar(tmp)-1),',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')
			phd[[i]]	<- write.tree(phd[[i]], digits = 10)
		}
		if(!Nnode(ph))
		{
			tmp			<- paste( '(',ph$tip.label,':',ph$root.edge,',NOEXIST_NA|NA|DOB_NA|',root.ctime,':0):0;', sep='')
			phd[[i]]	<- paste( '(',phd[[i]]$tip.label,':',phd[[i]]$root.edge, ');', sep='' )
		}
		df.ph[[i]]		<- data.table(ROOT_CALENDAR_TIME= root.ctime, IDCLU=node.stat[, unique(IDCLU)], NEWICK=tmp)
		#readline()
	}
	#	write cluster trees of each transmission network into a single newick file
	phd			<- sapply(phd,'[[',1)
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDCLUTREES.newick', sep='')
	cat(paste('\nWrite to file=',file))
	writeLines(phd, con=file)
	#	get multifurcating tree with brl in units of calendar time
	options(expressions=1e4)
	phd			<- eval(parse(text=paste('phd.plot[[',seq_along(phd.plot),']]', sep='',collapse='+')))
	options(expressions=5e3)
	phd			<- drop.tip(phd, which(grepl('DUMMY', phd$tip.label)), root.edge=1)
	phd			<- ladderize(phd)
	#	write multifurcating tree with brl in units of calendar time
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDTREE.newick', sep='')
	cat(paste('\nWrite to file=',file))
	write.tree(phd, file=file)
	#	write multifurcating tree with brl in units of subst/site
	options(expressions=1e4)
	phs			<- eval(parse(text=paste('phs[[',seq_along(phs),']]', sep='',collapse='+')))
	options(expressions=5e3)
	phs			<- drop.tip(phs, which(grepl('DUMMY', phs$tip.label)), root.edge=1)
	phs			<- ladderize(phs)
	file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'SUBSTTREE.newick', sep='')
	cat(paste('\nWrite to file=',file))
	write.tree(phs, file=file)
	#
	df.ph		<- do.call('rbind',df.ph)
	df.nodestat	<- do.call('rbind',df.nodestat)
	#	
	#	
	#	check that we have exactly one root edge with overall mean between host rate per cluster
	if(!is.na(root.edge.rate))
		stopifnot( df.nodestat[, length(unique(IDCLU))]==nrow(subset(df.nodestat, ER==root.edge.rate)) )	
	
	#tmp	<- unique(subset(df.inds, IDPOP>=-110 & IDPOP<0, IDCLU))
	#tmp	<- merge(df.inds, tmp, by='IDCLU')[, length(which(!is.na(TIME_SEQ)))] / df.inds[, length(which(!is.na(TIME_SEQ)))]
	#cat(paste('\nProportion of sequences descending from no import after baseline=', tmp))

	#df.nodestat[, length(which(!is.na(TIME_SEQ))), by='IDCLU']
	#subset(df.nodestat, IDCLU==72)
	#subset(df.ph, IDCLU==72)
	#subset(df.ph, IDCLU==47)	
	if(with.plot)
	{
		#subset(df.nodestat, ER!=root.edge.rate)[, table(BWM==1.)]
		#ggplot(subset(df.nodestat, ER!=root.edge.rate) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)
		#ggplot(subset(df.nodestat, ER!=root.edge.rate & BWM!=1.) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001)
		#ggplot(subset(df.nodestat, ER!=root.edge.rate), aes(x=ER, y=BWM)) + geom_point()	
		#	plot used within-host ERs
		tmp		<- root.edge.rate
		if(is.na(tmp))
			tmp	<- Inf
		ggplot(subset(df.nodestat, ER!=tmp), aes(x=ER/BWM)) + geom_histogram(binwidth=0.001)	+ labs(x='simulated within-host evolutionary rate') +
			scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_ER.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)
		#	plot used between host modifiers
		ggplot(subset(df.nodestat, ER!=tmp & BWM==1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.001) + labs(x='simulated within-host rate evolutionary rate\nwithout transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.002))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWM.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)
		#	plot used ERs along transmission edges
		ggplot(subset(df.nodestat, ER!=tmp & BWM!=1) , aes(x=ER/BWM)) + geom_histogram(binwidth=0.0001) + labs(x='simulated evolutionary rates along transmission edges') +
				scale_x_continuous(breaks= seq(0, 0.02, 0.0005))
		file	<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'INFO_sg_BWER.pdf', sep='')
		cat(paste('\nWrite to file=',file))
		ggsave(file, w=6, h=6)	
		
		if(grepl('fix',pipeline.args['index.starttime.mode',][,v]))
		{
			setkey(df.nodestat, IDPOP)
			tmp				<- merge(subset(unique(df.nodestat), select=c(IDCLU, LABEL)), unique(df.nodestat)[, list(IDCLU_N=length(which(!is.na(TIME_SEQ)))), by='IDCLU'], by='IDCLU')
			tmp[, LABEL_NEW:= tmp[, paste(IDCLU,'_',IDCLU_N,label.sep,LABEL,sep='')]]
			setkey(tmp, LABEL)			
			phd$tip.label	<- tmp[ phd$tip.label, ][, LABEL_NEW]
			file			<- paste(indir.epi, '/', substr(infile.epi,1,nchar(infile.epi)-6),'DATEDTREE.pdf', sep='')
			cat(paste('\nWrite to file=',file))
			pdf(file=file, w=10, h=Ntip(phd)*0.1)
			plot(phd, show.tip=TRUE, cex=0.5)			
			dev.off()
		}
	}
	#
	#	draw ancestral sequences and add to df.ph
	#
	if(pipeline.args['startseq.mode',][,v]=='many')
	{
		root.ctime		<- df.ph[, ROOT_CALENDAR_TIME]
		ancseq			<- rANCSEQ(root.ctime, rANCSEQ.args)
		ancseq			<- data.table(ANCSEQ= apply(as.character(ancseq),1,function(x) paste(x, collapse='')) )		
	}		
	if(pipeline.args['startseq.mode',][,v]=='one')
	{
		cat(paste('\nStartSeqModel=',pipeline.args['startseq.mode',][,v],'use first sampled starting sequence for all' ))
		stopifnot(max(abs(df.ph[1, ROOT_CALENDAR_TIME]-df.ph[, ROOT_CALENDAR_TIME]))<100*EPS)
		root.ctime		<- df.ph[1, ROOT_CALENDAR_TIME]
		tmp				<- rANCSEQ(root.ctime, rANCSEQ.args)
		tmp				<- data.table(ANCSEQ= apply(as.character(tmp),1,function(x) paste(x, collapse='')) )
		ancseq			<- data.table(ANCSEQ=rep(NA_character_, nrow(df.ph)))
		set(ancseq, NULL, 'ANCSEQ', tmp[1,ANCSEQ])
	}
	df.ph			<- cbind(df.ph, ancseq)
	#
	#	create SEQ-GEN input file
	#	
	partition.len	<- c( ncol(rANCSEQ.args$anc.seq.gag), ncol(rANCSEQ.args$anc.seq.pol), ncol(rANCSEQ.args$anc.seq.env) )
	#partition.er	<- c( 2.5, 4, 5 )
	#	split ancestral sequence into GENE and CODON_POS
	df.ph[, ANCSEQ.GAG:= substr(ANCSEQ, 1, partition.len[1])]
	df.ph[, ANCSEQ.POL:= substr(ANCSEQ, partition.len[1]+1, partition.len[1]+partition.len[2])]
	df.ph[, ANCSEQ.ENV:= substr(ANCSEQ, partition.len[1]+partition.len[2]+1, partition.len[1]+partition.len[2]+partition.len[3])]	
	stopifnot( all( strsplit( df.ph[1, ANCSEQ], '')[[1]] == strsplit( paste( df.ph[1, ANCSEQ.GAG], df.ph[1, ANCSEQ.POL], df.ph[1, ANCSEQ.ENV], sep=''), '' )[[1]] ) )
	df.ph[, ANCSEQ:=NULL]
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK'), value.name='ANCSEQ', variable.name='GENE', variable.factor=FALSE)
	set(df.ph, NULL, 'GENE', df.ph[, substr(GENE, 8, nchar(GENE))])	
	df.ph	<- df.ph[,  {
				tmp	<- strsplit(ANCSEQ, '')
				list( 	ANCSEQ.CP1= sapply(tmp, function(x) paste(x[seq.int(1,nchar(ANCSEQ[1]),3)],collapse='')  ),
						ANCSEQ.CP2= sapply(tmp, function(x) paste(x[seq.int(2,nchar(ANCSEQ[1]),3)],collapse='') ),
						ANCSEQ.CP3= sapply(tmp, function(x) paste(x[seq.int(3,nchar(ANCSEQ[1]),3)],collapse='')  ), ROOT_CALENDAR_TIME=ROOT_CALENDAR_TIME, IDCLU=IDCLU, NEWICK=NEWICK	)				
			}, by='GENE']
	df.ph	<- melt(df.ph, id.var=c('ROOT_CALENDAR_TIME','IDCLU','NEWICK','GENE'), value.name='ANCSEQ', variable.name='CODON_POS', variable.factor=FALSE)
	set(df.ph, NULL, 'CODON_POS', df.ph[, substr(CODON_POS, 8, nchar(CODON_POS))])
	#
	#	save to file all we need to call SeqGen
	#
	file	<- paste(outdir.sg,'/',infile.prefix, 'seqgen.R',sep='')
	cat(paste('\nsave to file=',file))
	df.seqgen	<- df.ph
	save(df.seqgen, log.df, df.nodestat, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Arguments for rPANGEAHIV simulation pipeline
#' @description Construct table with all input arguments to the rPANGEAHIV simulation pipeline.
#' The demographic within host coealescent model is N0tau*(1+exp(-r*T50))/(1+exp(-r*(T50-t))). Default parameters are set
#' so that the curve asymptotes at Ne*tau=3e5 and reaches half its asymptotic value 1 year post infection. Branch lengths are 
#' multiplied by within host evolutionary rates, and within host branch lengths ending in a transmission event are multiplied
#' in addition with a between-host evolutionary rate multiplier. 
#' @param yr.start				Start year of the epi simulation (default 1980)
#' @param yr.end				First year after the epi simulation (default 2020)
#' @param seed					Random number seed
#' @param s.MODEL				Sampling model to use
#' @param s.INC.recent			Proportion of incident cases sampled (not used)
#' @param s.INC.recent.len		Number of last year in which an exact proportion of incident cases is sampled (not used) 
#' @param s.PREV.min			Proportion of infected cases sampled at start of the simulation
#' @param s.PREV.max			Proportion of infected cases sampled at the end of the simulation
#' @param s.PREV.max.n			Number of infected cases sampled (usually NA, only used for Prop2Untreated)
#' @param s.INTERVENTION.prop	Proportion of infected cases that are sampled from after intervention start
#' @param s.INTERVENTION.start	Year in which the community intervention starts
#' @param s.INTERVENTION.mul	Multiplier of number of sequences sampled per year after start of the intervention
#' @param s.ARCHIVAL.n			Total number of sequences sampled before diagnosis
#' @param epi.model				The epi model used to create the epi simulation (default HPTN071, alternatively DSPS)
#' @param epi.dt				Time increment of the epi simulation (default 1/48)
#' @param epi.import			Proportion of imported cases of all cases (default 0.1)
#' @param root.edge.fixed		Boolean; fix evolutionary rate of root edges to mean rate if true 
#' @param v.N0tau				Parameter of BEAST::LogisticGrowthN0::N0 (default 3.58e4) 
#' @param v.r					Parameter of BEAST::LogisticGrowthN0::r (default 2)
#' @param v.T50					Parameter of BEAST::LogisticGrowthN0::T50
#' @param wher.mu				Mean within host evolutionary rate of log normal density
#' @param wher.sigma			Standard deviation in within host evolutionary rate of log normal density
#' @param bwerm.mu				Mean between host evolutionary rate multiplier of log normal density
#' @param bwerm.sigma			Standard deviation in between host evolutionary rate multiplier of log normal density
#' @param sp.prop.of.sexactive	Proportion of population sampled in seroprevalence survey
#' @param report.prop.recent	Proportion of individuals for whom recency of infection should be reported
#' @param dbg.GTRparam			debug flag
#' @param dbg.rER				debug flag
#' @param startseq.mode			Number of different starting sequences either 'many' or 'one'.
#' @param index.starttime.mode	distribution to sample times for starting sequence: Normal(1960,7) or Unif(1959.75, 1960.25)
#' @return data.table
#' @example example/ex.pipeline.HPTN071.R
#' @example example/ex.pipeline.DSPS.R
#' @export
rPANGEAHIVsim.pipeline.args<- function(	yr.start=1980, yr.end=2020, seed=42,
										s.INC.recent=0.1, s.INC.recent.len=5, s.PREV.min=0.01, s.PREV.max=0.25, s.PREV.base=exp(1), s.INTERVENTION.start=2015, 
										s.INTERVENTION.mul= 2, s.ARCHIVAL.n=50, s.MODEL='Prop2DiagB4I', s.PREV.max.n=NA, s.INTERVENTION.prop=NA,
										epi.model='HPTN071', epi.dt=1/48, epi.import=0.1, root.edge.fixed=0,
										v.N0tau=3.58e4, v.r=2, v.T50=-1,
										wher.mu=0.005, wher.sigma=0.8,
										bwerm.mu=1.5, bwerm.sigma=0.12, er.gamma=4, 
										sp.prop.of.sexactive= 0.05, 
										report.prop.recent=0.2,
										dbg.GTRparam=0, dbg.rER=0, 
										index.starttime.mode='normal', startseq.mode='many', seqtime.mode='Gamma9')									
{
	#explore within host Neff*tau model
	if(0)
	{
		#visualize dependence of size after 4 years on growth rate
		r		<- seq(1,16,0.1)
		plot(r, ( 1+exp(r) ) / ( 1+exp(-3*r) ), type='l')
		#estimate growth rates for desired sizes	T50=-1
		f		<- function(r, b){	abs(( 1+exp(r) ) / ( 1+exp(-3*r) )-b)	}		
		optimize(f=f, interval=c(2, 15), b=3e5 )	#12.61152
		optimize(f=f, interval=c(2, 15), b=3e4 )	#10.30891
		optimize(f=f, interval=c(2, 15), b=3e3 )	#8.006034
		optimize(f=f, interval=c(2, 15), b=3e2 )	#5.70044
		#estimate growth rates for desired sizes	T50=-2
		f		<- function(r, b){	abs(( 1+exp(2*r) ) / ( 1+exp(-2*r) )-b)	}		
		optimize(f=f, interval=c(2, 15), b=3e5 )	#6.305779
		optimize(f=f, interval=c(2, 15), b=3e4 )	#5.154461
		optimize(f=f, interval=c(2, 15), b=3e3 )	#4.003191
		optimize(f=f, interval=c(2, 15), b=3e2 )	#2.851904		
		f		<- function(r, b){	abs(( 1+exp(2*r) ) / ( 1+exp(-r*(-2+10)) )-b)	}
		optimize(f=f, interval=c(2, 15), b=3e2 )	#2.850242
		optimize(f=f, interval=c(2, 15), b=5e4 )	#5.409859
		optimize(f=f, interval=c(2, 15), b=2e5 )	#6.103021
		optimize(f=f, interval=c(2, 15), b=150 )	#2.501955
		optimize(f=f, interval=c(0.1, 15), b=10 )	#1.098685
		
		Net			<- function(t, N0tau, r, T50){  N0tau*(1+exp(-r*T50))/(1+exp(-r*(T50-t)))	}		
		x			<- seq(-10,0,0.001)
		#tmp			<- data.table(x=x, y5=Net(x, 1, 12.61152, -1), y4=Net(x, 1, 10.30891, -1), y3=Net(x, 1, 8.006034, -1), y2=Net(x, 1, 5.70044, -1))
		tmp			<- data.table(x=x, y5=Net(x, 1, 6.305779, -2), y4=Net(x, 1, 5.154461, -2), y3=Net(x, 1, 4.003191, -2), y2=Net(x, 1, 2.850242, -2), y1=Net(x, 1, 1.098685, -2))
		tmp			<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_y_log10(breaks=c(3e2,3e3,3e4,3e5)) + scale_x_continuous(breaks=seq(-20,10,1))
	}
	#	plot increasing sampling fraction
	if(0)
	{
		s.PREV.base	<- exp(1)
		yr.start	<- 1980
		yr.end		<- 2020
		s.PREV.min	<- 0.01
		df.sample	<- as.data.table(expand.grid(YR=seq.int(yr.start, yr.end), s.PREV.max=c(0.11, 0.15, 0.185)))		
		#	exponential rate of increasing s.TOTAL (total sampling rate) per year
		set(df.sample, NULL, 'r', df.sample[, log( s.PREV.max/s.PREV.min, base=s.PREV.base ) / diff(range(YR)) ] )
		set(df.sample, NULL, 's.fraction', df.sample[, s.PREV.base^( r*(YR-min(YR)) ) * s.PREV.min ] )
		set(df.sample, NULL, 's.PREV.max', df.sample[, factor(s.PREV.max, levels=c('0.11','0.15','0.185'), labels=c('A','B','C'))]) 			
		ggplot(df.sample, aes(x=YR, y=s.fraction, colour=s.PREV.max)) + geom_line() + scale_y_continuous(breaks=seq(0, 0.16, 0.02)) + labs(colour='scenario', x='year', y='sampling fraction')		
	}
	pipeline.args	<- data.table(	stat= 	c('yr.start','yr.end','s.MODEL','s.INC.recent','s.INC.recent.len', 's.PREV.min', 's.PREV.max', 's.PREV.max.n', 's.INTERVENTION.prop', 's.INTERVENTION.start', 's.INTERVENTION.mul', 's.ARCHIVAL.n', 's.seed', 'index.starttime.mode', 'startseq.mode', 'seqtime.mode','root.edge.fixed','epi.model', 'epi.dt', 'epi.import','v.N0tau','v.r','v.T50','wher.mu','wher.sigma','bwerm.mu','bwerm.sigma','er.gamma','sp.prop.of.sexactive','report.prop.recent','dbg.GTRparam','dbg.rER'), 
									v	=	c(yr.start, yr.end, s.MODEL, s.INC.recent, s.INC.recent.len, s.PREV.min, s.PREV.max, s.PREV.max.n, s.INTERVENTION.prop, s.INTERVENTION.start, s.INTERVENTION.mul, s.ARCHIVAL.n, seed, index.starttime.mode, startseq.mode, seqtime.mode, root.edge.fixed, epi.model, epi.dt, epi.import, v.N0tau, v.r, v.T50, wher.mu, wher.sigma, bwerm.mu, bwerm.sigma, er.gamma, sp.prop.of.sexactive, report.prop.recent, dbg.GTRparam, dbg.rER) )
	setkey(pipeline.args, stat)	
	pipeline.args
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 06-08-2015
##--------------------------------------------------------------------------------------------------------
pipeline.various<- function()
{
	if(1)	#align sequences in fasta file with Clustalo
	{
		cmd			<- cmd.various()
		cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=771, hpc.mem="5000mb")
		cat(cmd)		
		outdir		<- paste(HOME,"tmp",sep='/')
		outfile		<- paste("vrs",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
		cmd.hpccaller(outdir, outfile, cmd)
		quit("no")		
	}	
	if(0)
	{
		mfile		<- paste(DATA,'model_150816a.R',sep='/')
		indir.st	<- paste(DATA,'contigs_150408_wref_cutstat',sep='/')
		indir.al	<- paste(DATA,'contigs_150408_wref',sep='/')
		outdir		<- paste(DATA,'contigs_150408_model150816a',sep='/')
		trainfile	<- paste(DATA,'contigs_150408_trainingset_subsets.R',sep='/')
		batch.n		<- 200
		for(batch.id in seq.int(1,14))
		{			
			cmd			<- cmd.haircut.call(indir.st, indir.al, outdir, mfile, trainfile=trainfile, batch.n=batch.n, batch.id=batch.id, prog=PR.HAIRCUT.CALL )
			cmd			<- cmd.hpcwrapper(cmd, hpc.nproc= 1, hpc.q='pqeelab', hpc.walltime=21, hpc.mem="5000mb")
			cat(cmd)		
			outdir		<- paste(HOME,"tmp",sep='/')
			outfile		<- paste("cntcall",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')
			cmd.hpccaller(outdir, outfile, cmd)				
		}	
		quit("no")
	}
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 07-11-2015
##--------------------------------------------------------------------------------------------------------
prog.treecomparison.metrics<- function()
{
	file	<- '/work/or105/Gates_2014/tree_comparison/submitted_151101.rda'
	#treedist.quartets.add(file=file, with.save=1)
	treedist.billera.add(file=file, with.save=1)
	quit("no")
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title rPANGEAHIV simulation pipeline
#' @description Reads two files \code{infile.ind} and \code{infile.trm} from the epi simulator in directory \code{indir} and produces a UNIX batch
#' file that contains all the simulation steps in directory \code{outdir}. 
#' @param indir				Input directory
#' @param infile.ind		Input file with individual metavariables
#' @param infile.trm		Input file with transmission events
#' @param pipeline.args		Input arguments for the simulation
#' @return file name of qsub or UNIX batch file
#' @example example/ex.pipeline.HPTN071.R
#' @example example/ex.pipeline.DSPS.R
#' @export
rPANGEAHIVsim.pipeline<- function(indir, infile.ind, infile.trm, outdir, pipeline.args=rPANGEAHIVsim.pipeline.args() )
{
	verbose			<- 1
	infile.args		<- paste(outdir,'/',substr(infile.trm, 1, nchar(infile.trm)-7), 'PipeArgs.R',sep='')
	save(pipeline.args, file=infile.args)
	#
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, sep='\n'))
	}	
	#
	#	pipeline start
	#	
	##	sample sequences and draw imports 
	cmd				<- "#######################################################
#######################################################
#######################################################
#
# start: run rPANGEAHIVsim.pipeline
#
#######################################################
#######################################################
#######################################################"		
	outdir.TrChain	<- paste(outdir,'/TrChains',sep='')
	cmd				<- paste(cmd, '\nmkdir -p ', outdir.TrChain,sep='')
	if(pipeline.args['epi.model'][,v]=='HPTN071')
	{
		cmd			<- paste(cmd, cmd.HPTN071.input.parser.v4(indir, infile.trm, infile.ind, infile.args, outdir.TrChain,  infile.trm, infile.ind), sep='\n')	
	}	
	if(pipeline.args['epi.model'][,v]=='DSPS')
	{
		infile.ind	<- gsub('TRM','IND',infile.trm)
		cmd			<- paste(cmd, cmd.DSPS.input.parser.v2(indir, infile.trm, infile.args, outdir.TrChain,  infile.trm, infile.ind), sep='\n')		
	}
	##	run virus tree simulator
	outdir.VTS		<- paste(outdir,'/VirusTreeSimulator',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', outdir.VTS,sep='')
	outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
	prog.args		<- paste('-seed ',pipeline.args['s.seed',][, v],' -demoModel Logistic -N0 ',pipeline.args['v.N0tau',][,v] ,' -growthRate ', pipeline.args['v.r',][,v],' -t50 ',pipeline.args['v.T50',][,v], sep='')	
	cmd				<- paste(cmd, cmd.VirusTreeSimulator(outdir.TrChain, infile.trm, infile.ind, outdir.VTS, outfile, prog.args=prog.args), sep='\n')	
	##	create seq gen input files 
	outdir.SG		<- paste(outdir,'/SeqGen',sep='')
	cmd				<- paste(cmd, 'mkdir -p ', outdir.SG,sep='')
	infile.epi		<- paste( substr(infile.ind, 1, nchar(infile.ind)-7),'SAVE.R', sep='' )
	infile.vts		<- substr(infile.ind, 1, nchar(infile.ind)-7)
	cmd				<- paste(cmd, cmd.SeqGen.createInputFiles(outdir.TrChain, infile.epi, outdir.VTS, infile.vts, infile.args, outdir.SG), sep='\n')
	##	run SeqGen	
	outfile			<- substr(infile.ind, 1, nchar(infile.ind)-7)
	cmd				<- paste(cmd, cmd.SeqGen.run(outdir.TrChain, infile.epi, outdir.SG, outfile, infile.args, outdir), sep='')
	##	clean up
	cmd				<- paste(cmd,'rm -rf ',outdir.TrChain,' ', outdir.VTS,' ', outdir.SG,'\n',sep='')
	cmd				<- paste(cmd,"#######################################################
#######################################################
#######################################################
#
# end: run rPANGEAHIVsim.pipeline
#
#######################################################
#######################################################
#######################################################\n",sep='')
	if(verbose)
		cat(cmd)
	outfile			<- paste("pngea",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='.')	
	cmd.hpccaller(outdir, outfile, cmd)	
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 26-07-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 1)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v1.R
#' @export
prog.HPTN071.input.parser.v1<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- '/Users/Oliver/git/HPTN071sim/raw_trchain'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.txt'
	infile.args		<- NA
	outdir			<- '/Users/Oliver/git/HPTN071sim/sim_trchain'
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#	
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]		
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=' ', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))		
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(TIME_TR>pipeline.args['yr.start',][, as.numeric(v)] & length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(pipeline.args['yr.end',][, as.numeric(v)], z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')
	set(df.trm, NULL, 'YR', df.trm[, floor(TIME_TR)])
	#	check that all transmission times except baseline are unique
	tmp		<- subset(df.trm, TIME_TR>pipeline.args['yr.start',][, as.numeric(v)])
	stopifnot( nrow(tmp)==tmp[,length(unique(TIME_TR))] )
	
	
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c("Id","Gender","DoB","DateOfDeath","RiskGroup","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','CIRCM'))
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', pipeline.args['yr.end',][, as.numeric(v)]+1.)		
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	
	
	# compute prevalence and incidence by year	
	df.epi		<- df.trm[, list(INC=length(IDREC)), by='YR']
	tmp			<- df.epi[, 	{
				alive		<- which( floor(df.ind[['DOB']])<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  ceiling(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(alive), PREV=length(infected))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/POP])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/POP])			
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)] ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, exp( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])	
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences sampled=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences sampled=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- df.trm[, {
				tmp<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
				tmp<- sample(seq_along(IDREC), tmp)
				list( 	IDPOP=IDREC[tmp], TIME_TR=TIME_TR[tmp], 
						TIME_SEQ=TIME_TR[tmp]+rexp(length(tmp), rate=1/(3*30))/365, 
						INCIDENT_SEQ=rep('Y',length(tmp) ) )
			}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)			
	#	sample non-incident cases by year
	for(yr in df.sample[, YR][-1])
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr & floor(TIME_TR)<yr)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		tmp2	<- sample(seq_len(nrow(tmp)), tmp2)
		#	set variables in df.inds
		tmp		<- data.table(IDPOP= tmp[tmp2, IDPOP], TIME_SEQ=runif(length(tmp2), min=yr, max=yr+1), INCIDENT_SEQ=rep('N',length(tmp2) ))
		cat(paste('\nsampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- sapply(tmp[,IDPOP], function(x) df.inds[,which(IDPOP==x)])
		set(df.inds, tmp2, 'TIME_SEQ', tmp[,TIME_SEQ])
		set(df.inds, tmp2, 'INCIDENT_SEQ', tmp[,INCIDENT_SEQ])		
	}
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	#
	#	check that allocation OK
	#
	set(df.inds, NULL, 'TIME_SEQYR', df.inds[, floor(TIME_SEQ)])
	tmp	<- subset(df.inds, !is.na(TIME_SEQ))[, list(s.n.TOTAL=length(IDPOP)), by='TIME_SEQYR']
	setkey(tmp, TIME_SEQYR)
	set(tmp,NULL,'s.n.CUMTOTAL',tmp[, cumsum(s.n.TOTAL)])
	stopifnot(  tmp[,tail(s.n.CUMTOTAL,1)]==df.sample[, tail(s.n.CUMTOTAL,1)] ) 
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	TRANSMISSION NETWORKS
	#
	require(igraph)
	#	need a unique index number for every cluster
	setkey(df.trms, TIME_TR)
	tmp			<- df.trms[, which(IDTR<0)]
	set(df.trms, tmp, 'IDTR', rev(-seq_along(tmp)))
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 's.n.INC', df.sample[, as.integer(s.n.INC)])
		set(df.sample, NULL, 's.n.notINC', df.sample[, as.integer(s.n.notINC)])
		set(df.sample, NULL, 's.n.TOTAL', df.sample[, as.integer(s.n.TOTAL)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','s.n.TOTAL','s.n.INC','s.n.notINC'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','s.n.TOTAL','s.n.INC','s.n.notINC'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,pipeline.args['yr.end',][, as.numeric(v)],2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		ggsave(file=file, w=16, h=8)
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu, select=c(IDPOP, GENDER, TIME_SEQ))
					tmp[, IS_SEQ:= tmp[, factor(!is.na(TIME_SEQ), label=c('N','Y'), levels=c(FALSE, TRUE))]]
					clu.igr				<- graph.data.frame(subset(df.trms, IDCLU==clu & IDTR>=0, select=c(IDTR, IDREC)), directed=TRUE, vertices=subset(tmp, select=c(IDPOP, GENDER, IS_SEQ)))
					V(clu.igr)$color	<- ifelse( get.vertex.attribute(clu.igr, 'IS_SEQ')=='Y', 'green', 'grey90' )
					V(clu.igr)$shape	<- ifelse( get.vertex.attribute(clu.igr, 'GENDER')=='M', 'circle', 'square' )
					
					par(mai=c(0,0,1,0))
					plot(clu.igr, main=paste('IDCLU=',clu,sep=''), vertex.size=2, vertex.label.cex=0.25, edge.arrow.size=0.5, layout=layout.fruchterman.reingold(clu.igr, niter=1e3) )
					legend('bottomright', fill=c('green','grey90'), legend=c('sequence sampled','sequence not sampled'), bty='n')
					legend('bottomleft', legend=c('square= Female','circle= Male'), bty='n')				
				})
		dev.off()	
	}	
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample)
	#	save for virus tree simulator
	#	exclude columns that are not needed	
	df.inds	<- subset(df.inds, !is.na(TIME_TR))
	df.inds[, RISK:=NULL]
	df.inds[, INCIDENT_SEQ:=NULL]
	df.inds[, TIME_SEQYR:=NULL]	
	df.trms[, TR_ACUTE:=NULL]
	df.trms[, YR:=NULL]	
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title DSPS parser (version 2, includes simulation of imports)
#' @description Reads files from the DSPS epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.DSPS.v2.R
#' @export
prog.DSPS.input.parser.v2<- function()
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")	
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140911'
	infile.trm		<- '140911_DSPS_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140911_DSPS_RUN001_IND.csv'
	outfile.trm		<- '140911_DSPS_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	#infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmission data.table
	#
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=',', dec='.'))
	#	for now ignore SAMPLING: use sampling model consistent with HPTN071 sampling model
	df.trm	<- subset(df.trm, EventType%in%c('INFECTION','DEATH','BIRTH'))
	#	rescale time 
	set(df.trm, NULL, 'ActionTime', df.trm[, ActionTime+1950])
	#	since - is confusing with -1 change to _
	set(df.trm, NULL, 'FromDeme.FromHost', df.trm[, gsub('-','_',FromDeme.FromHost)])
	set(df.trm, NULL, 'ToDeme.ToHost', df.trm[, gsub('-','_',ToDeme.ToHost)])
	tmp		<- df.trm[, which(is.na(FromDeme.FromHost))]
	cat(paste('\nFound TRM entries with NA FromDeme.FromHost, n=', length(tmp)))
	set(df.trm, tmp, 'FromDeme.FromHost', paste('HouseNA_',-seq_along(tmp),'_NA',sep=''))
	#	get ID of transmitters
	set(df.trm, NULL, 'IDTR', as.numeric(df.trm[, sapply( strsplit(FromDeme.FromHost,'_'), '[[', 2 )]))
	#	get ID of recipients
	set(df.trm, NULL, 'IDREC', as.numeric(df.trm[, sapply( strsplit(ToDeme.ToHost,'_'), '[[', 2 )]))
	#	
	#	prepare metavariable data.table
	#
	df.ind	<- subset(df.trm, select=c(IDREC, ToDeme.ToHost, EventType, ActionTime))
	setnames(df.ind, c('IDREC','ToDeme.ToHost'), c('IDPOP','INFO'))
	tmp		<- subset(df.trm, select=c(IDTR, FromDeme.FromHost))
	setnames(tmp, c('IDTR','FromDeme.FromHost'), c('IDPOP','INFO'))
	df.ind	<- rbind(df.ind, tmp, useNames=TRUE, fill=TRUE)
	df.ind	<- df.ind[-nrow(df.ind), ]
	df.ind[, V1:=NULL]
	#	removing unrealistic IDPOP
	cat(paste('\nFound IDPOP>1e6, remove, n=',nrow(subset(df.ind, IDPOP>1e6)),'2Emma: these are birth events'))
	df.ind	<- subset(df.ind, IDPOP<1e6)
	#	get House ID
	set(df.ind, NULL, 'HOUSE', df.ind[, sapply( strsplit(INFO,'_'), '[[', 1 )])
	suppressWarnings( set(df.ind, NULL, 'HOUSE', df.ind[, as.integer(substr(HOUSE,6,nchar(HOUSE)))]) )
	#	get GENDER
	set(df.ind, NULL, 'GENDER', df.ind[, sapply( strsplit(INFO,'_'), '[[', 3 )])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER, levels=c('FEMALE','MALE'), labels=c('F','M'))])	
	#	add DOD	date of death if any
	tmp		<- subset(df.ind, EventType=='DEATH', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='DEATH'), unique(tmp), by='IDPOP', all.x=TRUE)	
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','DOD'))
	#	add DOB	date of birth if any
	tmp		<- subset(df.ind, EventType=='BIRTH', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='BIRTH'), unique(tmp), by='IDPOP', all.x=TRUE)
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','DOB'))
	#	add time of infection if any
	tmp		<- subset(df.ind, EventType=='INFECTION', select=c(IDPOP, ActionTime))
	setkey(tmp, IDPOP)
	#tmp		<- subset(df.trm, EventType=='INFECTION', select=c(IDREC, ActionTime) )
	#setnames(tmp, c('IDREC','ActionTime'), c('IDPOP', 'TIME_TR'))	
	#merge( subset(df.ind, is.na(EventType) | EventType!='INFECTION'), unique(tmp), by='IDPOP', all.x=TRUE, all.y=TRUE)
	df.ind	<- merge( subset(df.ind, is.na(EventType) | EventType!='INFECTION'), unique(tmp), by='IDPOP', all.x=TRUE, all.y=TRUE)
	setnames(df.ind, c('ActionTime.x','ActionTime.y'), c('ActionTime','TIME_TR'))	
	#	nothing else to process
	stopifnot( df.ind[, all(is.na(EventType))] )
	df.ind[, EventType:=NULL]
	df.ind[, ActionTime:=NULL]	
	#
	#	clean up df.ind
	#
	setkey(df.ind, IDPOP)	
	df.ind	<- unique(df.ind)
	stopifnot( nrow(subset(df.ind, DOD<=DOB))==0 )
	stopifnot( nrow(subset(df.ind, TIME_TR<=DOB))==0 )
	set(df.ind, df.ind[, which(TIME_TR<=DOB)], 'DOB', NA_real_)
	stopifnot( nrow(subset(df.ind, TIME_TR>=DOD))==0 )
	df.ind[, INFO:=NULL]
	set( df.ind, NULL, 'IDPOP', df.ind[, as.integer(IDPOP)] )
	df.ind	<- subset( df.ind, !is.na(DOB) | !is.na(DOD) | !is.na(TIME_TR) | IDPOP<0 )
	tmp		<- df.ind[, which(IDPOP>0 & is.na(DOD))]
	cat(paste('\nFound individuals alive at simulation end, n=', length(tmp)))	
	set(df.ind, tmp, 'DOD', df.ind[, ceiling(max(max(DOD, na.rm=TRUE), max(TIME_TR, na.rm=TRUE))+1)] )
	tmp		<- df.ind[, which(IDPOP>0 & is.na(DOB))]
	cat(paste('\nFound individuals with no birth date, n=', length(tmp)))
	set(df.ind, tmp, 'DOB', df.ind[, floor(min(TIME_TR, na.rm=TRUE))-1] )
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	cat(paste('\nFound index cases, n=', nrow(subset(df.ind,IDPOP<0))))
	#
	#	clean up df.trm
	#
	df.trm	<- subset(df.trm, EventType=='INFECTION')
	setnames(df.trm, c("ActionTime"), c('TIME_TR'))		
	df.trm[, EventType:=NULL]
	set( df.trm, NULL, 'IDTR', df.trm[, as.integer(IDTR)] )
	set( df.trm, NULL, 'IDREC', df.trm[, as.integer(IDREC)] )
	#	check that transmission happen at unique times;	this seems to be the case in the DSPS model
	stopifnot( nrow(df.trm)==df.trm[, length(unique(TIME_TR))] )	
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)	
	#	
	df.trm[, FromDeme.FromHost:=NULL]
	df.trm[, ToDeme.ToHost:=NULL]
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	#
	#	reduce to time frame of interest
	#
	tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( df.trm[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	optional: multiple baseline index cases
	if(0)
	{
		#	consider only transmissions from between yr.start and yr.end
		df.trm	<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.start',][, as.numeric(v)] ) & TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
		#
		#	all transmitters infected before yr.start are re-declared as from 'outside the study'. This is coded with negative POPIDs.
		#
		setkey(df.trm, TIME_TR)
		tmp		<- unique( subset( df.trm, IDTR_TIME_INFECTED<pipeline.args['yr.start',][, as.numeric(v)], IDTR) )	
		tmp[, IDTR.new:=-rev(seq_len(nrow(tmp)))]
		df.trm	<- merge(df.trm, tmp, by='IDTR', all.x=TRUE)
		setnames(tmp, 'IDTR', 'IDPOP')
		df.ind	<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
		tmp		<- df.ind[, which(!is.na(IDTR.new))]
		set(df.ind, tmp, 'IDPOP', df.ind[tmp, IDTR.new])
		tmp		<- df.trm[, which(!is.na(IDTR.new))]
		set(df.trm, tmp, 'IDTR', df.trm[tmp, IDTR.new]) 	
		df.trm[, IDTR.new:=NULL]
		df.ind[, IDTR.new:=NULL]
		#	delete all POPID that don t appear as IDTR or IDREC in df.trm
		tmp		<- melt(subset(df.trm,select=c(IDTR,IDREC)), measure.vars=c('IDTR','IDREC'), value.name='IDPOP')
		tmp		<- unique(subset(tmp, select=IDPOP))
		df.ind	<- merge(df.ind, tmp, by='IDPOP')
		#
		cat(paste('\nFound transmissions between',pipeline.args['yr.start',][, as.numeric(v)],'-',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
		cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))	
		cat(paste('\nTotal transmitters before',pipeline.args['yr.start',][, as.numeric(v)],', these are treated as index cases before study start, n=', subset(df.trm, IDTR<0)[, length(unique(IDTR))] ))		
	}
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for -1
	stopifnot( subset(df.trm, is.na(IDTR_TIME_INFECTED))[, unique(IDTR)]==-1 )
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(1, TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm)
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 08-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 2, includes simulation of imports)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v2.R
#' @export
prog.HPTN071.input.parser.v2<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmissions
	#	
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep='', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	stopifnot( df.trm[, !any(is.na(TR_ACUTE))] )
	set(df.trm, df.trm[, which(TR_ACUTE<0)], 'TR_ACUTE', NA_integer_)
	set(df.trm, NULL, 'TR_ACUTE', df.trm[, factor(TR_ACUTE, levels=c(0,1), labels=c('No','Yes'))])
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	tmp		<- df.trm[, range(TIME_TR)]
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(TIME_TR>tmp[1] & length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(ceiling(tmp[2]), z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')	
	#	set baseline cases as negative ID
	tmp		<- df.trm[, which(IDTR=='-1')]
	cat(paste('\nFound index cases, n=', length(tmp)))
	set(df.trm, tmp, 'IDTR', rev(-seq_along(tmp)))
	#	check that all transmission times except baseline are unique
	tmp		<- subset(df.trm, TIME_TR>min(TIME_TR))
	stopifnot( nrow(tmp)==tmp[,length(unique(TIME_TR))] )	
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	cat(paste('\nTotal index cases, n=', df.trm[, length(which(unique(IDTR)<0))]))
	#
	#	prepare patient metavariables
	#
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c("Id","Gender","DoB","DateOfDeath","RiskGroup","T1_CD4350","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','T1_CD4350','CIRCM'))
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(T1_CD4350==-10)],'T1_CD4350',NA_real_)
	set(df.ind, df.ind[, which(T1_CD4350%in%c(-8,-9,-11,-12))],'T1_CD4350',Inf)	
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	stopifnot( df.ind[, !any(is.na(DOB))] )
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)		
	#	simulate time individual ready for sequencing
	#	ignore T1_CD4350 for now
	df.ind[, T1_CD4350:=NULL]
	df.ind[, T1_SEQ:= df.ind[, rgamma(nrow(df.ind),shape=9,scale=0.25 ) + TIME_TR]]
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( subset(df.trm, IDTR>0)[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for baseline cases
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 26-01-2015
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 4, uses date of diagnosis)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v4.R
#' @export
prog.HPTN071.input.parser.v4<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	#indir			<- ifelse(indir=='','/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Jan15/regional/150125',indir)
	#outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	indir			<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/freeze_Jan15/regional/150125'
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp150126'
	infile.ind		<- '150107_RUN001_IND.csv'
	infile.trm		<- '150107_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '150107_RUN001_IND.csv'
	outfile.trm		<- '150107_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmissions
	#	
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep=',', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	stopifnot( df.trm[, !any(is.na(TR_ACUTE))] )
	set(df.trm, df.trm[, which(TR_ACUTE<0)], 'TR_ACUTE', NA_integer_)
	set(df.trm, NULL, 'TR_ACUTE', df.trm[, factor(TR_ACUTE, levels=c(0,1), labels=c('No','Yes'))])
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	tmp		<- df.trm[, range(TIME_TR)]
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(ceiling(tmp[2]), z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')	
	#	stop if not all transmission times are unique, except at the end which does not matter
	stopifnot( subset(df.trm, TIME_TR<max(TIME_TR))[, length(unique(TIME_TR))==length(TIME_TR)] )
	#	set baseline cases as negative ID
	tmp		<- df.trm[, which(IDTR=='-1')]
	cat(paste('\nFound index cases, n=', length(tmp)))
	set(df.trm, tmp, 'IDTR', rev(-seq_along(tmp)))
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	cat(paste('\nTotal index cases, n=', df.trm[, length(which(unique(IDTR)<0))]))
	#
	#	prepare patient metavariables
	#
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c(	"Id","Sex","DoB","DoD","RiskGp"), c('IDPOP','GENDER','DOB','DOD','RISK'))
	setnames(df.ind, 	c(	"HIV_pos", "t_diagnosed","cd4_diagnosis","cd4atfirstART","t_1stARTstart","t_1stVLsupp_start","t_1stVLsupp_stop"), 
						c( 'HIV', 'DIAG_T','DIAG_CD4','ART1_CD4','ART1_T',"VLS1_TS","VLS1_TE"))
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])
	stopifnot( df.ind[, !any(is.na(DOB))] )
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, NULL, 'HIV', df.ind[, factor(HIV, levels=c(0,1), labels=c('N','Y'))])
	cat(paste('\nFound HIV+, n=', df.ind[, length(which(HIV=='Y'))]))
	set(df.ind, df.ind[, which(DIAG_T=='ND')], 'DIAG_T', NA_character_)
	set(df.ind, NULL, 'DIAG_T', df.ind[, as.numeric(DIAG_T)])
	set(df.ind, df.ind[, which(DIAG_CD4<0)], 'DIAG_CD4', NA_real_)
	stopifnot( df.ind[, !any(HIV=='N' & !is.na(DIAG_T)) ] )
	stopifnot( df.ind[, !any(!is.na(DIAG_T) & is.na(DIAG_CD4))] )
	cat(paste('\nFound not % undiagnosed , diagnosed=', paste( subset(df.ind, HIV=='Y')[, round(table(!is.na(DIAG_T)) / length(DIAG_T), d=2)], collapse=', ' )))	
	set(df.ind, df.ind[, which(ART1_CD4<0)], 'ART1_CD4', NA_real_)
	set(df.ind, df.ind[, which(ART1_T<0)], 'ART1_T', NA_real_)
	stopifnot( df.ind[, all(DIAG_T<ART1_T, na.rm=TRUE)] )
	cat(paste('\nFound not % on ART , not on ART among diagnosed=', paste(subset(df.ind, !is.na(DIAG_T))[, round( table(is.na(ART1_T))/length(ART1_T), d=2)], collapse=', ') ))
	set(df.ind, df.ind[, which(VLS1_TS<0)], 'VLS1_TS', NA_real_)
	set(df.ind, df.ind[, which(VLS1_TE<0)], 'VLS1_TE', NA_real_)
	tmp			<- df.ind[, which(ART1_T>=VLS1_TS)]
	cat(paste('\nFound ART1_T<VLS1_TS, setting VLS1_TS and VLS1_TE to NA, n=', length(tmp)))
	set(df.ind, tmp, 'VLS1_TS', NA_real_)
	set(df.ind, tmp, 'VLS1_TE', NA_real_)
	stopifnot( df.ind[, all(ART1_T<VLS1_TS, na.rm=TRUE)] )
	stopifnot( df.ind[, all(VLS1_TS<VLS1_TE, na.rm=TRUE)] )
	cat(paste('\nFound not % reached viral suppression , did not reach viral suppression among treated=',  paste(subset(df.ind, !is.na(ART1_T))[, round( table(is.na(VLS1_TS))/length(VLS1_TS), d=2)], collapse=', ') ))
	tmp			<- df.ind[, which(DIAG_CD4<ART1_CD4-DIAG_CD4*0.5 & DIAG_CD4>250)]
	cat(paste('\nFound individuals whose CD4 at ART start is much higher than at diagnosis, n=', length(tmp)))
	#stopifnot( df.ind[, all(DIAG_CD4>ART1_CD4, na.rm=TRUE)] )
	stopifnot( df.ind[, all(DIAG_CD4>0, na.rm=TRUE)] )
	stopifnot( df.ind[, all(ART1_CD4>0, na.rm=TRUE)] )	
	#	add transmission time
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))	
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	stopifnot( df.ind[, all(TIME_TR<DIAG_T, na.rm=TRUE)] )
	stopifnot( df.ind[, !any(is.na(TIME_TR) & !is.na(DIAG_T))] )
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	reset times if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DIAG_T<TIME_TR)]
	set(df.ind, tmp, 'DIAG_T', df.ind[tmp, TIME_TR+(TIME_TR-DIAG_T)])
	tmp			<- df.ind[, which(ART1_T<DIAG_T)]
	set(df.ind, tmp, 'ART1_T', df.ind[tmp, DIAG_T+(DIAG_T-ART1_T)])
	stopifnot( df.ind[, all(ART1_T<VLS1_TS, na.rm=TRUE)] )	
	tmp			<- df.ind[, which(DOD<TIME_TR)]
	set(df.ind, tmp, 'DOD', df.ind[tmp, TIME_TR+(TIME_TR-DOD)])
	#	add if transmission in the last 6 mo	
	df.ind[, RECENT_TR:=df.ind[, factor(as.numeric((DIAG_T-TIME_TR)<.5), levels=c(0,1), labels=c('N','Y'))]]	
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)
	stopifnot( df.trm[, !any(TIME_TR<=IDTR_TIME_INFECTED, na.rm=TRUE)] )
	#	simulate time individual ready for sequencing
	df.ind	<- PANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2(df.ind, seqtime.mode= pipeline.args['seqtime.mode',][,v])	
	cat(paste('\nFound % sampled at or after ART start=', subset(df.ind, !is.na(DIAG_T))[, mean(!is.na(ART1_T) & T1_SEQ>=ART1_T, na.rm=TRUE)] ))
	cat(paste('\nFound % sampled after end of first viral suppression=', subset(df.ind, !is.na(DIAG_T))[, mean(!is.na(VLS1_TE) & T1_SEQ>=VLS1_TE, na.rm=TRUE)] ))
	# 
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters in sampling frame, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( subset(df.trm, IDTR>0)[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for baseline cases
	#tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	#set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase.v2(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler.v4(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}
##--------------------------------------------------------------------------------------------------------
##	olli originally written 23-10-2014
##--------------------------------------------------------------------------------------------------------
#' @title HPTN071 parser (version 3, uses CD4 counts)
#' @description Reads files from the epi simulator in directory \code{indir} and writes csv files
#' in directory \code{outdir} for the virus tree simulator. The program samples sequences according to
#' an exponentially increasing sampling fraction in the same way as \code{prog.HPTN071.input.parser.v1}.
#' In addition, transmissions are broken and treated as imported from outside the simulated population.
#' The infected of a broken transmission chain is considered a new index case of a transmission chain within the 
#' simulated population. All input arguments are specified via the \code{argv} 
#' string, see the Examples.
#' @return NULL. Saves simulations to file.
#' @example example/ex.seq.sampler.v2.R
#' @export
prog.HPTN071.input.parser.v3<- function()	
{
	require(data.table)
	verbose			<- 1
	with.plot		<- 1	
	pipeline.args	<- NULL
	indir			<- system.file(package="rPANGEAHIVsim", "misc")
	indir			<- ifelse(indir=='','/Users/Oliver/git/HPTN071sim/raw_trchain',indir)
	outdir			<- '/Users/Oliver/git/HPTN071sim/tmp140908'
	infile.ind		<- '140716_RUN001_IND.csv'
	infile.trm		<- '140716_RUN001_TRM.csv'
	infile.args		<- NA
	outfile.ind		<- '140716_RUN001_IND.csv'
	outfile.trm		<- '140716_RUN001_TRM.csv'
	#
	if(exists("argv"))
	{
		#	args input
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,6),
									indir= return(substr(arg,8,nchar(arg))),NA)	}))
		if(length(tmp)>0) indir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.ind= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,11),
									infile.trm= return(substr(arg,13,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.trm<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									infile.args= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) infile.args<- tmp[1]
		#	args output
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,7),
									outdir= return(substr(arg,9,nchar(arg))),NA)	}))
		if(length(tmp)>0) outdir<- tmp[1]		
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.ind= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.ind<- tmp[1]
		tmp<- na.omit(sapply(argv,function(arg)
						{	switch(substr(arg,2,12),
									outfile.trm= return(substr(arg,14,nchar(arg))),NA)	}))
		if(length(tmp)>0) outfile.trm<- tmp[1]		
	}
	if(verbose)
	{
		cat('\ninput args\n',paste(indir, infile.ind, infile.trm, outdir, outfile.ind, outfile.trm, sep='\n'))
	}	
	if(!is.na(infile.args))
	{
		load(infile.args)	#expect 'pipeline.args'
	}
	if(is.null(pipeline.args))
	{
		cat('\nCould not find pipeline.args, generating default')
		pipeline.args	<- rPANGEAHIVsim.pipeline.args()
	}
	stopifnot( all( c('yr.start', 'yr.end', 's.seed', 's.PREV.min', 's.PREV.max', 'epi.dt', 'epi.import')%in%pipeline.args[, stat] ) )
	#
	infile.ind	<- paste(indir, '/', infile.ind, sep='')
	infile.trm	<- paste(indir, '/', infile.trm, sep='')
	outfile.ind	<- paste(outdir, '/', outfile.ind, sep='')
	outfile.trm	<- paste(outdir, '/', outfile.trm, sep='')
	
	#	set seed
	set.seed( pipeline.args['s.seed',][, as.numeric(v)] )
	#
	#	prepare transmissions
	#	
	df.trm	<- as.data.table(read.csv(infile.trm, stringsAsFactors=FALSE, sep='', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	stopifnot( df.trm[, !any(is.na(TR_ACUTE))] )
	set(df.trm, df.trm[, which(TR_ACUTE<0)], 'TR_ACUTE', NA_integer_)
	set(df.trm, NULL, 'TR_ACUTE', df.trm[, factor(TR_ACUTE, levels=c(0,1), labels=c('No','Yes'))])
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	tmp		<- df.trm[, range(TIME_TR)]
	df.trm	<- df.trm[, {
				z<- TIME_TR
				if(length(IDTR)>1)
					z<- sort(runif(length(IDTR), z, min(ceiling(tmp[2]), z+pipeline.args['epi.dt',][, as.numeric(v)])))
				list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
			}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')	
	#	stop if not all transmission times are unique, except at the end which does not matter
	stopifnot( subset(df.trm, TIME_TR<max(TIME_TR))[, length(unique(TIME_TR))==length(TIME_TR)] )
	#	set baseline cases as negative ID
	tmp		<- df.trm[, which(IDTR=='-1')]
	cat(paste('\nFound index cases, n=', length(tmp)))
	set(df.trm, tmp, 'IDTR', rev(-seq_along(tmp)))
	cat(paste('\nFound transmissions, n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))
	cat(paste('\nTotal index cases, n=', df.trm[, length(which(unique(IDTR)<0))]))
	#
	#	prepare patient metavariables
	#
	df.ind	<- as.data.table(read.csv(infile.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c(	"Id","Gender","DoB","DateOfDeath","RiskGroup","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','CIRCM'))
	setnames(df.ind, c(	"t.cd4below350","t.cd4.500","t.cd4.350","t.cd4.200","t.aidsdeath","SPVL.category","t.sc"), c('T1_CD4_l350','T1_CD4_500','T1_CD4_350','T1_CD4_200',"DOAD","SPVL","DOSC"))		
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(DOAD==-10)], 'DOAD', ceiling(max(df.trm$TIME_TR))+1.)
	set(df.ind, df.ind[, which(DOSC==-10)], 'DOSC', NA_real_)
	set(df.ind, df.ind[, which(SPVL==-1)], 'SPVL', NA_integer_)
	set(df.ind, NULL, 'SPVL', df.ind[, factor(SPVL, levels=c(0, 1, 2, 3), labels=c('le40','le45','le50','g50'))])	
	set(df.ind, df.ind[, which(T1_CD4_500%in%c(-1,-10))],'T1_CD4_500',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_500<0, na.rm=TRUE)])
	set(df.ind, df.ind[, which(T1_CD4_350==-10)],'T1_CD4_350',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_350<0, na.rm=TRUE)] )
	set(df.ind, df.ind[, which(T1_CD4_200==-10)],'T1_CD4_200',NA_real_)
	stopifnot( df.ind[, !any(T1_CD4_200<0, na.rm=TRUE)] )
	set(df.ind, df.ind[, which(T1_CD4_l350==-10)],'T1_CD4_l350',NA_real_)
	set(df.ind, df.ind[, which(T1_CD4_l350%in%c(-8,-9,-11,-12))],'T1_CD4_l350',Inf)	
	stopifnot( df.ind[, !any(T1_CD4_l350<0, na.rm=TRUE)] )
	stopifnot( df.ind[, !any(DOD>DOAD+0.125)] )
	set(df.ind, NULL, c('T1_CD4_l350','DOAD'), NULL)
	#	add transmission time
	tmp			<- subset(df.trm, select=c(IDREC, TIME_TR))
	stopifnot( df.ind[, !any(is.na(DOB))] )
	setnames(tmp, 'IDREC','IDPOP')
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	cat(paste('\nFound individuals with a valid record, n=', nrow(df.ind)))
	cat(paste('\nFound individuals with an infection event, n=', nrow(subset(df.ind,!is.na(TIME_TR)))))
	#	reset DOSC if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DOSC<TIME_TR)]
	set(df.ind, tmp, 'DOSC', df.ind[tmp, TIME_TR+(TIME_TR-DOSC)])
	tmp			<- df.ind[, which(T1_CD4_500<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_500)])
	tmp			<- df.ind[, which(T1_CD4_350<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_350)])
	tmp			<- df.ind[, which(T1_CD4_200<TIME_TR)]
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, TIME_TR+(TIME_TR-T1_CD4_200)])
	#	reset DOD if needed, because TIME_TR got randomized by a small bit above
	tmp			<- df.ind[, which(DOD<TIME_TR)]
	set(df.ind, tmp, 'DOD', df.ind[tmp, TIME_TR+(TIME_TR-DOD)])
	#	set T1_CD4 if starting with lower count than 500
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_200) & is.na(T1_CD4_350))]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, TIME_TR])
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR-0.01])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_350) & is.na(T1_CD4_500))]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, TIME_TR])
	#	set T1_CD4 if ending with lower count than 200
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_350) & is.na(T1_CD4_200))]
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & !is.na(T1_CD4_500) & is.na(T1_CD4_350))]
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, DOD])
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD+0.01])
	tmp			<- df.ind[, which(!is.na(TIME_TR) & is.na(T1_CD4_500))]
	set(df.ind, tmp, 'T1_CD4_500', df.ind[tmp, DOD])
	set(df.ind, tmp, 'T1_CD4_350', df.ind[tmp, DOD+0.01])
	set(df.ind, tmp, 'T1_CD4_200', df.ind[tmp, DOD+0.02])
	#	T1_CD4 is only after acute phase, and before we interpolate from 1500. 
	#	add time of infection of transmitter to df.trm	
	tmp		<- subset(df.ind, select=c(IDPOP, TIME_TR))
	setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
	setkey(tmp, IDTR)
	df.trm	<- merge(df.trm, unique(tmp), by='IDTR', all.x=TRUE)
	stopifnot( df.trm[, !any(TIME_TR<=IDTR_TIME_INFECTED, na.rm=TRUE)] )
	#	simulate time individual ready for sequencing
	df.ind	<- PANGEA.Seqsampler.SimulateGuideToSamplingTimes(df.ind, seqtime.mode= pipeline.args['seqtime.mode',][,v])
	#
	#
	#	reduce to time frame of interest
	#
	#tmp		<- subset( df.trm, TIME_TR>=as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )[, IDREC]
	df.trm	<- subset( df.trm, TIME_TR<as.numeric( pipeline.args['yr.end',][, as.numeric(v)] ) )
	#df.ind	<- subset(df.ind, !IDPOP%in%tmp)
	df.ind	<- subset(df.ind, is.na(DOB) | DOB<pipeline.args['yr.end',][, as.numeric(v)] )
	df.ind	<- subset(df.ind, is.na(DOD) | DOD >= floor(min(df.trm$TIME_TR)) )
	cat(paste('\nFound individuals born before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.ind)))
	cat(paste('\nFound transmissions before',pipeline.args['yr.end',][, as.numeric(v)],', n=', nrow(df.trm)))
	cat(paste('\nTotal transmitters, n=', df.trm[, length(unique(IDTR))]))		
	stopifnot( length(setdiff( subset(df.trm, IDTR>0)[, IDTR], df.ind[, IDPOP] ))==0 )
	stopifnot( length(setdiff( df.trm[, IDREC], df.ind[, IDPOP] ))==0 )
	#	simulate a fraction of transmissions to be imports
	tmp		<- PANGEA.ImportSimulator.SimulateIndexCase(df.ind, df.trm, epi.import= pipeline.args['epi.import',][,as.numeric(v)])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind
	#
	#	set infection times for index case
	#
	#	add IDTR_TIME_INFECTED for baseline cases
	tmp		<- df.trm[, which(is.na(IDTR_TIME_INFECTED))]	
	set( df.trm, tmp, 'IDTR_TIME_INFECTED', df.trm[tmp, runif(length(tmp), TIME_TR-5, TIME_TR)] )	
	tmp		<- PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase(df.ind, df.trm, index.starttime.mode= pipeline.args['index.starttime.mode',][,v])
	df.trm	<- tmp$df.trm
	df.ind	<- tmp$df.ind	
	#
	#	sample sequences 
	#
	PANGEA.Seqsampler.v3(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=with.plot)
	#
	return(1)
}