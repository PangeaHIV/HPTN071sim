


dev.haircut<- function()	
{
	if(0)
	{
		#	read just one file
		#	determine statistics for each contig after LTR
		
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/InterestingContigAlignments'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut'
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',FILE))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]
		set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
		par		<- c('FRQx.quantile'=0.05, 'CNS_AGR.window'=200, 'GPS.window'=200)
		
		
		file	<- paste(indir, infiles[1, FILE], sep='/')
		#	read Contigs+Rrefs: cr
		cr		<- read.dna(file, format='fasta')			
		#	determine start of non-LTR position and cut 
		tmp		<- haircut.find.nonLTRstart(cr)
		cr		<- cr[, seq.int(tmp, ncol(cr))]
		#	determine reference sequences. 
		#	non-refs have the first part of the file name in their contig name and are at the top of the alignment
		tmp		<- strsplit(basename(file), '_')[[1]][1]
		tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
		stopifnot( all( tx[, which(CONTIG==1)] == seq.int(1, tx[, length(which(CONTIG==1))]) ) )
		cat(paste('\nFound contigs, n=', tx[, length(which(CONTIG==1))]))
		#	determine base frequencies at each site amongst references.
		tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
		cnsr	<- haircut.getconsensus(tmp, par, bases=c('a','c','g','t','-') )	#	CoNSensus of References: cnsr
		#	for each contig, determine %agreement with consensus on rolling window
		cnsc	<- rbind(cnsr, cr[subset(tx, CONTIG==1)[, TAXON],])
		#	determine first and last non-gap sites
		tx		<- data.table(	TAXON= rownames(cnsc), 
				FIRST= apply( as.character(cnsc), 1, function(x) which(x!='-')[1] ),
				LAST= ncol(cnsc)-apply( as.character(cnsc), 1, function(x) which(rev(x)!='-')[1] )		)
		#	get cut statistics
		cnsc.df	<- haircut.get.cut.statistics(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
		cnsc.df	<- haircut.get.cut.statistics(cnsc, tx, par, outdir=outdir, file=file, mode='rolling')		
	}
	if(0)
	{
		#	process several files
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/InterestingContigAlignments'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/interesting_150408'		
		par		<- c('FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir)
	}
	if(0)
	{
		#	run mafft --add to get contigs+ref
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/contigs_150408'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		reffile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta'

		haircutwrap.align.contigs.with.ref(indir, reffile, outdir)
	}
}
##--------------------------------------------------------------------------------------------------------
##	HAIRCUT program, version 15086 to: 
##	- align contigs to references
##	- calculate and save haircut statistics
##--------------------------------------------------------------------------------------------------------
prog.haircut.150806<- function()
{
	DATA	<<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut'
	#DATA	<<- '/work/or105/Gates_2014/2015_PANGEA_haircut'
	if(0)
	{		
		indir	<- paste(DATA, 'contigs_150408', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste('/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut', 'contigs_150408_wref', sep='/' )
		reffile	<- paste(DATA, 'HIV1_COM_2012_genome_DNA_WithExtraA1UG.fasta', sep='/' )	
		haircutwrap.align.contigs.with.ref(indir, reffile, outdir)	
	}
	if(1)
	{
		indir	<- paste(DATA, 'contigs_150408_wref', sep='/' )
		outdir	<- paste(DATA, 'contigs_150408_wref_cutstat', sep='/' )		
		par		<- c('FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircutwrap.get.cut.statistics(indir, par, outdir=outdir)
	}
	
}
##--------------------------------------------------------------------------------------------------------
##	process all files in indir with 'haircut.align.contigs.with.ref'
##--------------------------------------------------------------------------------------------------------
haircutwrap.align.contigs.with.ref<- function(indir, reffile, outdir)
{
	infiles	<- data.table(INFILE=list.files(indir, pattern='fasta$',recursive=T))
	infiles	<- subset(infiles, !grepl('Curated', INFILE))
	infiles[, OUTFILE:= gsub('hiv\\.fasta','raw_wRefs\\.fasta', gsub('hiv_cut','cut_wRefs',INFILE))]
	set(infiles, NULL, 'OUTFILE', infiles[, basename(OUTFILE)])
	invisible(infiles[1:20,][, {
						haircut.align.contigs.with.ref(paste(indir,'/',INFILE,sep=''), reffile, paste(outdir,'/',OUTFILE,sep=''))
					}, by='INFILE'])
	NULL
}
##--------------------------------------------------------------------------------------------------------
##	call to MAFFT to align contigs with reference compendium
##--------------------------------------------------------------------------------------------------------
haircut.align.contigs.with.ref<- function(infile, reffile, outfile)
{
	#mafft --reorder --anysymbol --add new_sequences --auto input
	tmp		<- c( 	gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',infile,fixed=T),fixed=T),fixed=T),
					gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',reffile,fixed=T),fixed=T),fixed=T),
					gsub(' ','\\ ',gsub('(','\\(',gsub(')','\\)',outfile,fixed=T),fixed=T),fixed=T)
					)
	cmd		<- paste('mafft --anysymbol --add ',tmp[1],' --auto ',tmp[2],' > ',tmp[3], sep='')
	cat('\ncalling')
	cat(cmd)
	system(cmd, ignore.stdout=TRUE, ignore.stderr=TRUE)	
	NULL
}
##--------------------------------------------------------------------------------------------------------
##	process all files in indir with 'haircut.get.cut.statistics'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.cut.statistics<- function(indir, par, outdir=indir)	
{
	require(zoo)
	#	read just one file
	#	determine statistics for each contig after LTR
	infiles		<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',FILE))]
	infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	#	check which contigs not yet processed
	infiles		<- merge(infiles, infiles[, {
				file	<- paste(indir, FILE, sep='/')
				tmp		<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
				options(show.error.messages = FALSE)		
				readAttempt		<-try(suppressWarnings(load(tmp)))
				list(	DONE=!inherits(readAttempt, "try-error")	)			
			}, by='FILE'], by='FILE')
	cat(paste('\nFound processed files, n=', infiles[, length(which(DONE))]))
	infiles		<- subset(infiles, !DONE)
	#
	#	infiles[, which(grepl('12559_1_5_cut',FILE))]	fls<- 41
	#	process files
	for(fls in infiles[, seq_along(FILE)])
	{
		#	see if not yet constructed
		file	<- paste(indir, infiles[fls, FILE], sep='/')
		cat(paste('\nProcess', file))
		tmp		<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
		options(show.error.messages = FALSE)		
		readAttempt		<-try(suppressWarnings(load(tmp)))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nFound results",tmp))			
		options(show.error.messages = TRUE)		
	}
	if(!resume || inherits(readAttempt, "try-error"))
		
		
		#	read Contigs+Rrefs: cr
		cr		<- read.dna(file, format='fasta')			
		#	determine start of non-LTR position and cut 
		tmp		<- haircut.find.nonLTRstart(cr)
		cat(paste('\nFound end of LTR at=', tmp-1))
		cr		<- cr[, seq.int(tmp, ncol(cr))]
		#	determine reference sequences. 
		#	non-refs have the first part of the file name in their contig name and are at the top of the alignment
		tmp		<- strsplit(basename(file), '_')[[1]][1]
		tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
		stopifnot( all( tx[, which(CONTIG==1)] == seq.int(1, tx[, length(which(CONTIG==1))]) ) )
		cat(paste('\nFound contigs, n=', tx[, length(which(CONTIG==1))]))
		#	determine base frequencies at each site amongst references.
		tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
		tmp		<- haircut.getconsensus(tmp, par, bases=c('a','c','g','t','-') )	#	CoNSensus of References: cnsr
		cnsr	<- tmp$DNAbin
		cnsr.df	<- tmp$DATATABLE
		#	for each contig, determine %agreement with consensus on rolling window
		cnsc	<- rbind(cnsr, cr[subset(tx, CONTIG==1)[, TAXON],])
		#	determine first and last non-gap sites
		tx		<- data.table(	TAXON= rownames(cnsc), 
				FIRST= apply( as.character(cnsc), 1, function(x) which(x!='-')[1] ),
				LAST= ncol(cnsc)-apply( as.character(cnsc), 1, function(x) which(rev(x)!='-')[1] )		)
		tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#	some contigs only map into LTR
		#	get cut statistics
		cnsc.df	<- haircut.get.cut.statistics(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
		#	get rolling CNS_FRQ
		cnsr.df[, CNS_FRQr:= rollapply( seq_len(nrow(cnsr.df)), width=par['CNS_FRQ.window'], FUN= function(z) mean(cnsr.df$CNS_FRQ[z]), align='center', partial=T )]		
		cnsc.df	<- merge(cnsc.df, subset(cnsr.df, select=c(SITE, CNS_FRQr)), by='SITE')
		cnsc.df[, PNG_ID:= infiles[fls, PNG_ID]]
		cnsc.df[, BLASTnCUT:= infiles[fls, BLASTnCUT]]
		cat(paste('\nSave contigs, n=', cnsc.df[, length(unique(TAXON))]))
		#	save
		file	<- paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(file)), sep='')
		cat(paste('\nSave to', file))
		save(cnsc, cnsc.df, file=file)		
	}	
	#	plot by PANGEA_ID
	invisible( infiles[,
			{
				cat(paste('\nPlot', PNG_ID))
				tmp		<- sapply(FILE, function(x)
						{
							paste(outdir, '/', gsub('\\.fasta',paste('_HAIRCUTSTAT_thr',100*par['FRQx.quantile'],'_aw',par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'],'_gw',par['GPS.window'],'.R',sep=''),basename(x)), sep='')
						}) 
				tmp		<- do.call('rbind',lapply(tmp, function(x)
								{
									load(x)
									cnsc.df
								}))
				tmp[, TAXONnCUT:= paste(TAXON,BLASTnCUT,sep='_BLASTnCUT:')]
				ggplot(tmp, aes(x=SITE, ymax=AGRpc, ymin=0, fill=BLASTnCUT, group=TAXONnCUT)) + 
						geom_ribbon(alpha=0.5) + 
						geom_line(aes(y=GPS), colour='grey30') +
						geom_line(aes(y=CNS_FRQr), colour='#FF7F00') +
						scale_x_continuous(breaks=seq(0,15e3, ifelse(tmp[,max(SITE)]>5e2, 5e2, floor(tmp[,max(SITE)/3])))) + 
						facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
						labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='black line: % gaps\norange line: % consensus call amongst references\nfill: % agreement with consensus')
				ggsave(w=10, h=2*cnsc.df[, length(unique(TAXON))], file= paste(outdir, '/', PNG_ID, '_CNSAGR_aw', par['CNS_AGR.window'],'_fw',par['CNS_FRQ.window'], '_gw',par['GPS.window'],'.pdf',sep=''))
				NULL
			}, by='PNG_ID'] )
}
##--------------------------------------------------------------------------------------------------------
##	calculate cut statistics for each contig
##--------------------------------------------------------------------------------------------------------
haircut.get.cut.statistics<- function(cnsc, tx, par, outdir=NA, file=NA, mode='rolling')
{
	require(zoo)
	stopifnot(mode%in%c('rolling','overall'))
	#	count overall disagreement with consensus		
	if(mode=='overall')
	{
		cnsc.df	<- merge(tx, subset(tx, TAXON!='consensus')[, 
						{
							tmp		<- cnsc[ c('consensus',TAXON), seq.int(FIRST,LAST)]
							
							list( DSGR= dist.dna(tmp, model='indel' ), GPS=mean( as.character( cnsc[ TAXON, seq.int(FIRST,LAST)] )=='-' ), LEN=LAST-FIRST+1 )
						},by='TAXON'], by='TAXON')	
	}
	#	count rolling agreement with consensus
	if(mode=='rolling')
	{
		cnsc.df	<- merge(tx, subset(tx, TAXON!='consensus')[, 
						{
							tmp		<- cnsc[ c('consensus',TAXON), seq.int(FIRST,LAST)]
							agrpc	<- 1-rollapply( seq_len(LAST-FIRST+1), width=par['CNS_AGR.window'], FUN= function(z) dist.dna(tmp[,z], model='indel' )/length(z), align='center', partial=T )
							tmp		<- as.character( cnsc[ TAXON, seq.int(FIRST,LAST)] )=='-'
							gps		<- rollapply( seq_len(LAST-FIRST+1), width=par['GPS.window'], FUN= function(z) mean(tmp[z]), align='center', partial=T )
							list( SITE=seq.int(FIRST,LAST),  AGRpc=agrpc, GPS=gps   )
						},by='TAXON'], by='TAXON')		
		if(!is.na(file) & !is.na(outdir))
		{
			ggplot(cnsc.df, aes(x=SITE, ymax=AGRpc, ymin=0, y=GPS, fill=TAXON, group=TAXON)) + 
					geom_ribbon(alpha=0.5) + geom_line(colour='grey30') +
					scale_x_continuous(breaks=seq(0,15e3, 5e2)) + 
					facet_grid(TAXON~.) + theme_bw() + theme(strip.text=element_blank(), legend.position='bottom') +
					labs(fill='Contig', x='position on consensus w/o LTR', y='line: % gaps\nfill: % agreement with consensus')
			ggsave(w=10, h=2*cnsc.df[, length(unique(TAXON))], file= paste(outdir, '/', gsub('\\.fasta',paste('_CNSAGR_aw',par['CNS_AGR.window'],'_gw',par['GPS.window'],'.pdf',sep=''),basename(file)), sep=''))	
		}
	}
	cnsc.df
}
##--------------------------------------------------------------------------------------------------------
##	find the first site outside the LTR in the contig + reference alignment 
##--------------------------------------------------------------------------------------------------------
haircut.find.nonLTRstart<- function(cr)
{
	ans				<- seq.find.pos.of.pattern(cr, pattern='a-*t-*g-*g-*g-*t-*g-*c-*g-*a-*g-*a-*g-*c-*g-*t-*c-*a', row.name='B.FR.83.HXB2_LA')
	stopifnot(length(ans)==1, ans>0)
	ans
}
##--------------------------------------------------------------------------------------------------------
##	determine the consensus sequence in the reference alignment
##	variable sites may have an NA consensus base call, depending on the consensus quantile cutoff 'FRQx.quantile'
##--------------------------------------------------------------------------------------------------------
haircut.getconsensus	<- function(seq, par, bases=c('a','c','g','t','-') )
{
	stopifnot('FRQx.quantile'%in%names(par))
	#	calculate base frequency Profile among References: rp
	rp		<- sapply(seq_len(ncol(seq)), function(i) base.freq(seq[,i], freq=F, all=T))
	rp		<- as.data.table( t(rp[bases,]) )
	rp[, SITE:= seq_len(nrow(rp))]
	rp		<- melt(rp, id.vars='SITE', variable.name='BASE', value.name='FRQ')
	set(rp, NULL, 'BASE', rp[, factor(BASE, levels=bases, labels=bases)])
	#	renormalize bases and ignore all other calls amongst references	
	rp		<- rp[ , list(BASE=BASE, FRQ=FRQ/sum(FRQ), FRQxu= max(FRQ), FRQx= max(FRQ/sum(FRQ))), by='SITE']
	#	get consensus threshold as lower quantile of a one-inflated Beta distribution, par['FRQx.quantile']	
	setkey(rp, SITE)
	tmp		<- unique(rp)
	#	ggplot( tmp, aes(x=FRQx)) + geom_histogram()	#	looks like BEINF
	#set(rp, NULL, 'DUMMY', rp[, as.integer(ceiling(SITE/500))])
	#ggplot(subset(rp, DUMMY<=2 ), aes(x=SITE, y=FRQ, colour=BASE, group=BASE)) + geom_line() + facet_wrap(~DUMMY, ncol=1, scales='free')
	#ggplot(subset(rp, BASE=='a'), aes(x=SITE, y=FRQx, colour=BASE, group=BASE)) + geom_line() + facet_wrap(~DUMMY, ncol=1, scales='free') + theme_bw() + theme(strip.text = element_blank())
	if(is.na(par['FRQx.thr']))
	{
		qu.m	<- gamlss(FRQx~1, data=tmp, family=BEINF)
		qu.par	<- c('mu'=as.double(predict(qu.m, type='response', what='mu')[1]), 'sigma'=as.double(predict(qu.m, type='response', what='sigma')[1]), 'nu'=as.double(predict(qu.m, type='response', what='nu')[1]), 'tau'=as.double(predict(qu.m, type='response', what='tau')[1]) ) 
		qu.par	<- qBEINF( par['FRQx.quantile'], mu=qu.par['mu'], sigma=qu.par['sigma'], nu=qu.par['nu'], tau=qu.par['tau'])		
	}
	if(!is.na(par['FRQx.thr']))
		qu.par	<- par['FRQx.thr']
	cat(paste('\nMinimum base frequency to call a consensus base=', round(qu.par, d=3)))
	set(rp, rp[, which(FRQx<qu.par)], 'BASE', NA_character_)
	set(rp, NULL, 'BASE', rp[, factor(as.character(BASE), levels=bases, labels=bases)])
	#	get consensus base above lower quantile. Some consensus bases will be NA
	crp		<- rp[, list(CNS_BASE= BASE[ which.max(FRQ) ], CNS_FRQ=FRQxu[1]), by='SITE']
	set(crp, crp[, which(is.na(CNS_BASE))], 'CNS_BASE','?')	
	tmp		<- crp[, as.DNAbin(t(as.matrix(CNS_BASE)))]
	rownames(tmp)	<- 'consensus'
	list(DNAbin=tmp, DATATABLE=crp)
}
##--------------------------------------------------------------------------------------------------------
##	find the first site with pattern in alignment 
##	TODO: put into hivclust
##--------------------------------------------------------------------------------------------------------
seq.find.pos.of.pattern<- function(seq, pattern, row.name)
{
	stopifnot(is.matrix(seq), !is.na(pattern), !is.na(row.name))
	tmp		<- which(grepl(row.name, rownames(seq), fixed=1))
	stopifnot(length(tmp)==1)	
	regexpr(pattern, paste(as.character( seq[tmp, ] ),collapse=''))	
}