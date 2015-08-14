


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
	if(0)
	{
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_cutsubraw.R'
		txe		<- haircutwrap.get.subset.among.raw.and.cut.contigs(indir, outfile)
	}
	if(0)	# get training data
	{	
		#	get contigs that are to be used for training: this determines the 1's and excludes some 0's
		#	read curated contigs 
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEA_data/contigs_150408/CuratedAlignmentsToRefs'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_curated.R'
		ctrain	<- haircut.get.curated.contigs(indir, outfile)		
		setnames(ctrain, c('FILE','LEN','FIRST','LAST'),c('ANS_FILE','ANS_LEN','ANS_FIRST','ANS_LAST'))
		#	remove cut & raw contigs that are identical as double 1's and as 0's
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_trainingset_subsets.R'
		ctrain	<- haircut.get.training.contigs(indir, outfile, ctrain)
		set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(ctrain, 'CUT', 'BLASTnCUT')		
		#	now create training data set: 
		#	expand the curated contigs in the training data to 1's per site, and
		#	expand all non-existing curated contigs to 0's, provided they are not excluded because identical to a 1
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref_cutstat'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- 'contigs_150408_train'
		par		<- c('FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200)
		haircut.get.training.data(indir, ctrain, par, outdir, outfile)
	}
	if(0)	#	call contigs on training data and plot
	{		
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train/model_150811a.R'
		#	get model coefficients across the chunks
		tmp						<- haircut.get.fitted.model.150811a(NULL, outfile)
		ctrmc					<- tmp$coef
		ctrev					<- tmp$ev
		model.150811a.predict	<- tmp$predict
		#	get info on which cut contig is a subset of the corresponding raw contig
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_cutsubraw.R'
		txe		<- haircutwrap.get.subset.among.raw.and.cut.contigs(indir, outfile)
		set(txe, NULL, 'CUT', txe[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(txe, 'CUT', 'BLASTnCUT')	
		#	get contigs that were used for training		
		outfile	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_trainingset_subsets.R'
		ctrain	<- haircut.get.training.contigs(NULL, outfile, NULL)
		set(ctrain, NULL, 'CUT', ctrain[, factor(CUT, levels=c('cut','raw'), labels=c('Y','N'))])
		setnames(ctrain, 'CUT', 'BLASTnCUT')		
		#	get covariates for all contigs
		indir.st<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref_cutstat'
		indir.al<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_wref'
		outdir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_model150811a'
		par		<- c(	'FRQx.quantile'=0.05, 'FRQx.thr'=0.566, 'CNS_FRQ.window'=100, 'CNS_AGR.window'=200, 'GPS.window'=200, 
						'PRCALL.thr'=0.8, 'PRCALL.cutprdcthair'=100, 'PRCALL.cutrawgrace'=100, 'PRCALL.rmintrnlgpsblw'=100 ,'PRCALL.rmintrnlgpsend'=9700)
		haircutwrap.get.call.for.PNG_ID(indir.st,indir.al,outdir,ctrmc,ctrev,predict.fun,txe,par,ctrain=ctrain)
	}
	if(0)	# evaluate training data
	{	
		#	
		indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_haircut/contigs_150408_train'
		outfile	<- 'model_150811a.R'
		#
		#	find 20 worst false pos
		#
		ctrfp	<- do.call('rbind', lapply(seq(1,10000,200), function(site)
				{
					ctr		<- haircut.calculate.training.data(indir, site)	
					tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
					ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
					do.call('rbind',lapply( ctr[, unique(CHUNK)], function(chunk){
								ctrch	<- subset(ctr, CHUNK==chunk)
								ctrchfp	<- subset( ctrch, ANS_CALL==0 & AGRpc>subset( ctrch, ANS_CALL==1 )[, min(AGRpc)])
								ctrchfp	<- ctrchfp[, list(ANS_CALL=mean(ANS_CALL), AGRpc=mean(AGRpc), GPS=mean(GPS), CNS_FRQr=mean(CNS_FRQr)), by=c('PNG_ID','TAXON')]
								setkey(ctrchfp, AGRpc)
								ctrchfp[ seq.int(max(1,nrow(ctrchfp)-20), nrow(ctrchfp)), ]	
							}))										
				}))		
		save(ctrfp, file= paste(indir,'/fp.R',sep=''))
		setkey(ctrfp, PNG_ID, TAXON)
		ctrfp[, length(unique(TAXON))]
		#
		#	show all three types of information across sites
		#
		invisible(do.call('rbind', lapply(seq(1,10000,200), function(site)
					{
						cat('\nProcess',site)
						ctr		<- haircut.calculate.training.data(indir, site)	
						tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
						ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
						ggplot(ctr, aes(x=GPS, y=AGRpc, colour=factor(ANS_CALL, levels=c(0,1), labels=c('N','Y')))) + 
								geom_hline(yintercept=ctrch$CNS_FRQr[1]) + geom_point(alpha=0.6, size=1) + 
								facet_wrap(~CHUNK, ncol=5) +
								theme_bw() + labs(x='%gappiness', y='%agreement with consensus\nline: %agreement of consensus with references', colour='In curated contigs') + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/',outfile,'_ANSCALLBYAGRpcGPS_SITE',site,'.pdf',sep=''), w=9, h=9)						
					})))
		#
		#	calculate model coefficients across the chunks
		#
		tmp						<- haircut.get.fitted.model.150811a(indir, outfile)
		ctrmc					<- tmp$coef
		ctrev					<- tmp$ev
		model.150811a.predict	<- tmp$predict
	
		
		
		
		#	select threshold
		ctrt		<- subset(ctrev, THR==0.8)
		ggplot(melt(ctrt, id.vars=c('CHUNK'), measure.vars=c('THR','SENS','SPEC','FDR','FOR')), aes(x=CHUNK, y=value, colour=variable)) + geom_line() +
					scale_x_continuous(breaks=seq(0,20e3,1e3), expand=c(0,0)) +
					scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
					theme_bw() + labs(x='base relative to HXB2') + theme(panel.grid.major=element_line(colour="grey50", size=0.2), panel.grid.minor=element_line(colour="grey70", size=0.1))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_THR_const80.pdf',sep=''), w=10, h=4)
		#
		setkey(ctrev, CHUNK, FDR)
		ctrt		<- ctrev[, {
					z		<- which(FDR<= max(0.01, min(FDR, na.rm=T)))[1]
					if(!is.finite(z) | THR[z]<0.8)
						z	<- which(THR==0.8)
					list(THR=THR[z], SENS=SENS[z], SPEC=SPEC[z], FDR=FDR[z], FOR=FOR[z])
				}, by='CHUNK']
		ggplot(melt(ctrt, id.vars=c('CHUNK'), measure.vars=c('THR','SENS','SPEC','FDR','FOR')), aes(x=CHUNK, y=value, colour=variable)) + geom_line() +
				scale_x_continuous(breaks=seq(0,20e3,1e3), expand=c(0,0)) +
				scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
				theme_bw() + labs(x='base relative to HXB2') + theme(panel.grid.major=element_line(colour="grey50", size=0.2), panel.grid.minor=element_line(colour="grey70", size=0.1))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_THR_BYFDR01.pdf',sep=''), w=10, h=4)
		#	just const threshold??



		tmp		<- ctrev[, list(FDR=min(FDR)), by='CHUNK']
		ggplot(tmp, aes(x=FDR)) + geom_histogram()
		subset(tmp, FDR>0.05)
		
		ctr		<- haircut.calculate.training.data(indir, 1)	
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
		ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
		
		ctrp	<- merge(ctr, ctrmc, by='CHUNK')
		ctrp[, PR_CALL:=model.150811a.predict(AGRpc, GPS, BETA0, BETA1, BETA2)]
		setkey(ctrp, CHUNK)
		ctrev	<- as.data.table(expand.grid(THR= seq(0.05,0.95,0.01), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
		ctrev	<- ctrev[, {
					tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=THR, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
					list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
				}, by=c('CHUNK','THR')]
		ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))
		#	plot Sensitivity & Specificity
		ggplot(melt(ctrev, id.vars=c('CHUNK','THR'), measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=100*value, colour=CHUNK)) + 
				geom_line() + labs(x='threshold on predicted call probability',y='%', colour='site\n(base relative to HXB2)') +
				theme(legend.position='bottom') + facet_wrap(~variable, ncol=2, scales='free') +
				guides(col = guide_legend(ncol=5, byrow=TRUE))
		ggsave(file=paste(outdir,'/',outfile,'_model.150811a_SensSpec_SITE',site,'.pdf',sep=''), w=9, h=9)
		#
		ctrev
		
		ctrp[CHUNK, table(ANS_CALL, PR_CALL>=THR)]
		

		#	calculate model coefficients across the chunks		
		ctrmc	<- do.call('rbind',lapply(ctr[, unique(CHUNK)], function(chunk)
				{
					ctrch	<- subset(ctr, CHUNK==chunk)
					ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
					ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3])					
				}))
		#	evaluate FN FP
		#	return the whole lot
		#	determine threshold based on FN FP, and make call
	
		
		
		#	add rolling mean by AGRpc
		site	<- 1610; chunk<- '1610'
		site	<- 9350; chunk<- '9350'
		site	<- 3650; chunk<- '3650'
		ctr		<- haircut.load.training.data(indir, site)		
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)],10)
		ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
		#tmp		<- seq(0,1.01,0.01)
		#ctr[, GPSC:= cut(GPS, breaks=tmp, labels=tmp[-length(tmp)])]
		#ctr[, AGRpcC:= cut(AGRpc, breaks=tmp, labels=tmp[-length(tmp)])]		
		ctrch	<- subset(ctr, CHUNK==chunk)
		ctrchs	<- ctrch[, {
					z	<- sample(length(AGRpc), min(length(AGRpc),50))
					list(AGRpc=AGRpc[z], GPS=GPS[z], ANS_CALL=ANS_CALL[z], CHUNK=CHUNK[z])	
				}, by=c('AGRpcC','GPSC')]
		
		
		
		ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrch, family=BI())
		ctrchms	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrchs, family=BI())
		ctrchm2	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrchs, family=BB())
		ctrchm3	<- gamlss(ANS_CALL~AGRpc+GPS, sigma.formula=~AGRpc+GPS, data=ctrchs, family=BB())
		#	predict sigma: how does it look?
		ctrchsp	<- as.data.table(expand.grid(AGRpc=seq(0.01,0.99,0.01), GPS=seq(0.1,0.9,0.1)))
		mu		<- predict(ctrchm2, newdata=ctrchsp, what='mu', type='response', se.fit=TRUE)
		
		predict(ctrchm2, what='mu', type='response', se.fit=TRUE)
		max(predict(ctrchm2, what='mu', type='response', se.fit=TRUE)$se.fit)
		
		ctrchsp[, mu:=mu]
		ctrchsp[, TYPE:='s']
		ctrchp	<- as.data.table(expand.grid(AGRpc=seq(0.01,0.99,0.01), GPS=seq(0.1,0.9,0.1)))
		mu		<- predict(ctrchm, newdata=ctrchp, what='mu', type='response')
		ctrchp[, mu:=mu]
		ctrchp[, TYPE:='a']
		ctrchp	<- rbind(ctrchp, ctrchsp)
		ggplot(ctrchp, aes(x=AGRpc, y=mu, group=TYPE, colour=TYPE)) + geom_line() + facet_grid(GPS~.)
		sigma	<- predict(ctrchm3, newdata=ctrchp, what='sigma', type='response')				
		ctrchp[, mu:=mu]
		ctrchp[, sigma:=sigma]
		ctrchp[, varn1:=mu*(1-mu)]
		
		
		setkey(ctrch, AGRpc)
		ctrch[, ANS_CALLrm:=ctrch[, rollapply(ANS_CALL, width=100, FUN=mean, align="center", partial=TRUE)]]	
		setkey(ctrch, GPSc, AGRpc)
		ctrch	<- merge(ctrch, ctrch[, list(TAXON=TAXON, SITE=SITE, BLASTnCUT=BLASTnCUT, ANS_CALLrm2=rollapply(ANS_CALL, width=100, FUN=mean, align="center", partial=TRUE)), by='GPSc'] , by=c('TAXON','SITE','BLASTnCUT','GPSc'))
		
		#ctrchm	<- gamlss(ANS_CALL~AGRpc+bs(GPSc,df=3), data=ctrch, family=BI())
		
		
		ctrch[, PR_CALL:= predict(ctrchm, type='response', what='mu')]
		setkey(ctrch, GPSc, AGRpc)
		ggplot(melt(ctrch, id.vars=c('AGRpc','GPSc'), measure.vars=c('ANS_CALLrm2','PR_CALL')), aes(x=AGRpc, y=value, colour=variable)) + geom_line() +facet_wrap(~GPSc, ncol=4)
		
		
		ctrchm	<- gamlss(ANS_CALL~AGRpc, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
		ctrchmc	<- data.table(CHUNK=chunk,BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2])
		#ctrch[, PR_CALL:= predict(ctrchm, type='response', what='mu')]
		ctrch	<- merge(ctrch, ctrchmc, by='CHUNK')
		ctrch[, PR_CALL2:= exp(BETA0+BETA1*AGRpc)/(exp(BETA0+BETA1*AGRpc)+1)]
		
		
		ggplot(melt(ctrch, id.vars='AGRpc', measure.vars=c('ANS_CALLrm','PR_CALL')), aes(x=AGRpc, y=value, colour=variable)) + geom_line()
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_ANS_CALLrmBYAGRpc.pdf',sep=''), w=6, h=4)
		
		
		plot(ctrchm)	#looks pretty good
		Rsq(ctrchm)		#60%
		AIC(ctrchm)
		
		
		ggplot(ctrch, aes(x=PR_CALL, y=ANS_CALL)) + geom_point()
		#	get sens and spec
		ctrev	<- data.table(THR= seq(0.05,0.95,0.01))
		ctrev	<- ctrev[, {
					tmp	<- ctrch[, table(ANS_CALL, PR_CALL>=THR)]
					list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
				}, by='THR']
		ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by='THR'], by='THR')
		ggplot(melt(ctrev, id.vars='THR', measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=value, colour=variable)) + geom_line() + theme(legend.position='bottom')
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_SensSpec.pdf',sep=''), w=4, h=4)
		ggplot(ctrev, aes(x=1-SPEC, y=SENS)) + geom_line() + scale_x_continuous(limit=c(0,1),expand=c(0,0)) + scale_y_continuous(limit=c(0,1),expand=c(0,0))
		ggsave(file=paste(outdir,'/',outfile,'_SITE',chunk,'_ROC.pdf',sep=''), w=4, h=4)
		#ctrchm2	<- gamlss(ANS_CALL~AGRpc+GPS+CNS_FRQr, data=ctrch, family=BB())	#BB does not always converge
	}
}
##--------------------------------------------------------------------------------------------------------
##	wrapper to call 'haircutwrap.get.call.for.PNG_ID'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.call.for.PNG_ID.150814<- function(indir.st,indir.al,outdir,ctrmc,ctrev,predict.fun,txe,par,ctrain=NULL)
{
	infiles	<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',INFILE))]
	infiles[, BLASTnCUT:= regmatches(INFILE,regexpr('cut|raw',INFILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	alfiles <- data.table(ALFILE=list.files(indir.al, pattern='\\.fasta$', recursive=T))
	alfiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',basename(ALFILE)))]
	alfiles[, BLASTnCUT:= regmatches(basename(ALFILE),regexpr('cut|raw',basename(ALFILE)))]
	set(alfiles, NULL, 'BLASTnCUT', alfiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	infiles	<- merge(infiles, alfiles, by=c('PNG_ID','BLASTnCUT'))
	
	#	predict by PANGEA_ID
	cnsc.info	<-  infiles[,
					{
						cat(paste('\nProcess', PNG_ID))
						if(0)	#devel
						{
							PNG_ID<- png_id	<- '15172_1_32'
							PNG_ID<- png_id	<- '12559_1_11'
							PNG_ID<- png_id	<- '14728_1_84'
							PNG_ID<- png_id	<- '14938_1_10'
							PNG_ID<- png_id	<- '14728_1_82'
							PNG_ID<- png_id	<- '12559_1_24'
							PNG_ID<- png_id	<- '12559_1_81'
							#PNG_ID<- png_id	<- '12559_1_87'
							files	<- subset(infiles, PNG_ID==png_id)[, INFILE]
							alfiles	<- subset(infiles, PNG_ID==png_id)[, ALFILE]
							bc		<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
							tmp		<- haircut.get.call.for.PNG_ID.150814(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)
						}
						#
						tmp		<- haircut.get.call.for.PNG_ID.150814(indir.str, indir.al, PNG_ID, INFILE, ALFILE, BLASTnCUT, par, ctrmc, predict.fun)
						crs		<- tmp$crs
						cnsc.df	<- tmp$cnsc.df	
						#	handle output
						if(any(grepl(PNG_ID,rownames(crs[['N']]))))
						{
							tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironraw.fasta',sep='')
							cat('\nWrite to file', tmp)
							write.dna(crs[['N']], file=tmp, format='fasta', colsep='', nbcol=-1)							
						}
						if(any(grepl(PNG_ID,rownames(crs[['Y']]))))
						{
							tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironcut.fasta',sep='')
							cat('\nWrite to file', tmp)
							write.dna(crs[['Y']], file=tmp, format='fasta', colsep='', nbcol=-1)							
						}
						#	see if there is curated contig available
						if(!is.null(ctrain))
						{
							cnsc.df	<- merge(cnsc.df, subset(ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, ANS_FIRST, ANS_LAST)), all.x=T, by=c('PNG_ID','TAXON','BLASTnCUT'))		
							cnsc.df[, CUR_CALL:=NA_integer_]
							set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & SITE>=ANS_FIRST & SITE<=ANS_LAST)], 'CUR_CALL', 1L)
							set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & (SITE<ANS_FIRST | SITE>ANS_LAST))], 'CUR_CALL', 0L)
							set(cnsc.df, cnsc.df[, which(is.na(CUR_CALL))], 'CUR_CALL', 0L)								
						}	
						if(is.null(ctrain))
							cnsc.df[, CUR_CALL:=NA_integer_]
						#	save as R
						tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.R',sep='')
						cat('\nSave to file', tmp)
						save(cnsc.df, crs, file=tmp)
						#	plot
						cnsc.df[, TAXONnCUT:= paste(TAXON,BLASTnCUT,sep='_BLASTnCUT:')]
						ggplot(cnsc.df, aes(x=SITE, fill=BLASTnCUT, group=TAXONnCUT)) +
								geom_ribbon(aes(ymax=CALL, ymin=0), alpha=0.5) +
								geom_line(aes(y=PR_CALL), colour='black') +
								geom_line(aes(y=CNS_FRQr), colour='blue') +
								geom_line(aes(y=CUR_CALL), colour='red') +
								scale_x_continuous(breaks=seq(0,15e3, ifelse(cnsc.df[,max(SITE)]>5e2, 5e2, floor(cnsc.df[,max(SITE)/3])))) + 
								facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
								labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='fill: predicted call\nblack line: predictive probability\nred line: curated call')
						tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.pdf',sep='')
						cat('\nPlot to file', tmp)
						ggsave(w=10, h=3*cnsc.df[, length(unique(TAXON))], file=tmp)
						#	report confidence score
						subset(cnsc.df, CALL==1)[, list(QUANTILE=c(0,0.01,0.05,0.1,0.2,0.5), PR_CALL=quantile(PR_CALL, p=c(0,0.01,0.05,0.1,0.2,0.5))), by=c('PNG_ID','TAXON','BLASTnCUT')]
					}, by='PNG_ID']
}	
##--------------------------------------------------------------------------------------------------------
##	wrapper to call 'haircutwrap.get.call.for.PNG_ID.150811'
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.call.for.PNG_ID.150811<- function(indir.st,indir.al,outdir,ctrmc,ctrev,predict.fun,par,ctrain=NULL)
{
	infiles	<- data.table(INFILE=list.files(indir.st, pattern='\\.R$', recursive=T))
	infiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',INFILE))]
	infiles[, BLASTnCUT:= regmatches(INFILE,regexpr('cut|raw',INFILE))]
	set(infiles, NULL, 'BLASTnCUT', infiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	alfiles <- data.table(ALFILE=list.files(indir.al, pattern='\\.fasta$', recursive=T))
	alfiles[, PNG_ID:= gsub('_wRefs.*','',gsub('_cut|_raw','',basename(ALFILE)))]
	alfiles[, BLASTnCUT:= regmatches(basename(ALFILE),regexpr('cut|raw',basename(ALFILE)))]
	set(alfiles, NULL, 'BLASTnCUT', alfiles[, factor(BLASTnCUT, levels=c('cut','raw'), labels=c('Y','N'))])
	infiles	<- merge(infiles, alfiles, by=c('PNG_ID','BLASTnCUT'))
	
	#	predict by PANGEA_ID
	invisible( infiles[,
					{
						cat(paste('\nProcess', PNG_ID))
						#PNG_ID<- png_id	<- '13548_1_21'
						#files	<- subset(infiles, PNG_ID==png_id)[, INFILE]
						#alfiles	<- subset(infiles, PNG_ID==png_id)[, ALFILE]
						#bc		<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]			
						#tmp		<- haircut.get.call.for.PNG_ID(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)
						tmp		<- haircut.get.call.for.PNG_ID(indir.str, indir.al, PNG_ID, INFILE, ALFILE, BLASTnCUT, par, ctrmc, predict.fun)
						crs		<- tmp$crs
						cnsc.df	<- tmp$cnsc.df	
						#	handle output
						tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironraw.fasta',sep='')
						cat('\nWrite to file', tmp)
						write.dna(crs[['N']], file=tmp, format='fasta', colsep='', nbcol=-1)
						tmp		<- paste(outdir,'/',PNG_ID,'_wref_nohaironcut.fasta',sep='')
						cat('\nWrite to file', tmp)
						write.dna(crs[['Y']], file=tmp, format='fasta', colsep='', nbcol=-1)
						#	see if there is curated contig available
						if(!is.null(ctrain))
						{
							cnsc.df	<- merge(cnsc.df, subset(ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, ANS_FIRST, ANS_LAST)), all.x=T, by=c('PNG_ID','TAXON','BLASTnCUT'))		
							cnsc.df[, CUR_CALL:=NA_integer_]
							set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & SITE>=ANS_FIRST & SITE<=ANS_LAST)], 'CUR_CALL', 1L)
							set(cnsc.df, cnsc.df[, which(!is.na(ANS_FIRST) & (SITE<ANS_FIRST | SITE>ANS_LAST))], 'CUR_CALL', 0L)
							set(cnsc.df, cnsc.df[, which(is.na(CUR_CALL))], 'CUR_CALL', 0L)								
						}	
						if(is.null(ctrain))
							cnsc.df[, CUR_CALL:=NA_integer_]
						#	save as R
						tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.R',sep='')
						cat('\nSave to file', tmp)
						save(cnsc.df, crs, file=tmp)
						#	plot
						cnsc.df[, TAXONnCUT:= paste(TAXON,BLASTnCUT,sep='_BLASTnCUT:')]
						ggplot(cnsc.df, aes(x=SITE, fill=BLASTnCUT, group=TAXONnCUT)) +
								geom_ribbon(aes(ymax=CALL, ymin=0), alpha=0.5) +
								geom_line(aes(y=PR_CALL), colour='black') +
								geom_line(aes(y=CUR_CALL), colour='red') +
								scale_x_continuous(breaks=seq(0,15e3, ifelse(cnsc.df[,max(SITE)]>5e2, 5e2, floor(cnsc.df[,max(SITE)/3])))) + 
								facet_wrap(~TAXONnCUT, ncol=1) + theme_bw() + theme(legend.position='bottom') +
								labs(fill='Contig BLASTnCUT', x='position on consensus w/o LTR', y='fill: predicted call\nblack line: predictive probability\nred line: curated call')
						tmp	<- paste(outdir, '/', PNG_ID, '_wref_nohaironcutraw.pdf',sep='')
						cat('\nPlot to file', tmp)
						ggsave(w=10, h=3*cnsc.df[, length(unique(TAXON))], file=tmp)
						NULL
					}, by='PNG_ID'] )
}	
##--------------------------------------------------------------------------------------------------------
##	predict Calls for contigs with same PANGEA ID, based on fitted model 
##--------------------------------------------------------------------------------------------------------
haircut.get.call.for.PNG_ID<- function(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)	
{
	#	load covariates
	cnsc.df	<- do.call('rbind',lapply(files, function(x)
					{
						load( paste(indir.st, '/', x, sep='') )
						cnsc.df
					})) 
	#	load alignment
	crs		<- lapply(alfiles, function(x){
				cr	<- read.dna(file=paste(indir.al,x,sep='/'), format='fasta')
				cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]										
			})
	names(crs)	<- bc
	#	predict
	tmp		<- seq(cnsc.df[, floor(min(SITE)/10)*10],cnsc.df[, max(SITE)+10],10)
	cnsc.df[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
	cnsc.df	<- merge(cnsc.df, ctrmc, by='CHUNK')
	cnsc.df[, PR_CALL:= predict.fun(AGRpc, GPS, BETA0, BETA1, BETA2)]							
	cnsc.df[, CALL:= as.integer(PR_CALL>=par['PRCALL.thr'])]	
	#	check if called contig has gaps of CALL=='0': if yes, return last non-gap before first CALL=='0
	if(!is.na(par['PRCALL.cutprdcthair']))
	{
		tmp		<- cnsc.df[, {
					z	<- gsub('0*$','',paste(CALL,collapse=''))
					#print(TAXON)
					#print(z)
					z	<- gregexpr('1+',z)[[1]]	
					list(CALL_POS=as.integer(z+min(SITE)-1L), CALL_LEN=attr(z, 'match.length'))
				}, by=c('PNG_ID','TAXON','BLASTnCUT')]
		tmp		<- subset( tmp, CALL_LEN<par['PRCALL.cutprdcthair'] )
		if(nrow(tmp))
		{
			cat('\nFound predicted extra hair of length <',par['PRCALL.cutprdcthair'],'delete, n=',tmp[,sum(CALL_LEN)])
			set(tmp, NULL, 'CALL_LEN', tmp[, CALL_POS+CALL_LEN-1])
			for(i in seq_len(nrow(tmp)))
			{				
				set(cnsc.df, cnsc.df[, which(TAXON==tmp$TAXON[i] & BLASTnCUT==tmp$BLASTnCUT[i] & SITE>=tmp$CALL_POS[i] & SITE<=tmp$CALL_LEN[i])], 'CALL', 0L)
			}				
		}										
	}
	#	produce fasta output:
	#	select cut and raw contigs with a call
	crs			<- lapply(crs, as.character)
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	tmp			<- rownames(crs[['N']])[ !grepl(png_id,rownames(crs[['N']])) | rownames(crs[['N']])%in%tmp ] 
	crs[['N']]	<- crs[['N']][tmp,]
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	tmp			<- rownames(crs[['Y']])[ !grepl(png_id,rownames(crs[['Y']])) | rownames(crs[['Y']])%in%tmp ] 
	crs[['Y']]	<- crs[['Y']][tmp,]
	#	set all characters with CALL==0 to -
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	for(tx in tmp)
	{
		crs[['N']][tx, subset(cnsc.df, BLASTnCUT=='N' & TAXON==tx & CALL==0)[, SITE]]	<- '-'
	}
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	for(tx in tmp)
	{
		crs[['Y']][tx, subset(cnsc.df, BLASTnCUT=='Y' & TAXON==tx & CALL==0)[, SITE]]	<- '-'
	}
	crs			<- lapply(crs, as.DNAbin)		
	
	list(crs=crs, cnsc.df=cnsc.df)
}
##--------------------------------------------------------------------------------------------------------
##	predict Calls for contigs with same PANGEA ID, based on fitted model
##	update: 
##		1)	do not return duplicate contigs (ie cut and raw, if both are to be kept)
##		2)	do not return raw contigs if cut exists and if raw extends into LTR
##--------------------------------------------------------------------------------------------------------
haircut.get.call.for.PNG_ID.150814<- function(indir.str, indir.al, png_id, files, alfiles, bc, par, ctrmc, predict.fun)	
{
	#	load covariates
	cnsc.df	<- do.call('rbind',lapply(files, function(x)
					{
						load( paste(indir.st, '/', x, sep='') )
						cnsc.df
					})) 
	#	load alignment	and cut LTR and anything that extends past references
	crs		<- lapply(alfiles, function(x)
					{
						cr	<- read.dna(file=paste(indir.al,x,sep='/'), format='fasta')
						cr	<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
						cr[, seq.int(1, haircut.find.lastRefSite(cr))]																
					})
	names(crs)	<- bc
	#
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)
					{
						tmp		<- rownames(crs[[i]])[grepl(png_id,rownames(crs[[i]]))]						
						data.table(	TAXON=tmp, 
								BLASTnCUT= factor(blastncut[i],levels=c('cut','raw'),labels=c('Y','N')), 
								FIRST= apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]][tmp,,drop=FALSE]), 1, function(x) which(rev(x)!='-')[1] ),
								CRS_ID=i)
					}))	
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
	tx[, OCNTG:= tx[, sapply(strsplit(CNTG,'.',fixed=T),'[[',1)]]
	tx[, CCNTG:= NA_character_]		
	tmp		<- tx[, which(grepl('.',CNTG,fixed=T))]
	if(length(tmp))
		set(tx, tmp, 'CCNTG', tx[tmp, sapply(strsplit(CNTG,'.',fixed=T),'[[',2)])	
	tmp		<- subset(tx, BLASTnCUT=='Y' & !is.na(CCNTG))[, list(CCNTGn=length(CCNTG)), by='OCNTG']
	tmp		<- subset(tmp, CCNTGn==1)[, OCNTG]	#check for cut contigs that should be present in multiple cuts by naming scheme, but after LTR removal there is only one cut
	if(length(tmp))
	{
		cat('\nFound lone cuts for which multiple cuts are expected by naming scheme, n=', length(tmp))
		tmp	<- tx[, which(BLASTnCUT=='Y' & OCNTG%in%tmp)]
		set(tx, tmp, 'CCNTG', NA_character_)
		set(tx, tmp, 'CNTG', tx[tmp, OCNTG])
	}
	#	predict
	tmp		<- seq(cnsc.df[, floor(min(SITE)/10)*10],cnsc.df[, max(SITE)+10],10)
	cnsc.df[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]	
	cnsc.df	<- merge(cnsc.df, ctrmc, by='CHUNK')
	cnsc.df[, PR_CALL:= predict.fun(AGRpc, GPS, BETA0, BETA1, BETA2)]							
	cnsc.df[, CALL:= as.integer(PR_CALL>=par['PRCALL.thr'])]	
	#	determine predicted sites
	cnsc.1s	<- cnsc.df[, {
				z	<- gsub('0*$','',paste(CALL,collapse=''))
				#print(TAXON)
				#print(z)
				z	<- gregexpr('1+',z)[[1]]	
				list(CALL_ID= seq_along(z), CALL_POS=as.integer(z+min(SITE)-1L), CALL_LEN=attr(z, 'match.length'))
			}, by=c('PNG_ID','TAXON','BLASTnCUT')]
	cnsc.1s[, CALL_LAST:= CALL_POS+CALL_LEN-1L]
	cnsc.1s	<- subset(cnsc.1s, CALL_LEN>0)
	#	fill internal predicted gaps		
	if(!is.na(par['PRCALL.rmintrnlgpsblw'] | !is.na(par['PRCALL.rmintrnlgpsend'])))
	{
		cnsc.g	<- cnsc.1s[, {
					if(length(CALL_ID)==1)
						ans	<- NA_integer_
					if(length(CALL_ID)>1)
						ans	<- CALL_POS[seq.int(2,length(CALL_POS))]-CALL_LAST[seq.int(1,length(CALL_LAST)-1)]-1L
					list(CALL_LAST=CALL_LAST[seq.int(1,length(CALL_LAST)-1)], GAP_LEN= ans)
				}, by=c('TAXON','BLASTnCUT')]
		cnsc.g	<- merge(cnsc.1s, cnsc.g,  by=c('TAXON','BLASTnCUT','CALL_LAST'), all.x=1)
		tmp		<- cnsc.g[, which(GAP_LEN<par['PRCALL.rmintrnlgpsblw'] | (CALL_LAST>par['PRCALL.rmintrnlgpsend'] & !is.na(GAP_LEN)))]
		for(i in tmp)	#	add ith called region to next call region
		{
			tmp2	<- cnsc.df[, which(TAXON==cnsc.g$TAXON[i] & BLASTnCUT==cnsc.g$BLASTnCUT[i] & SITE>cnsc.g$CALL_LAST[i] & SITE<=cnsc.g$CALL_LAST[i]+cnsc.g$GAP_LEN[i])]
			stopifnot( cnsc.df[tmp2,all(CALL==0)])
			cat('\nFound internal predicted gap and set to CALL=1, taxon=',cnsc.g[i,TAXON],', cut=',cnsc.g[i,as.character(BLASTnCUT)],', pos=',cnsc.g[i,CALL_LAST+1L],', len=', length(tmp2))
			set(cnsc.df, tmp2, 'CALL', 1L)
			set(cnsc.g, i+1L, 'CALL_POS', cnsc.g[i,CALL_POS])
			set(cnsc.g, i+1L, 'CALL_LEN', cnsc.g[i+1L,CALL_LEN]+cnsc.g[i,CALL_LEN]+cnsc.g[i,GAP_LEN])
			set(cnsc.g, i, 'CALL_ID', NA_integer_)
		}
		cnsc.1s		<- subset(cnsc.g, !is.na(CALL_ID))
	}
	#	check if all called chunks in cut and raw contigs correspond to each other
	if(!is.na(par['PRCALL.cutrawgrace']))
	{
		cnsc.1s	<- merge(cnsc.1s, tx, by=c('TAXON','BLASTnCUT'))
		tmp		<- subset(cnsc.1s, BLASTnCUT=='Y', select=c(OCNTG, TAXON, CALL_ID, CALL_POS, CALL_LEN ))
		setnames(tmp, c('TAXON','CALL_ID','CALL_POS','CALL_LEN'),c('TAXON_CUT','CALL_ID_CUT','CALL_POS_CUT','CALL_LEN_CUT'))
		tmp		<- merge(subset(cnsc.1s, BLASTnCUT=='N'), tmp, by='OCNTG', allow.cartesian=TRUE)
		tmp		<- subset(tmp, abs(CALL_POS-CALL_POS_CUT)<par['PRCALL.cutrawgrace'] )	
		#	of the corresponding calls, keep the longer one 
		#	dont extend shorter for now as alignments may not match
		tmp2	<- subset(tmp, CALL_LEN>CALL_LEN_CUT)
		for(i in seq_len(nrow(tmp2)))	#keep raw
		{			
			cat('\nkeep only raw:', tmp2[i,TAXON])
			z	<- cnsc.df[, which( TAXON==tmp2$TAXON_CUT[i] & BLASTnCUT=='Y' & SITE>=tmp2$CALL_POS_CUT[i] & SITE<(tmp2$CALL_POS_CUT[i]+tmp2$CALL_LEN_CUT[i]))]
			stopifnot( cnsc.df[z,all(CALL==1)])
			set(cnsc.df, z, 'CALL', 0L)		
		}
		if(length(tmp2))
			set(cnsc.1s, cnsc.1s[, which(TAXON%in%tmp2[, TAXON_CUT] & BLASTnCUT=='Y' & CALL_ID%in%tmp2[, CALL_ID_CUT])],'CALL_ID',NA_integer_)	
		tmp2	<- subset(tmp, CALL_LEN<=CALL_LEN_CUT)
		for(i in seq_len(nrow(tmp2)))	#keep cut
		{	
			cat('\nkeep only cut:', tmp2[i,TAXON_CUT])
			z	<- cnsc.df[, which( TAXON==tmp2$TAXON[i] & BLASTnCUT=='N' & SITE>=tmp2$CALL_POS[i] & SITE<=tmp2$CALL_LAST[i])]
			stopifnot( cnsc.df[z,all(CALL==1)])
			set(cnsc.df, z, 'CALL', 0L)
		}
		if(length(tmp2))
			set(cnsc.1s, cnsc.1s[, which(TAXON%in%tmp2[, TAXON] & BLASTnCUT=='N' & CALL_ID%in%tmp2[, CALL_ID])],'CALL_ID',NA_integer_)
		cnsc.1s	<- subset(cnsc.1s, !is.na(CALL_ID))		
	}	
	#	check if called contig has gaps of CALL=='0': if yes, return last non-gap before first CALL=='0
	if(!is.na(par['PRCALL.cutprdcthair']))
	{
		tmp		<- subset( cnsc.1s, CALL_LEN<par['PRCALL.cutprdcthair'] )
		if(nrow(tmp))
		{
			cat('\nFound predicted extra hair of length <',par['PRCALL.cutprdcthair'],'delete, n=',tmp[,sum(CALL_LEN)])
			set(tmp, NULL, 'CALL_LEN', tmp[, CALL_POS+CALL_LEN-1])
			for(i in seq_len(nrow(tmp)))
			{				
				set(cnsc.df, cnsc.df[, which(TAXON==tmp$TAXON[i] & BLASTnCUT==tmp$BLASTnCUT[i] & SITE>=tmp$CALL_POS[i] & SITE<=tmp$CALL_LEN[i])], 'CALL', 0L)
			}				
		}										
	}
	#	check there is no dust
	tmp			<- subset(cnsc.df, CALL==1)[, list(CALL_N= length(CALL)), by=c('TAXON','BLASTnCUT')][, CALL_N]
	if(length(tmp))
		stopifnot(min(tmp)>40)
	#	produce fasta output:
	#	select cut and raw contigs with a call, then set all characters with CALL==0 to -
	crs			<- lapply(crs, as.character)
	tmp			<- subset(cnsc.df, BLASTnCUT=='N' & CALL==1 )[, unique(TAXON)]
	tmp2		<- rownames(crs[['N']])[ !grepl(png_id,rownames(crs[['N']])) | rownames(crs[['N']])%in%tmp ] 
	crs[['N']]	<- crs[['N']][tmp2,]
	for(tx in tmp)
		crs[['N']][tx, subset(cnsc.df, BLASTnCUT=='N' & TAXON==tx & CALL==0 & SITE<=ncol(crs[['N']]))[, SITE]]	<- '-'
			
	tmp			<- subset(cnsc.df, BLASTnCUT=='Y' & CALL==1 )[, unique(TAXON)]
	tmp2		<- rownames(crs[['Y']])[ !grepl(png_id,rownames(crs[['Y']])) | rownames(crs[['Y']])%in%tmp ] 
	crs[['Y']]	<- crs[['Y']][tmp2,]
	for(tx in tmp)
		crs[['Y']][tx, subset(cnsc.df, BLASTnCUT=='Y' & TAXON==tx & CALL==0 & SITE<=ncol(crs[['Y']]))[, SITE]]	<- '-'
	crs			<- lapply(crs, as.DNAbin)			
	list(crs=crs, cnsc.df=cnsc.df)
}
##--------------------------------------------------------------------------------------------------------
##	fit Beta Binomial regression model to training data 
##--------------------------------------------------------------------------------------------------------
haircut.get.fitted.model.150814a<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		ctrmc	<- do.call('rbind', lapply(seq(1,10001,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							tmp		<- tmp[tmp<=10000 | tmp==floor(ctr[, max(SITE)+10]/10)*10]
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]							
							ctrmc	<- do.call('rbind',lapply(ctr[, unique(CHUNK)], function(chunk)
											{
												ctrch		<- subset(ctr, CHUNK==chunk)
												tmp			<- tryCatch( gamlss(ANS_CALL~AGRpc+GPS, sigma.formula=~AGRpc+GPS, data=ctrch, family=BB(), control=gamlss.control(trace=FALSE), i.control=glim.control(cyc=100,cc=1e-4)), warning=function(w) w, error=function(e) e)
												if(!is(tmp,'warning') & !is(tmp,'error'))
												{
													ctrchm	<- tmp
													tmp		<- tryCatch(vcov(ctrchm),warning=function(w) w, error=function(e) e)
													tmp2	<- !is(tmp,'warning') & !is(tmp,'error')
												}
												if(is(tmp,'warning') | is(tmp,'error'))
												{
													ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, sigma.formula=~AGRpc+GPS, data=ctrch, family=BI(), control=gamlss.control(trace=FALSE))
													tmp2	<- TRUE
												}													
												ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3], FAMILY=family(ctrchm)[1], CONVERGED=tmp2	)					
											}))
						}))
		#	deal with end of genome where little data is available: we estimated coefs for all sites > 1e4, now split up using CHUNK notation
		ctr		<- haircut.load.training.data(indir, 10001)
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
		tmp[tmp>10000]
		tmp		<- as.data.table(expand.grid(CHUNK=as.character(tmp[tmp>10000]), BETA0=subset(ctrmc, CHUNK==10000)[, BETA0], BETA1=subset(ctrmc, CHUNK==10000)[, BETA1], BETA2=subset(ctrmc, CHUNK==10000)[, BETA2], stringsAsFactors=F))
		ctrmc	<- rbind(ctrmc, tmp)
		#	model predict function, so we save mem by not having to call 'predict'
		model.150811a.predict<- function(agrpc, gps, b0, b1, b2)
		{	
			stopifnot(all(!is.na(agrpc)), all(!is.na(gps)))
			b0[which(is.na(b0))]	<- 0
			b1[which(is.na(b1))]	<- 0
			b2[which(is.na(b2))]	<- 0
			exp(b0+b1*agrpc+b2*gps)/(exp(b0+b1*agrpc+b2*gps)+1)	
		}
		#	calculate Sensitivity & Specificity on training data
		ctrev	<- do.call('rbind', lapply(seq(1,11000,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
							print(ctr)
							ctrp	<- merge(ctr, ctrmc, by='CHUNK')
							ctrp[, PR_CALL:=model.150811a.predict(AGRpc, GPS, BETA0, BETA1, BETA2)]
							setkey(ctrp, CHUNK)
							ctrev	<- as.data.table(expand.grid(THR= seq(0.05,0.95,0.01), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
							ctrev	<- ctrev[, {
										tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=THR, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
										list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
									}, by=c('CHUNK','THR')]
							ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))
							#	plot Sensitivity & Specificity
							ggplot(melt(ctrev, id.vars=c('CHUNK','THR'), measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=100*value, colour=CHUNK)) + 
									geom_line() + labs(x='threshold on predicted call probability',y='%', colour='site\n(base relative to HXB2)') +
									theme(legend.position='bottom') + facet_wrap(~variable, ncol=2, scales='free') +
									guides(col = guide_legend(ncol=5, byrow=TRUE))
							ggsave(file=paste(outdir,'/',outfile,'_model.150811a_SensSpec_SITE',site,'.pdf',sep=''), w=9, h=9)
							#
							ctrev							
						}))		
		set(ctrev, NULL, 'CHUNK', ctrev[, as.numeric(CHUNK)])
		
		save(ctrmc, ctrev, model.150811a.predict, file=outfile)
	}
	list(coef=ctrmc, ev=ctrev, predict=model.150811a.predict)
}
##--------------------------------------------------------------------------------------------------------
##	fit simple Binomial regression model to training data 
##--------------------------------------------------------------------------------------------------------
haircut.get.fitted.model.150811a<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		ctrmc	<- do.call('rbind', lapply(seq(1,10001,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							tmp		<- tmp[tmp<=10000 | tmp==floor(ctr[, max(SITE)+10]/10)*10]
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]							
							ctrmc	<- do.call('rbind',lapply(ctr[, unique(CHUNK)], function(chunk)
											{
												ctrch	<- subset(ctr, CHUNK==chunk)
												ctrchm	<- gamlss(ANS_CALL~AGRpc+GPS, data=ctrch, family=BI())				#as good as 'AGRpc+GPS+CNS_FRQr' in terms of FN, FP
												ctrchmc	<- data.table(CHUNK=chunk, BETA0=coef(ctrchm)[1], BETA1=coef(ctrchm)[2], BETA2=coef(ctrchm)[3])					
											}))
						}))
		#	deal with end of genome where little data is available: we estimated coefs for all sites > 1e4, now split up using CHUNK notation
		ctr		<- haircut.load.training.data(indir, 10001)
		tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
		tmp[tmp>10000]
		tmp		<- as.data.table(expand.grid(CHUNK=as.character(tmp[tmp>10000]), BETA0=subset(ctrmc, CHUNK==10000)[, BETA0], BETA1=subset(ctrmc, CHUNK==10000)[, BETA1], BETA2=subset(ctrmc, CHUNK==10000)[, BETA2], stringsAsFactors=F))
		ctrmc	<- rbind(ctrmc, tmp)
		#	model predict function, so we save mem by not having to call 'predict'
		model.150811a.predict<- function(agrpc, gps, b0, b1, b2)
			{	
				stopifnot(all(!is.na(agrpc)), all(!is.na(gps)))
				b0[which(is.na(b0))]	<- 0
				b1[which(is.na(b1))]	<- 0
				b2[which(is.na(b2))]	<- 0
				exp(b0+b1*agrpc+b2*gps)/(exp(b0+b1*agrpc+b2*gps)+1)	
			}
		#	calculate Sensitivity & Specificity on training data
		ctrev	<- do.call('rbind', lapply(seq(1,11000,200), function(site)
						{
							cat('\nProcess',site,'\n')							
							ctr		<- haircut.load.training.data(indir, site)	
							tmp		<- seq(ctr[, min(SITE)-1],ctr[, max(SITE)+10],10)
							ctr[, CHUNK:=cut(SITE, breaks=tmp, labels=tmp[-length(tmp)])]
							print(ctr)
							ctrp	<- merge(ctr, ctrmc, by='CHUNK')
							ctrp[, PR_CALL:=model.150811a.predict(AGRpc, GPS, BETA0, BETA1, BETA2)]
							setkey(ctrp, CHUNK)
							ctrev	<- as.data.table(expand.grid(THR= seq(0.05,0.95,0.01), CHUNK= as.character(ctrp[, unique(CHUNK)]), stringsAsFactors=FALSE))
							ctrev	<- ctrev[, {
										tmp	<- ctrp[CHUNK,][,table(ANS_CALL, factor(PR_CALL>=THR, levels=c('TRUE','FALSE'), labels=c('TRUE','FALSE')))]
										list(TP= tmp['1','TRUE'], FP= tmp['0','TRUE'], FN= tmp['1','FALSE'], TN= tmp['0','FALSE'] )					
									}, by=c('CHUNK','THR')]
							ctrev	<- merge(ctrev, ctrev[, list(SENS= TP/(TP+FN), SPEC=TN/(TN+FP), FDR=FP/(FP+TP), FOR=FN/(FN+TN)), by=c('CHUNK','THR')], by=c('CHUNK','THR'))
							#	plot Sensitivity & Specificity
							ggplot(melt(ctrev, id.vars=c('CHUNK','THR'), measure.vars=c('SENS','SPEC','FDR','FOR')), aes(x=THR, y=100*value, colour=CHUNK)) + 
									geom_line() + labs(x='threshold on predicted call probability',y='%', colour='site\n(base relative to HXB2)') +
									theme(legend.position='bottom') + facet_wrap(~variable, ncol=2, scales='free') +
									guides(col = guide_legend(ncol=5, byrow=TRUE))
							ggsave(file=paste(outdir,'/',outfile,'_model.150811a_SensSpec_SITE',site,'.pdf',sep=''), w=9, h=9)
							#
							ctrev							
						}))		
		set(ctrev, NULL, 'CHUNK', ctrev[, as.numeric(CHUNK)])
		
		save(ctrmc, ctrev, model.150811a.predict, file=outfile)
	}
	list(coef=ctrmc, ev=ctrev, predict=model.150811a.predict)
}
##--------------------------------------------------------------------------------------------------------
##	determine if cut contigs are identical to a stretch or possibly the whole raw contig with the same contig ID 
##--------------------------------------------------------------------------------------------------------
haircutwrap.get.subset.among.raw.and.cut.contigs<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs, and identify cut contigs that are part of the raw contig
		txe		<- infiles[,{
					#png_id	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID==png_id)[, FILE]
					#blastncut<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.subset.among.raw.and.cut.contigs(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		save(txe, file=outfile)
	}
	txe
}
##--------------------------------------------------------------------------------------------------------
##	get contigs to be used for training. this excludes duplicate raw/cut contigs for which a curated answer is available 
##--------------------------------------------------------------------------------------------------------
haircut.get.training.contigs<- function(indir, outfile, ctrain)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs, and identify cut contigs that are part of the raw contig
		txe		<- infiles[,{
					#png_id	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID==png_id)[, FILE]
					#blastncut<- subset(infiles, PNG_ID==png_id)[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.subset.among.raw.and.cut.contigs(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		#	consider identical contigs. merge cut contigs into equal contigs, so we can look at all equal pairs
		txec	<- merge(subset(txe, EQ), ctrain, all.x=TRUE, by=c('PNG_ID','TAXON'))		
		tmp		<- subset(txe, EQ & CUT=='cut', select=c(PNG_ID, OCNTG, CCNTG, TAXON, EQ_FIRST, EQ_LAST))
		setnames(tmp, c('TAXON','CCNTG','EQ_FIRST','EQ_LAST'),c('TAXON_CUT','CCNTG_CUT','EQ_FIRST_CUT','EQ_LAST_CUT'))
		txec	<- merge(txec, tmp, all.x=TRUE, allow.cartesian=TRUE,by=c('PNG_ID','OCNTG'))
		tmp		<- subset(ctrain, select=c(PNG_ID, TAXON, ANS_LEN))
		setnames(tmp, c('TAXON','ANS_LEN'), c('TAXON_CUT','ANS_LEN_CUT'))
		txec	<- merge(txec, tmp, all.x=1, by=c('PNG_ID','TAXON_CUT'))
		#	work out which contigs to use for training
		txec[, USE_IN_TRAIN:='Y']
		#	if raw and cut contigs in curated and identical: don t use raw (double counting).
		tmp		<- txec[, which(CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound raw and cut contigs in curated that are identical: don t use raw for training to avoid double counting, n=', length(tmp)))
		set(txec, tmp, 'ANS_FILE', NA_character_)
		set(txec, tmp, c('ANS_FIRST','ANS_LAST'), NA_integer_)	#need to set ANS_LEN to NA later
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#			
		#	if raw and concatenated cut contigs identical and only raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, USE_IN_TRAIN=='Y' & CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and cut contigs identical and only raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and only cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='raw' & is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and cut contigs identical and only cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#
		#	if cut contig subset of raw contig and both in curated: should not happen
		tmp		<- nrow(subset(txec, USE_IN_TRAIN=='Y' & CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT)))
		stopifnot(tmp==0)
		#	if cut contig subset of raw contig and only cut in curated: use cut only partial
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='cut' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: use cut only partial, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'P')
		#	if cut contig subset of raw contig and only cut in curated: dont use raw as 0's
		tmp		<- txec[, which(USE_IN_TRAIN=='Y' & CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: dont use raw as 0s, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	if cut contig subset of raw contig and only cut in curated: use cut only partial
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))[, TAXON_CUT]
		tmp		<- txec[, which(TAXON%in%tmp & CUT=='cut')]
		cat(paste('\nFound cut contig subset of raw contig and only cut in curated: dont use raw as 0s, n=', length(tmp)))		
		set(txec, tmp, 'USE_IN_TRAIN', 'P')
		#	can now set ANS_LEN to NA from above case 
		set(txec, txec[, which(is.na(ANS_FIRST) & !is.na(ANS_LEN))], 'ANS_LEN', NA_integer_)
		#	add contigs that are in curated and unique amongst cut and raw
		tmp		<- merge(subset(txe, !EQ), ctrain, by=c('PNG_ID','TAXON'))
		tmp[, USE_IN_TRAIN:='Y']
		stopifnot( length(intersect( tmp[, TAXON], txec[,TAXON] ))==0	)
		txec	<- rbind(txec, tmp, use.names=TRUE, fill=TRUE)
		#	if cut contig subset of raw contig and only raw in curated: dont use cut as 0's
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(EQ_FIRST_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT), c(PNG_ID, OCNTG))
		tmp		<- subset( merge(tmp, txe, by=c('PNG_ID','OCNTG')), CUT=='cut' )[, unique(TAXON)]
		cat(paste('\nFound cut contig subset of raw contig and only raw in curated: dont use cut as 0, n=', length(tmp)))
		tmp		<- txec[, which(TAXON%in%tmp & CUT=='cut')]
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	delete tmp cols
		set(txec, NULL, c('TAXON_CUT','CCNTG_CUT','ANS_LEN_CUT','EQ_FIRST_CUT','EQ_LAST_CUT'), NULL)
		setkey(txec, PNG_ID, TAXON, CUT)
		txec	<- unique(txec)					
		save(txec, file=outfile)
		#print( txec[, table(USE_IN_TRAIN)] )
	}
	txec
}

haircut.get.training.contigs.byidentical<- function(indir, outfile, ctrain)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)	
	if( inherits(readAttempt, "try-error")	)
	{
		infiles	<- data.table(FILE=list.files(indir, pattern='fasta$', recursive=T))
		infiles[, PNG_ID:= gsub('_wRefs\\.fasta','',gsub('_cut|_raw','',basename(FILE)))]
		infiles[, BLASTnCUT:= regmatches(FILE,regexpr('cut|raw',FILE))]		 
		#	identify identical raw and cut contigs
		txe		<- infiles[,{
					#PNG_ID	<- '12559_1_11'
					#files	<- subset(infiles, PNG_ID=='12559_1_10')[, FILE]
					#blastncut<- subset(infiles, PNG_ID=='12559_1_10')[, BLASTnCUT]
					cat('\nProcess', PNG_ID)
					tx	<- haircut.get.identical.among.raw.and.cut.contigs(indir, FILE, PNG_ID, BLASTnCUT)
					tx
				}, by='PNG_ID']
		#	consider identical contigs. merge cut contigs into equal contigs, so we can look at all equal pairs
		txec	<- merge(subset(txe, EQ), ctrain, all.x=TRUE, by=c('PNG_ID','TAXON'))		
		tmp		<- subset(txe, EQ & CUT=='cut', select=c(PNG_ID, OCNTG, CCNTG, TAXON))
		setnames(tmp, c('TAXON','CCNTG'),c('TAXON_CUT','CCNTG_CUT'))
		txec	<- merge(txec, tmp, all.x=TRUE, allow.cartesian=TRUE,by=c('PNG_ID','OCNTG'))
		tmp		<- subset(ctrain, select=c(PNG_ID, TAXON, ANS_LEN))
		setnames(tmp, c('TAXON','ANS_LEN'), c('TAXON_CUT','ANS_LEN_CUT'))
		txec	<- merge(txec, tmp, all.x=1, by=c('PNG_ID','TAXON_CUT'))
		#	work out which contigs to use for training
		txec[, USE_IN_TRAIN:='Y']
		#	if raw and cut contigs in curated and identical: don t use raw (double counting).
		tmp		<- txec[, which(CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		cat(paste('\nFound raw and cut contigs in curated that are identical: don t use raw for training to avoid double counting, n=', length(tmp)))
		set(txec, tmp, 'ANS_FILE', NA_character_)
		set(txec, tmp, c('ANS_FIRST','ANS_LAST'), NA_integer_)	#need to set ANS_LEN to NA later
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs in curated and identical: should not happen
		tmp		<- txec[, which(CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))]
		stopifnot(length(tmp)==0)
		#	if raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, CUT=='raw' & !is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(CUT=='raw' & !is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and concatenated cut contigs identical and cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')		
		#	if raw and cut contigs identical and cut in curated: don t use raw in training as 0
		tmp		<- txec[, which(CUT=='raw' & is.na(CCNTG_CUT) & is.na(ANS_LEN) & !is.na(ANS_LEN_CUT))] 
		cat(paste('\nFound raw and cut contigs identical and cut in curated: don t use raw in training as 0, n=', length(tmp)))
		set(txec, tmp, 'USE_IN_TRAIN', 'N')
		#	if raw and concatenated cut contigs identical and raw in curated: don t use cut in training as 0
		tmp		<- subset(txec, CUT=='raw' & is.na(CCNTG_CUT) & !is.na(ANS_LEN) & is.na(ANS_LEN_CUT))[, TAXON_CUT] 
		cat(paste('\nFound raw and cut contigs identical and raw in curated: don t use cut in training as 0, n=', length(tmp)))
		set(txec, txec[, which(TAXON%in%tmp & CUT=='cut')], 'USE_IN_TRAIN', 'N')
		#	can now set ANS_LEN to NA from above case 
		set(txec, txec[, which(is.na(ANS_FIRST) & !is.na(ANS_LEN))], 'ANS_LEN', NA_integer_)
		#	delete tmp rows for cut-cut
		#subset(txec, CUT=='raw')
		#tmp		<- subset(txec, CUT=='cut' )[, {
		#			stopifnot( all(USE_IN_TRAIN==USE_IN_TRAIN[1])	)
		#			list(CUT=CUT[1], FIRST=FIRST[1], LAST=LAST[1], CRS_ID=CRS_ID[1], CNTG=CNTG[1], OCNTG=OCNTG[1], CCNTG=CCNTG[1], EQ=EQ[1], INFILE=INFILE[1], ANS_FILE=ANS_FILE[1], ANS_LEN=ANS_LEN[1], ANS_FIRST=ANS_FIRST[1], ANS_LAST=ANS_LAST[1], TAXON_CUT=TAXON_CUT[1], CCNTG_CUT=CCNTG_CUT[1], ANS_LEN_CUT=ANS_LEN_CUT[1], USE_IN_TRAIN=USE_IN_TRAIN[1])
		#		}, by=c('PNG_ID','TAXON')]
		#subset(txec, CUT=='raw')		
		#	delete tmp cols
		set(txec, NULL, c('TAXON_CUT','CCNTG_CUT','ANS_LEN_CUT'), NULL)
		setkey(txec, PNG_ID, TAXON, CUT)
		txec	<- unique(txec)
		#
		#	add contigs that are in curated and unique amongst cut and raw
		#
		tmp		<- merge(subset(txe, !EQ), ctrain, by=c('PNG_ID','TAXON'))
		tmp[, USE_IN_TRAIN:='Y']
		stopifnot( length(intersect( tmp[, TAXON], txec[,TAXON] ))==0	)
		txec	<- rbind(txec, tmp, use.names=TRUE)
		save(txec, file=outfile)
		#print( txec[, table(USE_IN_TRAIN)] )
	}	
	txec
}
##--------------------------------------------------------------------------------------------------------
##	return EQ= 1 and EQ_FIRST EQ_LAST positions, if cut contig are subsets of the raw contig in that: 
##	- 	there is only once cut contig which disagrees on up to x positions, where x is the difference in alignment lengths in the cut and raw files
##	- 	the cut contigs matches the raw contig without gaps, where the raw contig is set to start from the start position of the cut contig
##--------------------------------------------------------------------------------------------------------
haircut.get.subset.among.raw.and.cut.contigs<- function(indir, files, png_id, blastncut)
{
	crs		<- lapply(files, function(x)
			{
				cr		<- read.dna(file=paste(indir,'/',x,sep=''), format='fasta')
				if(is.matrix(cr))
				{
					cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
					tmp		<- strsplit(basename(x), '_')[[1]][1]
					tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
					cr		<- cr[ subset(tx, CONTIG==1)[, TAXON], ]	
				}
				if(!is.matrix(cr))
					cr		<- as.DNAbin(matrix(vector('character',0), nrow=2, byrow=T, dimnames=list(c('DUMMY','DUMMY'),c())))
				cr					
			})
	names(crs)	<- blastncut
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)	data.table(	TAXON=rownames(crs[[i]]), 
								CUT= blastncut[i], 
								FIRST= apply( as.character(crs[[i]]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]]), 1, function(x) which(rev(x)!='-')[1] ),
								CRS_ID=i)	))
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
	tx[, OCNTG:= tx[, sapply(strsplit(CNTG,'.',fixed=T),'[[',1)]]
	tx[, CCNTG:= NA_character_]		
	tx[, EQ:=FALSE]
	tx[, EQ_FIRST:=NA_integer_]
	tx[, EQ_LAST:=NA_integer_]
	tmp		<- tx[, which(grepl('.',CNTG,fixed=T))]
	if(length(tmp))
		set(tx, tmp, 'CCNTG', tx[tmp, sapply(strsplit(CNTG,'.',fixed=T),'[[',2)])	
	tmp		<- subset(tx, CUT=='cut' & !is.na(CCNTG))[, list(CCNTGn=length(CCNTG)), by='OCNTG']
	tmp		<- subset(tmp, CCNTGn==1)[, OCNTG]	#check for cut contigs that should be present in multiple cuts by naming scheme, but after LTR removal there is only one cut
	if(length(tmp))
	{
		cat('\nFound lone cuts for which multiple cuts are expected by naming scheme, n=', length(tmp))
		tmp	<- tx[, which(CUT=='cut' & OCNTG%in%tmp)]
		set(tx, tmp, 'CCNTG', NA_character_)
		set(tx, tmp, 'CNTG', tx[tmp, OCNTG])
	}#	check if there are contigs with corresponding name in 'cut' and that are identical / subset
	#	differences in gaps are allowed: these will come from different alignment
	txe		<- dcast.data.table(tx, CNTG~CUT, value.var='TAXON')
	txe		<- subset( txe, !is.na(cut) & !is.na(raw) )
	if(nrow(txe) && c('cut','raw')%in%colnames(txe))
	{						
		txe		<- txe[, {
					x							<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['cut']][cut,])),collapse='')))
					y							<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][raw,])),collapse='')))
					#print(CNTG)
					#print(x)
					#print(y)
					z							<- nchar(x)
					if(nchar(x)==nchar(y))
					{						
						z2	<- x==y
						z3	<- 1L
						attr(z3, 'match.length')<- nchar(y)
					}						
					if(nchar(x)!=nchar(y))
					{
						z	<- gsub('-','',x)
						z2	<- z==gsub('-','',y)
						z3	<- haircut.large.regexpr(z,gsub('-','',y))
						z	<- nchar(x)
					}						
					list(EQtmp=z2, EQ_FIRSTtmp=as.integer(z3), EQ_LASTtmp=as.integer(z3+attr(z3, 'match.length')-1L+nchar(x)-z))					
				}, by='CNTG']
		tx		<- merge(tx, txe, all.x=TRUE, by='CNTG')
		tmp		<- tx[, which(!is.na(EQ_FIRSTtmp) & EQ_FIRSTtmp>0)]
		set(tx, tmp, 'EQ_FIRST', tx[tmp,FIRST]+tx[tmp,EQ_FIRSTtmp]-1L)
		set(tx, tmp, 'EQ_LAST', tx[tmp,FIRST]+tx[tmp,EQ_LASTtmp]-1L-tx[tmp,EQ_FIRSTtmp])
		set(tx, tmp, 'EQ', tx[tmp,EQtmp])
		set(tx, NULL, c('EQ_FIRSTtmp','EQ_LASTtmp','EQtmp'),NULL)				
	}
	#	see if cut contigs are subset of the raw contig
	txe		<- subset(tx, !is.na(CCNTG))
	if(nrow(txe))
	{
		tmp		<- subset(tx, CUT=='raw')
		setnames(tmp, colnames(tmp)[ !grepl('OCNTG',colnames(tmp))], paste( colnames(tmp)[ !grepl('OCNTG',colnames(tmp))], '_raw',sep='' ))
		txe		<- merge(txe, tmp, by='OCNTG')
		txe		<- txe[, {
					z		<- sub('-*$','',paste(as.vector(as.character(crs[['cut']][TAXON,])),collapse=''))
					x		<- sub('^-*','',z)
					tmp		<- nchar(z)-nchar(x)	#number '-' clipped up to start of cut contig
					tmp		<- max(1L, 1L+tmp-abs(diff(c(ncol(crs[['cut']]), ncol(crs[['raw']])))))	#read pos in raw contig
					y		<- sub('-*$','',paste(as.vector(as.character(crs[['raw']][TAXON_raw, seq.int(tmp, ncol(crs[['raw']])) ])),collapse=''))
					#	y starts just a few sites before x should fit, need to determine exact position
				#print(CNTG)
				#print(x)
				#print(y)
				#print(tmp)				
					z		<- gsub('-','',x)	#rm internal '-' because cut and raw are not in same alignment
					z2		<- gsub('-','',y)	#rm internal '-' because cut and raw are not in same alignment
					z3		<- regexpr( substr(z,1,min(9,nchar(z))),  z2 )
				#print(z3)					
					if(z3>0)
					{
						z2	<- substring(z2, z3)
						z2	<- haircut.large.regexpr(z, z2)
					}
					if(	z3<0 || z3> (abs(diff(c(ncol(crs[['cut']]), ncol(crs[['raw']]))))+1)  )	#positioned raw to fit cut within alignment accuracy - reject match if that s not the case
					{
						z3	<- -1L					
					}					
				#print(z2)	
					list(EQ_FIRSTtmp=as.integer(z3), EQ_LASTtmp=as.integer(z3+attr(z2, 'match.length')-1L+nchar(x)-nchar(z)) )					
				}, by=c('CNTG','OCNTG')]
		txe		<- subset(txe, EQ_FIRSTtmp>0)
		tx		<- merge(tx,txe,all.x=TRUE,by=c('CNTG','OCNTG'))
		tmp		<- tx[, which(!is.na(EQ_FIRSTtmp) & EQ_FIRSTtmp>0)]
		set(tx, tmp, 'EQ_FIRST', tx[tmp,FIRST]+tx[tmp,EQ_FIRSTtmp]-1L)
		set(tx, tmp, 'EQ_LAST', tx[tmp,FIRST]+tx[tmp,EQ_LASTtmp]-1L-tx[tmp,EQ_FIRSTtmp])
		set(tx, NULL, c('EQ_FIRSTtmp','EQ_LASTtmp'),NULL)
		set(tx, tx[, which( CUT=='raw' & OCNTG%in%txe[, unique(OCNTG)] )], 'EQ', TRUE)
		set(tx, tx[, which( !is.na(EQ_FIRST))], 'EQ', TRUE)
	} 					
	tx		<- merge(tx, data.table(INFILE=files, CUT=blastncut), by='CUT')
	subset(tx, select=c(INFILE, TAXON, CUT, CNTG, OCNTG, CCNTG, FIRST, LAST, EQ, EQ_FIRST, EQ_LAST))
}
##--------------------------------------------------------------------------------------------------------
##	deal with large expressions by performing them subsequently
##--------------------------------------------------------------------------------------------------------
haircut.large.regexpr<- function(pattern, x, n=500, ...)
{
	tmp		<- seq_len(ceiling( nchar(pattern)/500 ))*500
	tmp		<- matrix( c(1, tmp[-length(tmp)]+1, tmp[-length(tmp)], nchar(pattern), rep(NA, 2*length(tmp))), byrow=T, nrow=4  )
	for(j in seq_len(ncol(tmp)))
	{
		z			<- regexpr(substr(pattern, tmp[1,j], tmp[2,j]), x, ...)
		tmp[4,j]	<- attr(z,'match.length')
		tmp[3,j]	<- as.numeric(z)
	}
	tmp		<- tmp[, as.logical(cummin( tmp[4,]>0 )), drop=FALSE]
	z		<- all( tmp[3, ]+tmp[4, ]-1==tmp[2, ] )		#check if complete match
	if(ncol(tmp) & z)
	{
		ans							<- tmp[3,1]
		attr(ans,'match.length')	<- sum(tmp[4,])
	}
	if(!ncol(tmp) | !z)
	{
		ans							<- -1
		attr(ans,'match.length')	<- -1
	}
	ans	
}
##--------------------------------------------------------------------------------------------------------
##	return EQ= 1 if raw and contigs are identical in that: 
##	- 	there is only once cut contig which disagrees on up x positions, where x is the difference in alignment lengths in the cut and raw files
##	- 	the concatenated cut contigs disagree on up x positions, where x is the difference in alignment lengths in the cut and raw files
##		gaps between cut contigs are not counted
##--------------------------------------------------------------------------------------------------------
haircut.get.identical.among.raw.and.cut.contigs <- function(indir, files, png_id, blastncut)
{
	crs		<- lapply(files, function(x)
			{
				cr		<- read.dna(file=paste(indir,'/',x,sep=''), format='fasta')
				if(is.matrix(cr))
				{
					cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
					tmp		<- strsplit(basename(x), '_')[[1]][1]
					tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )
					cr		<- cr[ subset(tx, CONTIG==1)[, TAXON], ]	
				}
				if(!is.matrix(cr))
					cr		<- as.DNAbin(matrix(vector('character',0), nrow=2, byrow=T, dimnames=list(c('DUMMY','DUMMY'),c())))
				cr					
			})
	names(crs)	<- blastncut
	#	get contig table
	tx		<- do.call('rbind',lapply(seq_along(crs), function(i)	data.table(	TAXON=rownames(crs[[i]]), 
								CUT= blastncut[i], 
								FIRST= apply( as.character(crs[[i]]), 1, function(x) which(x!='-')[1] ),
								LAST= ncol(crs[[i]])-apply( as.character(crs[[i]]), 1, function(x) which(rev(x)!='-')[1] ),
								CRS_ID=i)	))
	tx		<- subset(tx, !is.na(FIRST) & !is.na(LAST))	#some contigs may just be in LTR
	tx[, CNTG:=tx[, gsub(paste(png_id,'.',sep=''),'',substring(TAXON, regexpr(png_id, TAXON)))]]
	tx[, OCNTG:= tx[, sapply(strsplit(CNTG,'.',fixed=T),'[[',1)]]
	tx[, CCNTG:= NA_character_]		
	tx[, EQ:=FALSE]
	tmp		<- tx[, which(grepl('.',CNTG,fixed=T))]
	if(length(tmp))
		set(tx, tmp, 'CCNTG', tx[tmp, sapply(strsplit(CNTG,'.',fixed=T),'[[',2)])	
	tmp		<- subset(tx, CUT=='cut' & !is.na(CCNTG))[, list(CCNTGn=length(CCNTG)), by='OCNTG']
	tmp		<- subset(tmp, CCNTGn==1)[, OCNTG]	#check for cut contigs that should be present in multiple cuts by naming scheme, but after LTR removal there is only one cut
	if(length(tmp))
	{
		cat('\nFound lone cuts for which multiple cuts are expected by naming scheme, n=', length(tmp))
		tmp	<- tx[, which(CUT=='cut' & OCNTG%in%tmp)]
		set(tx, tmp, 'CCNTG', NA_character_)
		set(tx, tmp, 'CNTG', tx[tmp, OCNTG])
	}
	#	check if there are contigs with corresponding name in 'cut' and that are identical
	#	differences in gaps are allowed: these will come from different alignment
	txe		<- dcast.data.table(tx, CNTG~CUT, value.var='TAXON')
	txe		<- subset( txe, !is.na(cut) & !is.na(raw) )
	if(nrow(txe) && c('cut','raw')%in%colnames(txe))
	{						
		txe		<- txe[, {
					x<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['cut']][cut,])),collapse='')))
					y<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][raw,])),collapse='')))
					if(nchar(x)==nchar(y))
						z	<- x==y
					if(nchar(x)!=nchar(y))
						z	<- gsub('-','',x)==gsub('-','',y)
					list(EQ=z)					
				}, by='CNTG']
		set(tx, which( tx[, CNTG]%in%subset(txe, EQ)[, CNTG] ), 'EQ', TRUE)						
	}
	#	concatenate cut contigs 
	txe		<- subset(tx, !is.na(CCNTG))
	if(nrow(txe))
	{
		txe		<- txe[, {
					z	<- sub('^-*','',paste(as.vector(as.character(crs[['cut']][TAXON[1],])),collapse=''))
					tmp	<- ncol(crs[['cut']])-nchar(z)	#number of initial gaps removed
					#print(tmp)
					z	<- sub('-*$','',z)
					tmp	<- tmp+nchar(z)					#number of sites to remove from next contig
					#print(tmp)
					#print(seq_len(length(TAXON))[-1])
					for(k in seq_len(length(TAXON))[-1])
					{
						if( tmp+1 < ncol(crs[['cut']])	)	#in some cases eg 15065_1_10, reversed contigs and non-reversed contigs are prodived. these cannot be concatenated to reconstitute the original raw contig
						{
							z2	<- paste(as.vector(as.character(crs[['cut']][TAXON[k], seq.int(tmp+1, ncol(crs[['cut']]))])), collapse='')
							z2	<- sub('-*$','',z2)
							tmp	<- tmp+nchar(z2)	
							#print(tmp)
							z	<- paste(z, z2, sep='')	
						}										
					}
					#print(nchar(z))
					list( CSEQ=z )	
				}, by='OCNTG']
		#	check if concatenated contig equals raw contig
		txe		<- subset(merge(tx, txe, by='OCNTG'), CUT=='raw')
		txe		<- txe[, {
					x<- sub('^-*','',sub('-*$','',paste(as.vector(as.character(crs[['raw']][TAXON,])),collapse='')))
					z	<- FALSE
					if(nchar(x)==nchar(CSEQ))
						z	<- x==CSEQ
					#print(abs(diff(c(nchar(x),nchar(CSEQ))))<=abs(diff(sapply(crs,ncol))))
					#print( gsub('-','',x) )
					#print( gsub('-','',CSEQ) )
					if(	abs(diff(c(nchar(x),nchar(CSEQ))))<=abs(diff(sapply(crs,ncol)))	)	#contigs may have cut and been put elsewhere, so only allow for a difference in # gaps that is not larger than the difference in alignment length
						z	<- gsub('-','',x)==gsub('-','',CSEQ)
					list(EQ=z)					
				}, by='OCNTG']
		set(tx, which( tx[, OCNTG]%in%subset(txe, EQ)[, OCNTG] ), 'EQ', TRUE)						
	}					
	tx		<- merge(tx, data.table(INFILE=files, CUT=blastncut), by='CUT')
	tx
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
	invisible(infiles[, {
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
				options(show.error.messages = TRUE)	
				list(	DONE=!inherits(readAttempt, "try-error")	)			
			}, by='FILE'], by='FILE')
	cat(paste('\nFound processed files, n=', infiles[, length(which(DONE))]))
	infiles		<- subset(infiles, !DONE)
	#
	#	infiles[, which(grepl('12559_1_5_cut',FILE))]	fls<- 41
	#	process files
	for(fls in infiles[, seq_along(FILE)])
	{
		file	<- paste(indir, infiles[fls, FILE], sep='/')
		cat(paste('\nProcess', file))
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
		#	get rolling CNS_GPS
		tmp		<- subset(tx, TAXON=='consensus')[, {
					tmp		<- as.character( cnsr.df$CNS_BASE[seq.int(FIRST,LAST)] )=='-'
					list(SITE=seq.int(FIRST,LAST), CNS_GPSr=rollapply( seq_len(LAST-FIRST+1), width=par['GPS.window'], FUN= function(z) mean(tmp[z]), align='center', partial=T ))
				}, by='TAXON']
		cnsr.df	<- merge(cnsr.df, subset(tmp, select=c(SITE,CNS_GPSr)), all.x=1, by='SITE')
		cnsc.df	<- merge(cnsc.df, subset(cnsr.df, select=c(SITE, CNS_FRQr,CNS_GPSr)), all.x=TRUE, by='SITE')
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
##	load training data set for sites 
##--------------------------------------------------------------------------------------------------------
haircut.load.training.data<- function(indir, site)
{
	tmp		<- cut(site, breaks=c(-1,seq.int(200, 10200, 200),Inf), labels=c(paste('<',seq.int(200, 10200, 200),sep=''),'>10200'))
	tmp		<- gsub('<|>','',tmp)
	tmp		<- list.files(indir, pattern=paste('SITES',tmp,'\\.R',sep=''), full.names=T)
	stopifnot(length(tmp)==1)
	load(tmp)
	ctr
}
##--------------------------------------------------------------------------------------------------------
##	get training data set: contig statistics aligned to sites where curated contigs were not cut 
##--------------------------------------------------------------------------------------------------------
haircut.get.training.data<- function(indir, ctrain, par, outdir, outfile)
{
	require(plyr)
	infiles	<- data.table(INFILE=list.files(indir, pattern='\\.R$', recursive=T))
	infiles[, BATCH:= ceiling(seq_len(nrow(infiles))/100)]
	for(b in infiles[, unique(BATCH)])
	{
		#	add curated answers to data by batch
		tmp		<- subset(infiles, BATCH==b)[, {
					cat(paste('\nProcess', INFILE ))
					#	INFILE	<- infiles[12,INFILE]
					load( paste(indir,INFILE,sep='/'))		
					#	select contig+refs for which curated answer available				
					cm		<- merge( unique(subset(cnsc.df, select=c(PNG_ID, TAXON, BLASTnCUT))), subset( ctrain, select=c(PNG_ID, TAXON, BLASTnCUT, USE_IN_TRAIN, EQ_FIRST, EQ_LAST, ANS_FILE, ANS_LEN, ANS_FIRST, ANS_LAST) ), all.x=1, by=c('PNG_ID','TAXON','BLASTnCUT') )
					#	code below identifies ANS_CALL=1
					#	we also need ANS_CALL=0; these are all those that don t have a curated file + all those that have USE_IN_TRAIN=='Y'					
					cm		<- subset(cm, is.na(USE_IN_TRAIN) |  USE_IN_TRAIN=='Y' |  USE_IN_TRAIN=='P')
					if( cm[, all(is.na(ANS_LEN))] )
					{
						cat(paste('\nNo data contigs matched in training data set', INFILE ))
						ca	<- merge(subset(cnsc.df, select=c(TAXON, SITE, AGRpc, GPS, CNS_FRQr, PNG_ID, BLASTnCUT)), subset(cm, select=TAXON), by='TAXON')
						ca[, ANS_CALL:=0L]							
					}
					if( cm[, !all(is.na(ANS_LEN))] && ncol(cnsc)!=subset(cm, !is.na(ANS_LEN))[1, ANS_LEN] )
					{			
						cat(paste('\nMatching contigs exist in training data set, but are not necessarily aligned', INFILE )) 
						cr		<- read.dna(subset(cm, !is.na(ANS_LEN))[1,  ANS_FILE], format='fasta')
						#	determine start of non-LTR position and cut 
						cr		<- cr[, seq.int(haircut.find.nonLTRstart(cr), ncol(cr))]
						#	determine consensus to calculate offset
						tmp		<- strsplit(basename(subset(cm, !is.na(ANS_LEN))[1,  ANS_FILE]), '_')[[1]][1]
						tx		<- data.table(TAXON= rownames(cr), CONTIG=as.integer(grepl(tmp, rownames(cr))) )			
						tmp		<- cr[subset(tx, CONTIG==0)[, TAXON],]
						crc		<- haircut.getconsensus(tmp, par, bases=c('a','c','g','t','-') )$DNAbin	
						cdc		<- cnsc['consensus', ]
						#	calculate site offset in the curated sequence (crc)
						offset	<- haircut.calculate.offset(crc,cdc)
						#	update ANS_FIRST and ANS_LAST according to offset of curated contig in data contigs
						cm		<- cm[, list(ANS_LEN=ANS_LEN+offset[ANS_LEN], ANS_FIRST=ANS_FIRST+offset[ANS_FIRST], ANS_LAST=ANS_LAST+offset[ANS_LAST]), by=c('PNG_ID','TAXON','BLASTnCUT','USE_IN_TRAIN','EQ_FIRST','EQ_LAST')]
						if(any(offset!=0) && ncol(cnsc)!=subset(cm, !is.na(ANS_LEN))[1, ANS_LEN])
							warning(paste('\nPerhaps check: Found unequal lengths for', INFILE))
					}
					if( cm[, !all(is.na(ANS_LEN))] )
					{
						#	expand to SITEs that were called manually
						ca		<- subset(cm, !is.na(ANS_LEN))[, list(SITE=seq.int(ANS_FIRST,ANS_LAST), ANS_CALL=1L), by='TAXON']
						#	keep only those contigs that are not in USE_IN_TRAIN=='N'
						tmp		<- merge(subset(cnsc.df, select=c(TAXON, SITE, AGRpc, GPS, CNS_FRQr, PNG_ID, BLASTnCUT)), subset(cm, select=c(TAXON, USE_IN_TRAIN, EQ_FIRST, EQ_LAST)), by='TAXON')
						ca		<- merge(tmp, ca, all.x=1, by=c('TAXON','SITE'))
						ca		<- subset(ca, USE_IN_TRAIN=='Y' | is.na(USE_IN_TRAIN) | (USE_IN_TRAIN=='P' & SITE>=EQ_FIRST & SITE<=EQ_LAST))
						set(ca, ca[,which(is.na(ANS_CALL))],'ANS_CALL',0L)
						set(ca, NULL, c('USE_IN_TRAIN','EQ_FIRST','EQ_LAST'), NULL)
					}							
					ca
				}, by='INFILE']
		#	save batches in chunks of sites
		cat(paste('\nSave batch to file', b ))
		tmp[, CHUNK:= cut(SITE, breaks=c(-1,seq.int(200, 10200, 200),Inf), labels=c(paste('<',seq.int(200, 10200, 200),sep=''),'>10200' ))]
		for(c in tmp[, unique(CHUNK)])
		{
			ctr	<- subset(tmp, CHUNK==c)				
			save(ctr, file=paste(outdir, '/', outfile,'_sites',gsub('<|>','',c),'_batch',b,'.R',sep=''))
		}		
		cat(paste('\nSaved batch to file', b ))
	}	
	#	load batches and save chunks
	ofiles	<- data.table(FILE=list.files(outdir, pattern='\\.R$', recursive=T))
	set(ofiles, NULL, 'SITES', ofiles[, substring(regmatches(FILE, regexpr('sites[0-9]+',FILE)), 6)])
	ofiles[, {
				ctr	<- do.call('rbind',lapply(FILE, function(x)
								{
									load(paste(outdir,x,sep='/'))
									ctr
								}))
				save(ctr, file=paste(outdir, '/', outfile,'_SITES',SITES,'.R',sep=''))
			}, by='SITES']
	#	rm intermediates files
	invisible( file.remove(ofiles[, paste(outdir,FILE,sep='/')]) )
	#Perhaps check: Found unequal lengths for 15065_1_10_cut_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1367
	#Perhaps check: Found unequal lengths for 15070_1_9_raw_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1522
	#Perhaps check: Found unequal lengths for 15099_1_87_raw_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	1924
	#Perhaps check: Found unequal lengths for 15172_1_43_cut_wRefs_HAIRCUTSTAT_thr5_aw200_fw100_gw200.R	2325
	NULL
}
##--------------------------------------------------------------------------------------------------------
##	get curated contigs: names of curated contigs and first/last sites where they were not cut 
##--------------------------------------------------------------------------------------------------------
haircut.get.curated.contigs<- function(indir, outfile)
{
	options(show.error.messages = FALSE)		
	readAttempt		<-try(suppressWarnings(load(outfile)))
	options(show.error.messages = TRUE)
	tmp				<- !inherits(readAttempt, "try-error")
	if(!tmp)
	{
		infiles	<- data.table(FILE=list.files(indir, recursive=T, pattern='fasta$'))
		infiles[, PNG_ID:= gsub('\\.fasta','',basename(FILE))]
		stopifnot( nrow(infiles)==infiles[, length(unique(PNG_ID))] )					
		ctrain	<- infiles[, {
					file	<- paste(indir, FILE, sep='/')
					cat(paste('\nProcess', file))
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
					cns		<- cr[subset(tx, CONTIG==1)[, TAXON],]
					#	determine first and last non-gap sites
					tx		<- data.table(	TAXON= rownames(cns), PNG_ID=PNG_ID, LEN= ncol(cr),
							FIRST= apply( as.character(cns), 1, function(x) which(x!='-')[1] ),
							LAST= ncol(cns)-apply( as.character(cns), 1, function(x) which(rev(x)!='-')[1] )		)
					subset(tx, !is.na(FIRST) & !is.na(LAST))	#	some contigs only map into LTR
				}, by='FILE']
		set(ctrain, NULL, 'FILE', ctrain[, paste(indir, FILE, sep='/')])
		save(ctrain, file=outfile)
	}
	ctrain
}
##--------------------------------------------------------------------------------------------------------
##	calculate offset of the first sequence to the second sequence, assuming that the only difference are gap characters
##--------------------------------------------------------------------------------------------------------
haircut.calculate.offset<- function(x,y)
{
	x		<- as.character(x)	
	y		<- as.character(y)		
	stopifnot( gsub('-','',paste(as.vector(x), collapse=''))==gsub('-','',paste(as.vector(y), collapse='')) )
	z		<- rbind.fill.matrix(x,y)
	offset	<- rep(0, ncol(z))			
	k		<- seq_len(ncol(z))+offset
	k		<- k[ k>0 & k<=ncol(z)]
	k		<- which( z[1, seq_along(k)]!=z[2, k] )[1]
	while(!is.na(k))
	{
		if( z[1,k]!='-' )
			offset[ seq.int(k,length(offset)) ] <- offset[ seq.int(k,length(offset)) ]+1
		if( z[1,k]=='-' )
			offset[ seq.int(k,length(offset)) ] <- offset[ seq.int(k,length(offset)) ]-1
		k		<- seq_len(ncol(z))+offset
		k		<- k[ k>0 & k<=ncol(z)]
		k		<- which( z[1, seq_along(k)]!=z[2, k] )[1]		
	}
	if(length(offset)>length(x))
		offset	<- offset[ seq_len(length(x)) ]
	offset
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
##	find the last site inside the reference alignment 
##--------------------------------------------------------------------------------------------------------
haircut.find.lastRefSite<- function(cr)
{
	
	ans				<- seq.find.pos.of.pattern(cr, pattern='t-*t-*t-*t-*a-*g-*t-*c-*a-*g-*t-*g-*t-*g-*g-*a-*a-*a-*a-*t-*c-*t-*c-*t-*a-*g-*c-*a', row.name='B.FR.83.HXB2_LA')
	stopifnot(length(ans)==1, ans>0)
	as.integer(ans+attr(ans,'match.length'))
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