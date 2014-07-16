


prog.hello<- function()	
{
	print('hello')
}
##--------------------------------------------------------------------------------------------------------
##	simple first program to generate files for Matt s phylo simulator
##--------------------------------------------------------------------------------------------------------
prog.HPTN071.parser.v1<- function()	
{
	require(data.table)
	fin.ind		<- '/Users/Oliver/git/HPTN071sim/raw/140716_RUN001_IND.csv'
	fin.trm		<- '/Users/Oliver/git/HPTN071sim/raw/140716_RUN001_TRM.txt'
	fout.ind	<- '/Users/Oliver/git/HPTN071sim/sim/140716_RUN001_IND.csv'
	fout.trm	<- '/Users/Oliver/git/HPTN071sim/sim/140716_RUN001_TRM.csv'
	
	
	setup.df<- data.table(stat= c('yr.start','yr.end','s.INC.recent','s.INC.recent.len', 's.PREV.min', 's.PREV.max'), v=c(1980, 2020, 0.1, 5, 0.01, 0.25) )
	setkey(setup.df, stat)	
	
	df.trm	<- as.data.table(read.csv(fin.trm, stringsAsFactors=FALSE, sep=' ', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))
	set(df.trm, NULL, 'YR', df.trm[, floor(TIME_TR)])	
	
	df.ind	<- as.data.table(read.csv(fin.ind, stringsAsFactors=FALSE))		
	setnames(df.ind, c("Id","Gender","DoB","DateOfDeath","RiskGroup","Circumcised"), c('IDPOP','GENDER','DOB','DOD','RISK','CIRCM'))
	set(df.ind, df.ind[, which(CIRCM=='')], 'CIRCM', NA_character_)
	set(df.ind, NULL, 'CIRCM', df.ind[, factor(CIRCM)])
	set(df.ind, NULL, 'GENDER', df.ind[, factor(GENDER)])
	set(df.ind, NULL, 'RISK', df.ind[, factor(RISK)])	
	set(df.ind, df.ind[, which(DOD==-1)], 'DOD', setup.df['yr.end',][,v]+1.)		
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
	df.sample	<- subset( df.epi, YR>= setup.df['yr.start',][,v] & YR<setup.df['yr.end',][,v] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( 1+setup.df['s.PREV.max',][,v]-setup.df['s.PREV.min',][,v] ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, exp( tmp*(YR-min(YR)) ) - 1 + setup.df['s.PREV.min',][,v] ]
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
									TIME_SEQ=runif(length(tmp),min=TIME_TR[tmp],max=TIME_TR[tmp]), 
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
	stopifnot( all( tmp[,s.n.TOTAL]==df.sample[, s.n.TOTAL] ) )
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
	tmp			<- subset(df.trms, IDTR>=0, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	#
	#	PLOTS
	#
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
			scale_x_continuous(name='year', breaks=seq(1980,2020,2)) + scale_y_continuous(name='total')	+
			facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
	file<- paste(substr(fin.ind, 1, nchar(fin.ind)-7),'INFO_Totals.pdf',sep='')
	ggsave(file=file, w=16, h=8)
	#	plot distribution between transmission time and sequencing time
	tmp	<- subset(df.inds, !is.na(TIME_SEQ))
	set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
	ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
			scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
	file<- paste(substr(fin.ind, 1, nchar(fin.ind)-7),'INFO_Time2Seq.pdf',sep='')
	ggsave(file=file, w=8, h=8)
	#	plot transmission network
	file		<- paste(substr(fin.ind, 1, nchar(fin.ind)-7),'INFO_TrNetworks.pdf',sep='')
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
	#
	#	SAVE SAMPLED RECIPIENTS AND TRANSMISSIONS TO SAMPLED RECIPIENTS
	#
	#	save for us
	file		<- paste(substr(fin.ind, 1, nchar(fin.ind)-7),'SAVE.R',sep='')
	save(file=file, df.epi, df.trms, df.inds, df.sample)
	#	save for Matt
	#	exclude columns that are not needed	
	df.inds	<- subset(df.inds, !is.na(TIME_TR))
	df.inds[, RISK:=NULL]
	df.inds[, INCIDENT_SEQ:=NULL]
	df.inds[, TIME_SEQYR:=NULL]	
	df.trms[, TR_ACUTE:=NULL]
	df.trms[, YR:=NULL]	
	cat(paste('\nwrite to file',fout.ind))
	write.csv(file=fout.ind, df.inds)
	cat(paste('\nwrite to file',fout.trm))
	write.csv(file=fout.trm, df.trms)
	
}

