######################################################################################
#	return GAG POL ENV ancestral sequences from BEAST PARSER output	
#	olli originally written 06-08-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	bseq		data.table containing original sequences. only needed for BEAST decompression.
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq.withDecompression<- function(tree, node.stat, bseq, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=2)
{
	require(data.table)
	require(ape)
	
	tree.id				<- names(tree)
	#	add calendar time for inner node at NODE_ID to node.stat
	node.stat[, CALENDAR_TIME:=NA_real_]		
	setkey(node.stat, TREE_ID, NODE_ID)
	for(i in seq_along(tree.id))
	{
		label.ctime			<- sapply( strsplit(tree[[i]]$tip.label, label.sep, fixed=TRUE), '[[', label.idx.ctime ) 
		label.ctime			<- as.numeric(label.ctime)			
		depth				<- node.depth.edgelength( tree[[ i ]] )
		tmp					<- which.max(depth)
		depth				<- depth-depth[tmp]+label.ctime[tmp]
		for(j in seq_along(depth))
			set(node.stat, node.stat[, which(TREE_ID==tree.id[i] & NODE_ID==j)], 'CALENDAR_TIME', depth[j])					
	}
	tmp			<- node.stat[, length(which(is.na(CALENDAR_TIME)))]
	cat(paste('\nTotal node statistics with no CALENDAR_TIME [should be zero], n=', tmp  ))
	#	keep only inner nodes
	tmp			<- lapply(seq_along(tree.id), function(i)
			{
				subset(node.stat, TREE_ID==tree.id[i] & NODE_ID>Ntip(tree[[i]]))
			})
	node.stat	<- do.call('rbind',tmp)
	set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
	#
	#	reconstruct ancestral sequences, need to decompress patterns that were compressed with BEAST
	#	TODO ** this results in duplicate columns and should be removed at a later point **
	#
	bseq			<- merge(bseq, bseq[, list(SEQ_N=nchar(SEQ)), by=c('GENE','TAXON_ID')], by=c('GENE','TAXON_ID'))
	bseq			<- bseq[, {
				tmp<- unlist(strsplit(SEQ,''))
				list(	CP1= paste(tmp[seq.int(1,length(tmp),by=3)], collapse='',sep=''), 
						CP2= paste(tmp[seq.int(2,length(tmp),by=3)], collapse='',sep=''), 
						CP3= paste(tmp[seq.int(3,length(tmp),by=3)], collapse='',sep='') 	)
			}, by=c('TAXON_ID','GENE')]		
	bseq			<- melt(bseq, measure.var=c('CP1','CP2','CP3'), variable.name='CODON_POS', value.name='SEQ')
	#	get index of orginal patterns and duplicate patterns
	bseq.decompress	<- bseq[, {
				#print(GENE)
				#print(CODON_POS)
				tmp		<- do.call('rbind',strsplit(SEQ,''))
				tmp2	<- apply( tmp, 2, function(x) paste(x,sep='',collapse=''))	#identical patterns?
				tmp2	<- which(duplicated(tmp2))
				#for each duplicate, work out index of original
				tmp3	<- sapply(tmp2, function(j1) which(apply( tmp[,seq.int(1,j1-1,1), drop=FALSE] == tmp[,j1], 2, all))[1]	 )
				list(NSEQ=ncol(tmp), DUPLICATE_PATTERN=tmp2, MOTHER_PATTERN=tmp3)
			}, by=c('GENE','CODON_POS')]		
	set(bseq.decompress, bseq.decompress[, which(GENE=='env')], 'GENE', 'ENV')
	set(bseq.decompress, bseq.decompress[, which(GENE=='gag')], 'GENE', 'GAG')
	set(bseq.decompress, bseq.decompress[, which(GENE=='pol')], 'GENE', 'POL')
	set(bseq.decompress, NULL, 'xSTAT', bseq.decompress[, paste(GENE,CODON_POS,sep='.')])		
	#	reconstruct ancestral genes sequences - decompress patterns		
	ancseq	<- node.stat[,  {													
				#print(STAT)
				#TREE_ID<- 'STATE_0'
				#STAT<- 'GAG.CP3'
				tmp								<- subset(bseq.decompress, xSTAT==STAT)											
				seq								<- matrix(data=NA_character_, nrow=length(VALUE), ncol=tmp[,NSEQ])
				seq.compressed					<- setdiff( seq_len(ncol(seq)), tmp[, DUPLICATE_PATTERN])	
				#print(dim(seq))										
				tmp2							<- do.call('rbind',strsplit(VALUE,''))
				#print(dim(tmp2))
				stopifnot(length(seq.compressed)==ncol(tmp2))
				seq[, seq.compressed]			<- tmp2
				#print(seq[1,])
				seq[, tmp[, DUPLICATE_PATTERN]] <- seq[, tmp[, MOTHER_PATTERN]]
				#print(seq[1,])
				#stop()
				seq								<- apply(seq, 1, function(x) paste(x, sep='',collapse=''))
				list(TREE_ID=TREE_ID, NODE_ID=NODE_ID, CALENDAR_TIME=CALENDAR_TIME, SEQ=seq) 
			}, by=c('STAT')]
	#	checks of ancseq before we proceed
	tmp		<- ancseq[, list(NSEQ= nchar(SEQ)), by=c('TREE_ID', 'NODE_ID', 'STAT')]		
	stopifnot( tmp[, list(CHECK= all(NSEQ==NSEQ[1])), by='STAT'][, all(CHECK)] )
	set(tmp, NULL, 'GENE', tmp[, sapply(strsplit(STAT,'\\.'),'[[',1)])
	set(tmp, NULL, 'CODON_POS', tmp[, sapply(strsplit(STAT,'\\.'),'[[',2)])
	stopifnot( tmp[, list(CHECK=all(NSEQ==NSEQ[1])), by='GENE'][, all(CHECK)] )
	#	reconstruct genes from codon positions
	ancseq		<- dcast.data.table(ancseq, TREE_ID + NODE_ID + CALENDAR_TIME ~ STAT, value.var="SEQ")
	ancseq		<- ancseq[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, CALENDAR_TIME=CALENDAR_TIME)
			}, by=c('TREE_ID','NODE_ID')]
	set(ancseq, NULL, 'LABEL', ancseq[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(ancseq, NULL, 'BEAST_MCMC_IT', ancseq[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	ancseq		<- subset(ancseq, BEAST_MCMC_IT>tree.id.burnin)
	#
	#	return DNAbin
	#
	ancseq.gag				<- tolower(do.call('rbind',strsplit(ancseq[, GAG],'')))
	rownames(ancseq.gag)	<- ancseq[, LABEL]
	ancseq.gag				<- as.DNAbin(ancseq.gag)		
	ancseq.pol				<- tolower(do.call('rbind',strsplit(ancseq[, POL],'')))
	rownames(ancseq.pol)	<- ancseq[, LABEL]
	ancseq.pol				<- as.DNAbin(ancseq.pol)		
	ancseq.env				<- tolower(do.call('rbind',strsplit(ancseq[, ENV],'')))
	rownames(ancseq.env)	<- ancseq[, LABEL]
	ancseq.env				<- as.DNAbin(ancseq.env)				
	#ancseq					<- cbind(ancseq.gag, ancseq.pol, ancseq.env)
	#
	list(GAG=ancseq.gag, POL=ancseq.pol, ENV=ancseq.env)
}
######################################################################################
#	return GAG POL ENV ancestral sequences from BEAST PARSER output	
#	olli originally written 13-08-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq<- function(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=5)
{
	require(data.table)
	require(ape)
	
	tree.id				<- names(tree)
	#	add calendar time for inner node at NODE_ID to node.stat
	node.stat[, CALENDAR_TIME:=NA_real_]		
	setkey(node.stat, TREE_ID, NODE_ID)
	for(i in seq_along(tree.id))
	{
		cat(paste('\nProcess CALENDAR_TIME for tree id', tree.id[i]  ))
		label.ctime			<- sapply( strsplit(tree[[i]]$tip.label, label.sep, fixed=TRUE), '[[', label.idx.ctime ) 
		label.ctime			<- as.numeric(label.ctime)			
		depth				<- node.depth.edgelength( tree[[ i ]] )
		tmp					<- which.max(depth)
		depth				<- depth-depth[tmp]+label.ctime[tmp]
		tmp					<- node.stat[, which(TREE_ID==tree.id[i])]
		for(j in seq_along(depth))
		{
			set(node.stat, tmp[ node.stat[tmp, which(NODE_ID==j)] ], 'CALENDAR_TIME', depth[j])
		}								
	}
	tmp			<- node.stat[, length(which(is.na(CALENDAR_TIME)))]
	cat(paste('\nTotal node statistics with no CALENDAR_TIME [should be zero], n=', tmp  ))
	#	keep only inner nodes
	tmp			<- sapply(tree, Ntip)
	stopifnot(all(tmp==tmp[1]))
	node.stat	<- subset(node.stat, NODE_ID>tmp[1])	
	#
	set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
	#	checks of ancseq before we proceed
	tmp			<- node.stat[, list(NSEQ= nchar(VALUE)), by=c('TREE_ID', 'NODE_ID', 'STAT')]		
	stopifnot( tmp[, list(CHECK= all(NSEQ==NSEQ[1])), by='STAT'][, all(CHECK)] )
	set(tmp, NULL, 'GENE', tmp[, sapply(strsplit(STAT,'\\.'),'[[',1)])
	set(tmp, NULL, 'CODON_POS', tmp[, sapply(strsplit(STAT,'\\.'),'[[',2)])	
	tmp			<- dcast.data.table(tmp, TREE_ID + NODE_ID + GENE ~ CODON_POS, value.var='NSEQ')
	tmp			<- tmp[, list(CPM=min(CP1, CP2, CP3)), by=c('TREE_ID','NODE_ID','GENE')]
	stopifnot( tmp[, list(CHECK=all(CPM==CPM[1])), by='GENE'][, all(CHECK)] )
	setkey(tmp, GENE)
	#	truncate to following size of coding regions (if necessary)
	tmp			<- unique(tmp)[, list(STAT=paste(GENE,'.CP',1:3,sep=''), CPM=CPM), by='GENE']
	node.stat	<- merge(node.stat, subset(tmp, select=c(STAT, CPM)), by='STAT')
	set(node.stat, NULL, 'VALUE', node.stat[, substr(VALUE,1,CPM)])
	set(node.stat, NULL, 'CPM', NULL)
	#	reconstruct genes from codon positions
	node.stat	<- dcast.data.table(node.stat, TREE_ID + NODE_ID + CALENDAR_TIME ~ STAT, value.var="VALUE")
	node.stat	<- node.stat[, {
				tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
				env		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(GAG.CP1,GAG.CP2,GAG.CP3), strsplit, ''))
				gag		<- paste(as.vector(tmp), collapse='')
				tmp		<- do.call('rbind',sapply(list(POL.CP1,POL.CP2,POL.CP3), strsplit, ''))
				pol		<- paste(as.vector(tmp), collapse='')
				list(GAG=gag, POL=pol, ENV=env, CALENDAR_TIME=CALENDAR_TIME)
			}, by=c('TREE_ID','NODE_ID')]
	
	set(node.stat, NULL, 'LABEL', node.stat[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(node.stat, NULL, 'BEAST_MCMC_IT', node.stat[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	node.stat			<- subset(node.stat, BEAST_MCMC_IT>tree.id.burnin)
	cat(paste('\nFound ancestral sequences, n=', nrow(node.stat)  ))
	node.stat
}