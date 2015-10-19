##--------------------------------------------------------------------------------------------------------
##	olli 05.10.15
##--------------------------------------------------------------------------------------------------------
treecomparison.submissions.161015<- function()	
{
	require(data.table)
	require(ape)
	require(phangorn)
	#
	#	get true trees
	#
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/freeze_July15'
	tfiles	<- list.files(indir, pattern='newick$', full.names=TRUE)
	tfiles	<- data.table( FILE_T=tfiles[ grepl('SUBSTTREE', tfiles) | grepl('Vill_99', tfiles) | grepl('Vill.*DATEDTREE', tfiles) ] )
	tfiles[, SC:= toupper(gsub('_SUBSTTREE|_DATEDTREE','',gsub('.newick','',basename(FILE_T))))]
	tmp		<- rbind( subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15'), subset(tfiles, SC=='VILL_99_APR15') )
	set(tmp, NULL, 'SC', c('150701_VILL_SCENARIO-C','150701_VILL_SCENARIO-D','150701_VILL_SCENARIO-E'))
	tfiles	<- rbind(tfiles, tmp)
	ttrs	<- lapply(tfiles[, FILE_T], function(x)	read.tree(file=x) )
	names(ttrs)	<- tfiles[, SC]
	tfiles[, IDX_T:=seq_along(ttrs)]
	tfiles[, TAXAN_T:= sapply(ttrs, Ntip)]
	#	info on true trees
	tinfo	<- merge(tfiles, do.call('rbind',lapply(seq_along(ttrs), function(i) data.table(TAXA=ttrs[[i]]$tip.label, IDX_T=i))), by='IDX_T')	
	tinfo[, IDPOP:=NA_character_]
	tmp		<- tinfo[, which(grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp,regmatches(TAXA, regexpr('IDPOP_[0-9]+',TAXA))])
	tmp		<- tinfo[, which(!grepl('REGIONAL',SC))]
	set(tinfo, tmp, 'IDPOP', tinfo[tmp, regmatches(TAXA, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',TAXA))])		
	stopifnot(subset(tinfo, grepl('VILL',SC))[, length(which(substring(TAXA,1,10)!=substring(IDPOP,1,10)))]==0)	
	stopifnot( tinfo[, length(which(is.na(IDPOP)))==0] )	
	set(tinfo, NULL, 'IDPOP', tinfo[,toupper(IDPOP)])
	set(tinfo, NULL, 'TAXA', tinfo[,toupper(TAXA)])
	#	read cluster membership from DATEDCLUTREES	
	tmp		<- list.files(indir, pattern='DATEDCLUTREES', full.names=TRUE)
	tmp		<- data.table( FILE_CLU_T= tmp, SC= toupper(gsub('_DATEDCLUTREES','',gsub('.newick','',basename(tmp))))) 
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)	
	tmp		<- subset(tfiles, !is.na(FILE_CLU_T))[, {
				z		<- read.tree(FILE_CLU_T)
				do.call('rbind',lapply(seq_along(z), function(i) data.table(IDCLU=i, TAXA=z[[i]]$tip.label)))				
			}, by='SC']	
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all=1)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N= length(IDPOP)), by=c('SC','IDCLU')]
	tinfo	<- merge(tinfo, tmp, by=c('SC','IDCLU'), all=1)
	#
	#	get submitted trees
	#	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))
	infiles	<- c(infiles, list.files(indir, pattern="best_tree.newick", recursive=1, full.names=1))
	infiles	<- data.table(FILE=infiles)
	strs	<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(strs)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) paste(names(x)[1],'_use',sep='')), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	names(MetaPIGA.trees)	<- gsub("'",'',names(MetaPIGA.trees), fixed=1)	
	strs					<- c(strs, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(strs))
	#
	#
	#	
	submitted.info[, IDX:=seq_along(strs)]	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML|RAxML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA|Consensus pruning|Best individual of population',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')	
	stopifnot( submitted.info[, length(which(is.na(TEAM)))==0] )
	#
	#	scenario
	#	
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	stopifnot( submitted.info[, length(which(is.na(SC)))] )
	#
	#	set covariates of scenarios
	#
	tmp		<- data.table(	SC=		c("150701_REGIONAL_TRAIN1","150701_REGIONAL_TRAIN2","150701_REGIONAL_TRAIN3","150701_REGIONAL_TRAIN4" ,"150701_REGIONAL_TRAIN5", "150701_VILL_SCENARIO-A", "150701_VILL_SCENARIO-B", "VILL_99_APR15","150701_VILL_SCENARIO-C", "150701_VILL_SCENARIO-D", "150701_VILL_SCENARIO-E"),
				MODEL=	c('R','R','R','R','R','V','V','V','V','V','V'),
				SEQCOV= c(0.16, 0.16, 0.16, 0.16, 0.16, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
				ACUTE=	c('low', 'low', 'high', 'low', 'high', 'high', 'high', 'high', 'high', 'high', 'high'),
				GAPS=	c('none', 'low', 'low', 'high', 'high', 'low', 'high', 'none', 'none', 'low', 'high'), 
				ART=	c('none', 'none', 'none', 'none', 'none', 'none', 'none', 'fast', 'fast', 'fast', 'fast'),
				EXT= 	c('5pc', '5pc', '5pc', '5pc', '5pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc', '~0pc')
				)
	submitted.info	<- merge(submitted.info, tmp, by='SC')
	#
	#	best tree for each scenario
	#
	submitted.info[, BEST:='N']
	set(submitted.info, submitted.info[, which(grepl('RAxML', FILE) & grepl('best_tree', FILE))], 'BEST', 'Y')	
	#	copied from ListOfBestTrees_IQTree150818.txt
	#	there are several best trees for some scenarios
	tmp	<- c( '150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_07',	
	'150701_Vill_SCENARIO-A_IQTree150814_partition_12_3_04.',	
	'150701_Vill_SCENARIO-B_IQTree150814_partition_12_3_03.',	
	'Vill_99_Apr15_IQTree150814_partition_123.',			
	'150701_Vill_SCENARIO-D_IQTree150814_partition_12_3.',	
	'150701_Vill_SCENARIO-E_IQTree150814_partition_12_3.',
	'150701_Vill_SCENARIO-A_IQTree150814_pol_partition_12_3.',	
	'150701_Vill_SCENARIO-B_IQTree150814_pol_partition_12_3_05.',
	'Vill_99_Apr15_IQTree150814_pol_partition_12_3_09.',		
	'Vill_99_Apr15_IQTree150814_pol_partition_12_3_10.',		
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_05.',
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_06.',
	'150701_Vill_SCENARIO-D_IQTree150814_pol_partition_12_3_09.',
	'150701_Vill_SCENARIO-E_IQTree150814_pol_partition_12_3_06.',
	'150701_Regional_TRAIN1_IQTree150818_partition_123_03.',	
	'150701_Regional_TRAIN2_IQTree150814_partition_123_03.',	
	'150701_Regional_TRAIN3_IQTree150814_partition_123_01.',	
	'150701_Regional_TRAIN4_IQTree150814_partition_123_02.',	
	'150701_Regional_TRAIN5_IQTree150814_partition_123.',
	'150701_Regional_TRAIN1_IQTree150818_pol_partition_123_05.',
	'150701_Regional_TRAIN2_IQTree150814_pol_partition_123_10.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_05.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_06.',
	'150701_Regional_TRAIN3_IQTree150814_pol_partition_123_08.',
	'150701_Regional_TRAIN4_IQTree150814_pol_partition_123_10.',
	'150701_Regional_TRAIN5_IQTree150814_pol_partition_123_05.')
	tmp	<- sapply(tmp, function(x) submitted.info[, which((grepl('IQTree150814/', FILE, fixed=1) | grepl('IQTree150818/', FILE, fixed=1)) & grepl(x, FILE, fixed=1))] )
	set(submitted.info, tmp, 'BEST', 'Y')
	#	PhyML no replicates: all files best
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'BEST', 'Y')		
	#
	#	set OTHER (ie old or some preliminary/unknown tree)
	#
	submitted.info[, OTHER:='N']
	#	MetaPIGA tree to be used is first in nexus list (which was tagged with best above)
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA' & !grepl('use', FILE))], 'OTHER', 'Y')
	#	IQTree did several uploads, use only most recent in main analysis
	set(submitted.info, submitted.info[, which(grepl('150701_Regional_TRAIN1_IQTree150814', FILE))], 'OTHER', 'Y')
	#
	#	set which gene used to construct tree (either pol or concatenated gag+pol+env)
	#
	submitted.info[, GENE:=NA_character_]
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('full', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='RAXML' & grepl('pol', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='RAXML' & is.na(GENE)))==0)
	set(submitted.info, submitted.info[, which(TEAM=='PhyML')], 'GENE', 'POL')
	set(submitted.info, submitted.info[, which(TEAM=='MetaPIGA')], 'GENE', 'GAG+POL+ENV')	
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_partition', FILE))], 'GENE', 'GAG+POL+ENV')
	set(submitted.info, submitted.info[, which(TEAM=='IQTree' & grepl('[0-9]_pol_partition', FILE))], 'GENE', 'POL')
	stopifnot(nrow(subset(submitted.info, TEAM=='IQTree' & is.na(GENE)))==0)
	#
	#	number taxa in tree
	#
	setkey(submitted.info, IDX)
	submitted.info[, TAXAN:= sapply(strs, Ntip)]
	#
	#	are trees rooted?
	#
	setkey(submitted.info, IDX)
	submitted.info[, ROOTED:=factor(sapply(strs, is.rooted),levels=c(TRUE,FALSE),labels=c('Y','N'))]
	#
	#	add index of true tree
	#
	require(phangorn)
	submitted.info	<- merge(submitted.info, subset(tfiles, select=c('SC','IDX_T','TAXAN_T')), by='SC')
	stopifnot(nrow(subset(submitted.info, TAXAN>TAXAN_T))==0)
	#
	#	fix taxa names that teams have changed
	#
	tmp		<- subset(submitted.info, TEAM=='IQTree' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		strs[[i]]$tip.label	<- sapply(strsplit(strs[[i]]$tip.label,'_'), function(x)	paste(x[1],'_',x[2],'|',x[3],'|',x[4],'_',x[5],'|',x[6],sep='')	)
	}
	for(i in seq_along(strs))
	{
		strs[[i]]$tip.label	<- toupper(strs[[i]]$tip.label)
	}
	for(i in seq_along(ttrs))
	{
		ttrs[[i]]$tip.label	<- toupper(ttrs[[i]]$tip.label)
	}	
	###
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='R')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('IDPOP_[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	tmp		<- subset(submitted.info, TEAM=='PhyML' & MODEL=='V')[, IDX]
	for(i in tmp)
	{
		
		z	<- data.table(IDX=seq_along(strs[[i]]$tip.label), IDPOP=regmatches(strs[[i]]$tip.label, regexpr('HOUSE[0-9]+-[0-9]+|House[0-9]+-[0-9]+',strs[[i]]$tip.label)), SC=subset(submitted.info, IDX==i)[,SC])
		z	<- merge(subset(tinfo, select=c(IDPOP,SC,TAXA)), z, by=c('IDPOP','SC'))
		stopifnot(nrow(z)==length(strs[[i]]$tip.label))
		setkey(z, IDX)
		strs[[i]]$tip.label	<- z[, TAXA]
	}
	#
	#	compute Robinson Fould of complete tree
	#
	setkey(submitted.info, IDX)
	#tmp				<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 463
	#IDX_T<-6
	tmp				<- submitted.info[, {
				#cat('\nAt IDX', IDX)
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(multi2di(ttrs[[IDX_T]]))				
				if(!is.binary.tree(stree))
				{
					cat('\nFound non-binary tree at IDX',IDX)
					stree	<- multi2di(stree)
				}
				#print(stree)
				#print(otree)
				z			<- setdiff(otree$tip.label, stree$tip.label)
				stopifnot( length(z)==abs(diff(c(TAXAN, TAXAN_T))) )
				if(length(z))
					otree	<- unroot(drop.tip(otree, z))				
				#https://groups.google.com/forum/#!topic/raxml/JgvxgknTeqw
				#normalize with 2n-6		
				rf			<- RF.dist(otree, stree, check.labels=TRUE)
				list(RF=rf, NRF=rf/(2*Ntip(otree)-6))
			}, by='IDX']
	submitted.info	<- merge(submitted.info, tmp, by='IDX')
	#
	#	compute Robinson Fould of clusters, then take sum
	#
	tmp		<- subset(submitted.info, MODEL=='Model: Regional')[1,]	
	IDX<- 1
	IDX_T<-1
	setkey(tinfo, IDX_T)
	tmp		<- subset(submitted.info, MODEL=='R')[, {
				cat('\nAt IDX', IDX)
				stree		<- unroot(strs[[IDX]])
				otree		<- unroot(multi2di(ttrs[[IDX_T]]))				
				if(!is.binary.tree(stree))
				{
					cat('\nFound non-binary tree at IDX',IDX)
					stree	<- multi2di(stree)
				}
				z			<- IDX_T
				z			<- subset(tinfo, CLU_N>3 & IDX_T==z)
				z			<- merge(z, data.table(TAXA=stree$tip.label, IN_STREE=1), by='TAXA', all.x=1)
				z			<- merge(z, z[, list(CLU_NS= length(which(IN_STREE==1))), by='IDCLU'], by='IDCLU')
				z			<- subset(z, CLU_NS>3)
				if(nrow(z))
				{
					#IDCLU	<- 6
					#TAXA	<- subset(z, IDCLU==6)[, TAXA]
					ans		<- z[, {								
										sclu	<- unroot(drop.tip(stree, setdiff(stree$tip.label,TAXA)))
										oclu	<- unroot(drop.tip(otree, union( setdiff(otree$tip.label, stree$tip.label), setdiff(otree$tip.label,TAXA))))
										rf		<- RF.dist(oclu, sclu, check.labels=TRUE)
										list(TAXA_NC=Ntip(oclu), RFC=rf, NRFC=rf/(2*Ntip(oclu)-6))
									}, by='IDCLU']	
				}
				if(!nrow(z))
					ans		<- data.table(IDCLU=NA_integer_, TAXA_NC=NA_integer_, RFC=NA_integer_, NRFC=NA_real_)
				ans			
			}, by='IDX']	
	sclu.info	<- merge(submitted.info, tmp, by='IDX')	
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151016.rda'
	save(strs, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
}
treecomparison.ana.151019<- function()
{
	require(ggplot2)
	edir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation'
	file	<- paste(edir,'/','submitted_151016.rda',sep='')
	load(file)
	sa		<- copy(submitted.info)
	set(sa, NULL, 'MODEL', sa[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sa, sa[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sa, NULL, 'SC', sa[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
										labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sa, NULL, 'GAPS', sa[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sa, NULL, 'BEST', sa[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sa, NULL, 'GENE', sa[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sa, NULL, 'TEAM', sa[, factor(TEAM, levels=sa[, sort(unique(TEAM))],labels=sa[, sort(unique(TEAM))])])
	sa		<- subset(sa, OTHER=='N')
	
	#
	#	polvsall by gaps
	#	-->
	#	all leads to improvements throughout
	#	with gaps, the topology without branch lengths is increasingly difficult to estimate
	#	regional overall more difficult!
	ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=NRF, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
			geom_jitter(position = position_jitter(height = .01, width=0.2)) +			
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			scale_fill_brewer(palette='Paired') +
			scale_colour_brewer(palette='Paired') +
			facet_wrap(MODEL~GAPS, scales='free_x') +	
			labs(x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
			theme_bw() 
	ggsave(w=10, h=6, file=paste(edir,'/151016_RF_polvsall_by_gaps.pdf',sep=''))

	
	#	teams
	#	--> 
	#	MetaPIGA fairly bad in terms of RF
	#	PhyML IQTree RAxML similar in terms of RF,
	#	but PhyML behaved poorly when many gaps
	ggplot( sa, aes(y=NRF, x=SC, shape=TEAM, colour=TEAM, fill=TEAM, size=BEST) ) + 
			geom_jitter(position = position_jitter(height = .01, width=0.2), alpha=0.7) +
			scale_size_manual(values=c(4, 1)) +
			scale_shape_manual(values=c(21,22,23,24)) +
			scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
			facet_grid(~GENE) +			
			labs(x='\nsimulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='Method', colour='Method') +
			theme_bw() 
	ggsave(w=10, h=5, file=paste(edir,'/151016_RF_team_by_scenarioandgene.pdf',sep=''))
	
	#	effect of not using all sequences?
	#	effect of acute in terms of RF? --> potentially Yes
	#	effect of ART roll out in terms of RF? --> No
	


}
treecomparison.submissions.300915<- function()	
{
	require(data.table)
	require(ape)
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/IQTree/IQTree201507'
	infiles	<- list.files(indir, pattern='treefile$', recursive=1, full.names=1)
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/PhyML'
	infiles	<- c(infiles, list.files(indir, pattern='*tree*', recursive=1, full.names=1))
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/RAxML'
	infiles	<- c(infiles, list.files(indir, pattern='*RAxML_bestTree*', recursive=1, full.names=1))	
	infiles	<- data.table(FILE=infiles)
	submitted.trees			<- lapply(infiles[, FILE], function(x)
			{
				cat(x)
				read.tree(file=x)	
			})
	names(submitted.trees)	<- infiles[, FILE]
	
	indir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/MetaPIGA'
	tmp		<-  list.files(indir, pattern='*result*', recursive=1, full.names=1)
	tmp		<- data.table(FILE=tmp)
	
	tmp.trees			<- lapply(tmp[, FILE], function(x)
			{
				cat(x)
				read.nexus(file=x)	
			})
	sapply(tmp.trees, length)
	MetaPIGA.trees			<- c(lapply(tmp.trees, '[[', 1), lapply(tmp.trees, '[[', 2), lapply(tmp.trees, '[[', 3), lapply(tmp.trees, '[[', 4))
	names(MetaPIGA.trees)	<- c(sapply(tmp.trees, function(x) names(x)[1]), sapply(tmp.trees, function(x) names(x)[2]), sapply(tmp.trees, function(x) names(x)[3]), sapply(tmp.trees, function(x) names(x)[4]))	
	submitted.trees			<- c(submitted.trees, MetaPIGA.trees)	
	submitted.info			<- data.table(FILE=names(submitted.trees))
	
	submitted.info[, TEAM:=NA_character_]
	set(submitted.info, submitted.info[, which(grepl('RAXML',FILE))], 'TEAM', 'RAXML')
	set(submitted.info, submitted.info[, which(grepl('IQTree',FILE))], 'TEAM', 'IQTree')
	set(submitted.info, submitted.info[, which(grepl('MetaPIGA',FILE))], 'TEAM', 'MetaPIGA')
	set(submitted.info, submitted.info[, which(grepl('PhyML',FILE))], 'TEAM', 'PhyML')
	submitted.info[, SC:=NA_character_]
	tmp		<- submitted.info[, which(grepl('150701_Regional_TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Regional_TRAIN[0-9]',FILE))])
	tmp		<- submitted.info[, which(grepl('150701_Vill_SCENARIO-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_Vill_SCENARIO-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('TRAIN[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Regional_',regmatches(FILE, regexpr('TRAIN[0-9]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('scenario[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, paste('150701_Vill_',regmatches(FILE, regexpr('scenario[A-Z]',FILE)),sep='')])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_regional_train[0-9]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_regional_train[0-9]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('150701_vill_scenario-[A-Z]', FILE))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, regmatches(FILE, regexpr('150701_vill_scenario-[A-Z]',FILE))])
	tmp		<- submitted.info[, which(is.na(SC) & grepl('Vill_99_Apr15', FILE))]
	set(submitted.info, tmp, 'SC', 'Vill_99_Apr15')
	
	set(submitted.info, NULL, 'SC', submitted.info[, toupper(SC)])
	tmp		<- submitted.info[, which(grepl('150701_VILL_SCENARIO[A-Z]', SC))]
	set(submitted.info, tmp, 'SC', submitted.info[tmp, gsub('150701_VILL_SCENARIO','150701_VILL_SCENARIO-',SC)])
	
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/submitted_150911.rda'
	save(submitted.trees, submitted.info, file=outfile)
}