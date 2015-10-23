##--------------------------------------------------------------------------------------------------------
##	olli 23.10.15
##--------------------------------------------------------------------------------------------------------
treedist.pathdifference<- function(otree,stree, lambda=1)
{
	l				<- length(otree$tip.label)
	dt1				<- dist.nodes(otree)[1:l, 1:l]
	dt2				<- dist.nodes(stree)[1:l, 1:l]
	rownames(dt1)	<- colnames(dt1)	<- otree$tip.label
	rownames(dt2)	<- colnames(dt2)	<- stree$tip.label
	ct1				<- cophenetic.phylo(otree)[1:l, 1:l]
	ct2				<- cophenetic.phylo(stree)[1:l, 1:l]
	dt2				<- dt2[rownames(dt1),colnames(dt1)]
	ct2				<- ct2[rownames(ct1),colnames(ct1)]
	ind				<- lower.tri(dt1)
	pd				<- sum((dt1[ind] - dt2[ind])^2)
	cd				<- sum((ct1[ind] - ct2[ind])^2)
	ld				<- sum((lambda*dt1[ind]+(1-lambda)*ct1[ind] - lambda*dt2[ind]-(1-lambda)*ct2[ind])^2)
	c('path'=sqrt(pd), 'path.std'=sqrt(pd/choose(l,2)), 'pathl'=sqrt(ld))
}
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
	#	read sequences and determine %gappiness
	tmp		<- list.files(indir, pattern='fa$|fasta$', full.names=TRUE)
	tmp		<- data.table( FILE_SEQ_T= tmp, SC= toupper(gsub('_SIMULATED','',gsub('.fa','',basename(tmp)))))
	z		<- subset(tmp, SC=='VILL_99_APR15')
	set(z, NULL, 'SC', '150701_VILL_SCENARIO-C')
	tmp		<- rbind( tmp, z )	
	tfiles	<- merge(tfiles, tmp, by='SC', all=1)
	tmp		<- subset(tfiles, !is.na(FILE_SEQ_T))[, {
				z		<- read.dna(FILE_SEQ_T, format='fasta')	
				ans		<- sapply(seq_len(nrow(z)), function(i) base.freq(z[i,], all=1))
				ans		<- apply(ans[c('n','-','?'),], 2, sum)
				list(TAXA=rownames(z), GPS=ans)				
			}, by='SC']
	tinfo	<- merge(tinfo, tmp, by=c('SC','TAXA'), all.x=1)
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
	#IDX_T<-10
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
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
	#	compute path differences on complete trees
	#
	setkey(submitted.info, IDX)
	#tmp				<- subset(submitted.info, IDX==463)[1,]
	#IDX<- 1; IDX_T<- 1
	#IDX<- 822; IDX_T<- 11
	tmp				<- submitted.info[, {
				cat('\nAt IDX', IDX)
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
				#normalize with choose(n,2)		
				tmp			<- treedist.pathdifference(otree, stree, lambda=0)
				list(PD=tmp['path'], NPD=tmp['path.std'])
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
	outfile	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim/201507_TreeReconstruction/evaluation/submitted_151023.rda'
	save(strs, ttrs, tinfo, submitted.info, sclu.info, file=outfile)
}
treecomparison.ana.151019<- function()
{
	require(ggplot2)
	require(gamlss)
	
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
	set(sa, NULL, 'EXT', sa[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sa, NULL, 'ACUTE', sa[, factor(ACUTE, levels=c('low','high'),labels=c('10%','40%'))])
	set(sa, NULL, 'ART', sa[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sa		<- subset(sa, OTHER=='N')
	
	sc		<- copy(sclu.info)
	tmp		<- subset(tinfo, !is.na(IDCLU))[, list(CLU_N=CLU_N[1], MXGPS_CLU= max(GPS), MDGPS_CLU=median(GPS)), by=c('SC','IDCLU')]
	sc		<- merge(sc, tmp, by=c('SC','IDCLU'))	
	set(sc, NULL, 'MODEL', sc[, factor(MODEL, levels=c('V','R'),labels=c('Model: Village','Model: Regional'))])
	set(sc, sc[, which(SC=="VILL_99_APR15")],'SC',"150701_VILL_SCENARIO-C")	
	set(sc, NULL, 'SC', sc[, factor(SC,	levels=c("150701_REGIONAL_TRAIN1", "150701_REGIONAL_TRAIN2", "150701_REGIONAL_TRAIN3", "150701_REGIONAL_TRAIN4","150701_REGIONAL_TRAIN5","150701_VILL_SCENARIO-A","150701_VILL_SCENARIO-B","150701_VILL_SCENARIO-C","150701_VILL_SCENARIO-D","150701_VILL_SCENARIO-E"), 
							labels=c('sc 1','sc 2','sc 3','sc 4','sc 5','sc A','sc B','sc C','sc D','sc E'))])
	set(sc, NULL, 'GAPS', sc[, factor(GAPS, levels=c('none','low','high'),labels=c('Gaps: none','Gaps: low','Gaps: high'))])
	set(sc, NULL, 'BEST', sc[, factor(BEST, levels=c('Y','N'),labels=c('best tree','replicate tree'))])									
	set(sc, NULL, 'GENE', sc[, factor(GENE, levels=c('POL','GAG+POL+ENV'),labels=c('pol','gag+pol+env'))])	
	set(sc, NULL, 'TEAM', sc[, factor(TEAM, levels=sc[, sort(unique(TEAM))],labels=sc[, sort(unique(TEAM))])])
	set(sc, NULL, 'EXT', sc[, factor(EXT, levels=c('~0pc','5pc'),labels=c('~ 0%/year','5%/year'))])
	set(sc, NULL, 'ART', sc[, factor(ART, levels=c('none','fast'),labels=c('none','fast'))])
	sc		<- subset(sc, OTHER=='N')
	
	
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
	ggplot( subset(sa, TEAM!='MetaPIGA'), aes(y=NPD, x=SC, shape=TEAM, fill=GENE, colour=GENE, size=BEST) ) + 
			geom_jitter(position = position_jitter(height = .01, width=0.2)) +			
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			scale_fill_brewer(palette='Paired') +
			scale_colour_brewer(palette='Paired') +
			facet_wrap(MODEL~GAPS, scales='free') +	
			labs(x='\nsimulated data set', y='Path difference\n(standardized)\n', size='', shape='Method', fill='part of genome', colour='part of genome') +
			theme_bw() 
	ggsave(w=10, h=6, file=paste(edir,'/151023_PD_polvsall_by_gaps.pdf',sep=''))
	
	#	RF may be confounded by size of data set when evaluating the extent that regional is more difficult
	#	-->
	#	hard to extrapolate how standardized RF grows with size of data set,
	#	but regression extrapolation suggests there is an effect
	z		<- subset(sa, TEAM!='MetaPIGA' & !grepl('Reg',MODEL))
	mo		<- gamlss(NRF~TAXAN+GENE+GAPS, sigma.formula=~TAXAN+GENE+GAPS, family=BE(mu.link='cauchit'), data=z)
	tmp		<- subset(sa, TEAM!='MetaPIGA')
	tmp[, NRFP:=predict(mo, data=z, newdata=subset(tmp, select=c(TAXAN,GENE,GAPS)), what='mu',type='response')]
	ggplot( tmp, aes(x=TAXAN) ) + 
			geom_jitter(aes(y=NRF, shape=TEAM, colour=EXT, fill=EXT,  size=BEST), position = position_jitter(height = .01, width=20), alpha=0.7) +
			geom_line(aes(y=NRFP), colour='black', size=0.5) +
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,23,24)) +
			scale_fill_brewer(palette='Set1') + scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +
			facet_grid(GAPS~GENE) +
			labs(x='\nsize of simulated data set', y='Robinson-Fould\n(standardized)\n', size='', shape='Method', fill='trms/outside', colour='trms/outside') +
			theme_bw()
	ggsave(w=10, h=8, file=paste(edir,'/151020_RF_trmsoutside.pdf',sep=''))
	
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
	
	
	#	taxa excluded:	plot cluster RF as a function of cluster size
	#	-->
	#	excluding taxa did not lead to noticeably lower RFs
	tmp		<- sc[, list(NRFC=median(NRFC)), by=c('CLU_N','SC','GENE','TEAM')]
	ggplot( subset(sc, GENE=='gag+pol+env'), aes(y=NRFC, x=TAXA_NC, size=BEST, shape=TEAM, fill=TAXA_NC<CLU_N, colour=TAXA_NC<CLU_N)) +
			geom_jitter(position = position_jitter(height = .01, width=0.1), alpha=0.7) +
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,22,24), guide=FALSE) +
			scale_fill_brewer(palette='Set1') +
			scale_colour_brewer(palette='Set1') +
			scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.1)) +			
			#coord_trans(x='log10') +
			scale_x_log10(breaks=c(1,2,3,4,5,6,8,10,20,50,100,200,300), minor_breaks=NULL) +
			labs(x='\nsize of sampled transmission cluster', y='Robinson-Fould\n(standardized per transmission cluster)\n', size='', shape='Method', fill='taxa excluded\nprior to reconstruction', colour='taxa excluded\nprior to reconstruction') +
			facet_grid(SC~TEAM) +
			theme_bw() 
	ggsave(w=16, h=8, file=paste(edir,'/151020_RFCLU_iftaxaexcludedbeforetreereconstruction.pdf',sep=''))
	
		
	#	effect of acute in terms of RF? --> Yes
	ggplot( subset(sa, TEAM!='MetaPIGA' & TEAM!='PhyML' & grepl('Reg',MODEL) & !grepl('none',GAPS)), aes(y=NRF, x=ACUTE, shape=TEAM, fill=ACUTE, colour=ACUTE) ) + 
			geom_jitter(aes(size=BEST), position = position_jitter(height = .01, width=0.2), alpha=0.8) +
			geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
			scale_size_manual(values=c(3, 1)) +
			scale_shape_manual(values=c(21,24), guide=FALSE) +
			scale_fill_brewer(palette='Set1', guide=FALSE) +	scale_colour_brewer(palette='Set1', guide=FALSE) +
			facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
			labs(x='\ntransmissions from those in acute infection', y='Robinson-Fould\n(standardized)\n', size='') +
			theme_bw() 
	ggsave(w=10, h=8, file=paste(edir,'/151020_RF_impactAcute.pdf',sep=''))

	#	effect of ART roll out in terms of RF? --> No
	ggplot( subset(sa, TEAM!='MetaPIGA' & grepl('Vill',MODEL) & !grepl('none',GAPS)), aes(y=NRF, x=ART, shape=TEAM, fill=ART, colour=ART) ) + 
		geom_jitter(aes(size=BEST), position = position_jitter(height = .01, width=0.2), alpha=0.8) +
		geom_boxplot(outlier.shape=NA, colour='black', alpha=0.3) +
		scale_size_manual(values=c(3, 1)) +
		scale_shape_manual(values=c(21,23,24), guide=FALSE) +
		scale_fill_brewer(palette='Set2', guide=FALSE) +	scale_colour_brewer(palette='Set2', guide=FALSE) +
		facet_grid(GAPS~TEAM+GENE, scales='free_x') +	
		labs(x='\nART roll-out', y='Robinson-Fould\n(standardized)\n', size='') +
		theme_bw()
	ggsave(w=10, h=8, file=paste(edir,'/151020_RF_impactART.pdf',sep=''))

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