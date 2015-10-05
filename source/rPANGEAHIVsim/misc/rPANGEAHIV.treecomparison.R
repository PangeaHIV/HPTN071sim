##--------------------------------------------------------------------------------------------------------
##	olli 05.10.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.treecomparison.submissions<- function()	
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