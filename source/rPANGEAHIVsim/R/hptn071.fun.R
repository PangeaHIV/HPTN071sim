##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.unique<- function(seq.DNAbin.matrix)
{
	x<- as.character(seq.DNAbin.matrix)
	x<- apply(x, 1, function(z) paste(z,collapse=''))
	seq.DNAbin.matrix[!duplicated(x),]			
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.singleton2bifurcatingtree<- function(ph.s, dummy.label=NA)
{	
	if(!Nnode(ph.s))
	{
		stopifnot(!nrow(ph.s$edge))
		if(is.na(dummy.label))
			dummy.label		<- paste('DUMMY',ph.s$tip.label, sep='_')
		ph.s$edge			<- matrix(c(3,1,3,2), nrow=2, ncol=2, byrow=TRUE) 	
		ph.s$edge.length	<- c(0,0)
		ph.s$tip.label		<- c(ph.s$tip.label, dummy.label)
		ph.s$Nnode			<- 1
	}
	ph.s
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.collapse.singles<- function (tree) 
{
	elen 		<- tree$edge.length
	xmat 		<- tree$edge
	node.lab 	<- tree$node.label
	nnode 		<- tree$Nnode
	ntip 		<- length(tree$tip.label)
	root		<- 0
	singles 	<- NA
	while (length(singles) > 0) 
	{
		tx <- tabulate(xmat[, 1])
		singles <- which(tx == 1)
		if (length(singles) > 0) 
		{
			i 					<- singles[1]
			prev.node 			<- which(xmat[, 2] == i)
			next.node 			<- which(xmat[, 1] == i)
			xmat[prev.node, 2] 	<- xmat[next.node, 2]
			xmat 				<- xmat[xmat[, 1] != i, , drop=0]
			xmat[xmat > i] 		<- xmat[xmat > i] - 1L
			if(!length(prev.node))
				root			<- root + elen[next.node]
			if(length(prev.node))
				elen[prev.node] <- elen[prev.node] + elen[next.node]
			if (!is.null(node.lab)) 
				node.lab <- node.lab[-c(i - ntip)]
			nnode <- nnode - 1L
			elen <- elen[-next.node]
		}
	}
	tree$edge 			<- xmat
	tree$edge.length 	<- elen
	tree$node.label 	<- node.lab
	tree$Nnode 			<- nnode
	tree$root.edge		<- root
	tree
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
seq.read.newick<- function (file = "", text) 
{
	if (file != "") 
		text <- scan(file, sep = "\n", what = "character")	
	Nnode		<- length(gregexpr("\\(", text)[[1]])
	Ntip		<- 1+length(gregexpr(",", text)[[1]])	
	tree 		<- unlist(strsplit(text, NULL))
	tip.label 	<- vector(mode = "character")
	node.label	<- vector(mode = 'character')
	edge 		<- matrix(data = 0, Nnode + Ntip - 1, 2)
	edge.length <- rep(0, Nnode + Ntip - 1)
	ntip 		<- vector(mode = "numeric")
	currnode 	<- Ntip + 1
	nodecount 	<- currnode
	i 	<- 1
	j 	<- 1
	k 	<- 1
	while(tree[i] != ";") 
	{
		if(tree[i] == "(") 
		{
			edge[j, 1] <- currnode
			i <- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1					
				}
			}
			else if(tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
		else if(tree[i] == ")") 
		{
			i <- i + 1			
			if(is.na(match(tree[i], c(":", ")")))) 
			{
				l	<- gregexpr(":|;|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				node.label[currnode-Ntip]	<- substr(text, i, i+l-2)
				i	<- i+l-1				
			}
			if(tree[i] == ":") 
			{
				i 	<- i + 1
				l	<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				edge.length[match(currnode, edge[, 2])] <- as.numeric(substr(text, i, i+l-2))
				i	<- i+l-1	
			}
			ntip[match(currnode, edge[, 2])] 	<- sum(ntip[which(edge[, 1] == currnode)])
			currnode 							<- edge[match(currnode, edge[, 2]), 1]
		}
		else if(tree[i]==",") 
		{
			edge[j, 1] 	<- currnode
			i 			<- i + 1
			if(is.na(match(tree[i], c("(", ")", ",", ":", ";")))) 
			{
				l				<- gregexpr(",|:|\\)", substr(text, i, nchar(text)))[[1]][1]
				stopifnot(l>0)
				tip.label[k] 	<- substr(text, i, i+l-2)	
				i				<- i+l-1
				edge[j, 2] 		<- k
				k 				<- k + 1
				ntip[j] 		<- 1
				if (tree[i] == ":") 
				{
					i				<- i + 1
					l				<- gregexpr(",|\\)", substr(text, i, nchar(text)))[[1]][1]
					stopifnot(l>0)
					edge.length[j]	<- as.numeric(substr(text, i, i+l-2))
					i				<- i+l-1
				}
			}
			else if (tree[i] == "(") 
			{
				nodecount 	<- nodecount + 1
				currnode 	<- nodecount
				edge[j, 2] 	<- currnode
			}
			j <- j + 1
		}
	}
	tmp	<- which(edge[,1]==0)
	if(length(tmp))
	{
		edge		<- edge[-tmp,]
		edge.length	<- edge.length[-tmp]		
		tmp			<- sort( unique( as.numeric( edge ) ) )
		tmp			<- rbind(tmp, seq_along(tmp))		
		tmp			<- sapply( as.numeric( edge ), function(j)	tmp[2, match(j, tmp[1,])] )
		edge		<- matrix(tmp, ncol=2)
	}
	phy <- list(edge = edge, Nnode = as.integer(Nnode), tip.label = tip.label, Ndesc = ntip)
	if(sum(edge.length) > 1e-08) 
		phy$edge.length	<- edge.length
	if(length(node.label))
		phy$node.label	<- node.label
	class(phy) <- "phylo"
	return(phy)
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nodestats <- function(bstr) 
{
	#	remove anything before first '('
	bstr	<- regmatches(bstr, regexpr('\\(.*',bstr))
	# 	store meta info for inner nodes that is given in [], and not in :[] which is meta info for edges	
	tmp		<- unlist(regmatches(bstr,gregexpr('[^:]\\[[^]]+',bstr)))
	stopifnot(length(tmp)>0)
	tmp		<- sapply( tmp, function(x) substr(x, 4, nchar(x)) ) 
	#	for each inner node, extract stats
	tmp		<- strsplit(tmp, ',')
	tmp		<- lapply(seq_along(tmp), function(i)
			{
				z<- strsplit(tmp[[i]],'=')				
				data.table(NODE_PARSE_ID=i, STAT=sapply(z,'[',1), VALUE=sapply(z,'[',2))
			})
	node.stat	<- do.call('rbind', tmp)
	tmp			<- node.stat[, unique(STAT)]
	cat(paste('\nFound node statistics=',paste(tmp,collapse=' ')))
	tmp			<- node.stat[, list(has.all.stats= !length(setdiff(tmp, STAT))  ) , by='NODE_PARSE_ID']
	tmp			<- subset(tmp, !has.all.stats)[, NODE_PARSE_ID]
	cat(paste('\nSome statistics missing for nodes=',paste(tmp,collapse=' ')))
	node.stat 
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nodeidtree <- function(bstr, method.node.stat='any.node') 
{
	# strip all meta variables and ; at end
	bstr		<- gsub("\\[[^]]*\\]", "", bstr)
	bstr		<- gsub(';','',bstr)
	# for each node, add a dummy node label NODE_PARSE_IDxx	
	dummy.tree	<- unlist(strsplit(bstr, ":"))
	if(method.node.stat=='inner.node')
	{
		#	interior branch length: 	previous index ends in ). so tmp is the index of the dummy.tree chunks that gives the start of a branch length of an inner node
		tmp			<- which( c(FALSE, grepl(')$',dummy.tree)[-length(dummy.tree)]) )
		#	prepend NODE_PARSE_IDxx before the branch length of an inner node
		tmp			<- tmp-1			
	}
	if(method.node.stat=='any.node')
		tmp			<- seq_along(dummy.tree)
	dummy.tree	<- sapply(seq_along(dummy.tree), function(i)
			{
				z<- which(i==tmp)
				ifelse(length(z),	paste(dummy.tree[i],'NODE_PARSE_ID',z,sep=''),	dummy.tree[i] )
			}) 			
	dummy.tree	<- paste(dummy.tree, collapse=':',sep='')
	dummy.tree	<- regmatches(dummy.tree, regexpr('\\(.*',dummy.tree))
	dummy.tree	<- paste(dummy.tree, ';', sep='')	
	ph			<-  seq.read.newick(text=dummy.tree)
	ph
}
##--------------------------------------------------------------------------------------------------------
#	olli copied from hivclust
##--------------------------------------------------------------------------------------------------------
hivc.beast2out.read.nexus.and.stats<- function(file, tree.id=NA, method.node.stat='any.node') 
{	
	stopifnot(method.node.stat%in%c('any.node','inner.node'))
	
	X				<- scan(file = file, what = "", sep = "\n", quiet = TRUE)	
	#	read TRANSLATE chunk
	X.endblock		<- grep("END;|ENDBLOCK;|End;", X, ignore.case = TRUE)
	X.semico 		<- grep(";", X)
	X.i1 			<- grep("BEGIN TREES;|Begin trees;", X, ignore.case = TRUE)
	X.i2 			<- grep("TRANSLATE|Translate", X, ignore.case = TRUE)	
	tmp 			<- X.semico[X.semico > X.i2][1]
	tmp 			<- X[(X.i2 + 1):tmp]
	tmp				<- gsub('[,;]$','',gsub('^\\s+','',tmp))
	tmp				<- tmp[nzchar(tmp)]
	tmp				<- strsplit(tmp, ' ')
	df.translate	<- data.table(NEXUS_ID= sapply(tmp, '[[', 1), NEXUS_LABEL=sapply(tmp, '[[', 2) )
	set(df.translate, NULL, 'NEXUS_LABEL', df.translate[, gsub("\'","",NEXUS_LABEL)])
	cat(paste('\nFound taxa, n=', nrow(df.translate)))
	
	if(!is.na(tree.id))
	{
		#	read one newick tree with id 'tree.id'
		bstr		<- X[grep(paste(tree.id,"[[:space:]]+",sep=''), X)]
		node.stat	<- hivc.beast2out.read.nodestats(bstr)
		cat(paste('\nFound node statistics, n=', nrow(node.stat)))
		set(node.stat, NULL, 'tree.id', tree.id[i] )		
		btree		<- hivc.beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat) 
		#
		# link node.stats with tree nodes (tip + inner node)
		# NODE_ID is index of node in 'btree' phylo object
		#
		tmp			<- strsplit( btree$tip.label, 'NODE_PARSE_ID' )
		df.link		<- data.table(NODE_ID=seq_along(btree$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
		df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
		cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
		tmp			<- strsplit( btree$node.label, 'NODE_PARSE_ID' )
		tmp			<- data.table(NODE_ID=Ntip(btree)+seq_along(btree$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
		df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
		set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
		set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])
		set(df.link,NULL,'TREE_ID',tree.id)
		node.stat	<- merge( node.stat, subset(df.link, select=c(NODE_PARSE_ID, NODE_ID, TREE_ID)), by='NODE_PARSE_ID' )
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
		cat(paste('\nLinked node statistics to tree nodes, n=', nrow(node.stat)))
		#
		# set tip.labels and rm node.labels
		#
		setkey(df.link, NODE_ID)
		btree$tip.label		<- df.link[seq_len(Ntip(btree)),][,NEXUS_LABEL]
		btree$node.label	<- NULL		
	}
	if(is.na(tree.id))
	{
		#	read all newick trees in nexus file
		tmp			<- regexpr('^tree\\s\\S+',X)
		tree.id		<- sapply( regmatches(X,tmp), function(x) substr(x, 5, nchar(x)))
		tree.id		<- gsub('\\s','',tree.id)		
		cat(paste('\nFound tree id=', paste(tree.id, collapse=' ')))
		X			<- X[ which(tmp>0) ]
		cat(paste('\nFound trees, n=',length(tree.id)))
		node.stat	<- lapply(seq_along(tree.id), function(i)
				{
					
					bstr	<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
					cat(paste('\nGet node statistics for tree id=',tree.id[i]))
					tmp		<- hivc.beast2out.read.nodestats(bstr)
					set(tmp, NULL, 'TREE_ID', tree.id[i] )
					tmp
				})
		if(length(node.stat)>1)
			node.stat	<- do.call('rbind',node.stat)
		if(length(node.stat)==1)
			node.stat	<- node.stat[[1]]
		suppressWarnings({ node.stat[, NODE_ID:=NA_integer_] })				
		setkey(node.stat, TREE_ID, NODE_PARSE_ID)
		btree		<- vector('list',length(tree.id))
		for(i in seq_along(tree.id))
		{
			bstr		<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
			cat(paste('\nRead tree for tree id=',tree.id[i]))
			btree.i		<- hivc.beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat)
			#
			# link node.stats with tree nodes (tip + inner node)
			# NODE_ID is index of node in 'btree.i' phylo object
			#
			tmp			<- strsplit( btree.i$tip.label, 'NODE_PARSE_ID' )
			df.link		<- data.table(NODE_ID=seq_along(btree.i$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
			df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
			cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
			tmp			<- strsplit( btree.i$node.label, 'NODE_PARSE_ID' )
			tmp			<- data.table(NODE_ID=Ntip(btree.i)+seq_along(btree.i$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
			df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
			set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
			set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])		
			for(j in seq_len(nrow(df.link)))
				set(node.stat, node.stat[, which(TREE_ID==tree.id[i] & NODE_PARSE_ID==df.link[j,NODE_PARSE_ID])], 'NODE_ID', df.link[j,NODE_ID])		
			tmp			<- node.stat[, length(which(!is.na(NODE_ID)))]
			cat(paste('\nTotal linked node statistics to tree nodes, n=', tmp  ))
			#
			# set tip.labels and rm node.labels
			#
			setkey(df.link, NODE_ID)
			btree.i$tip.label	<- df.link[seq_len(Ntip(btree.i)),][,NEXUS_LABEL]
			btree.i$node.label	<- NULL	
			btree[[i]]			<- btree.i
		}
		if(length(btree)>=2)
		{
			names(btree)	<- tree.id
			class(btree)	<- "multiPhylo"			
		}
		if(length(btree)<2)
			btree	<- btree[[1]]
		tmp				<- node.stat[, length(which(is.na(NODE_ID)))]
		cat(paste('\nTotal unlinked node statistics [should be zero], n=', tmp  ))
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
	}
	list(tree=btree, node.stat=node.stat)	 
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 20-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create data.table of GTR parameters
#' @description Returns a data.table of GTR parameters. 
#' @return data.table
#' @export
PANGEA.GTR.params<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='ENV' & FILE=='pool3'))
	#	exclude cols
	log.df[, ucld.mean:=NULL]
	log.df[, ucld.stdev:=NULL]
	log.df[, coefficientOfVariation:=NULL]
	log.df[, treeModel.rootHeight:=NULL]
	#	check that relative rates have mean 1
	stopifnot( log.df[, list(CHECK=sum(mu)), by=c('FILE','GENE','state')][, !any(abs(CHECK-3)>2*1e-12)] )
	#	get mean rate. need to be a bit careful here: rate is per site, so length of gene does not matter
	#	but we may have different numbers of samples for each gene
	tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
	#	ignore variation in meanRate by gene, keep variation across codon_pos for each state
	log.df	<- merge(log.df, log.df[, list(mu.gene= mean(meanRate)/tmp), by='GENE'], by='GENE')	
	set(log.df, NULL, 'mu', log.df[, mu*mu.gene])
	set(log.df, NULL, 'meanRate', tmp)
	set(log.df, NULL, 'mu.gene', NULL)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 14-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v3<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='ENV' & FILE=='pool3'))
	#	exclude cols
	log.df[, ucld.mean:=NULL]
	log.df[, ucld.stdev:=NULL]
	log.df[, coefficientOfVariation:=NULL]
	log.df[, treeModel.rootHeight:=NULL]
	#	get mean rate. need to be a bit careful here: rate is per site, so length of gene does not matter
	#	but we may have different numbers of samples for each gene
	tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
	#	set mean meanRate and put all variation into the mu's
	set(log.df, NULL, 'mu', log.df[, mu * meanRate / tmp])
	set(log.df, NULL, 'meanRate', tmp)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v2<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='GAG' & FILE=='pool1'))
	log.df		<- subset(log.df, !(GENE=='POL' & FILE=='pool2'))
	#	all ENV are a bit odd ... 
	#	exclude cols
	log.df[, ucld.mean:=NULL]
	log.df[, ucld.stdev:=NULL]
	log.df[, coefficientOfVariation:=NULL]
	log.df[, treeModel.rootHeight:=NULL]
	#	set mean meanRate and put all variation into the mu's
	tmp		<- log.df[, mean(meanRate)]
	set(log.df, NULL, 'mu', log.df[, mu * meanRate / tmp])
	set(log.df, NULL, 'meanRate', tmp)
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	return distribution of GTR parameters	
#	olli originally written 09-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.GTR.params.v1<- function()
{	
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140811_n390_BEASTlog.R')		
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	log.df[, state:=NULL]
	log.df[, ucldmean:=NULL]
	log.df[, ucldstdev:=NULL]
	log.df[, treeLikelihood:=NULL]
	log.df[, FILE:=NULL]
	log.df
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase.v2<- function(df.ind, df.trm, index.starttime.mode='normal')
{
	stopifnot(grepl('normal|fix|shift', index.starttime.mode))
	#	add transmission time for index case -- this is 40 years back in time so we can sample a starting sequence 
	#	and then generate a long branch to the transmission chain in the population. No hack. :-)
	tmp			<- subset( df.trm, IDTR<0, select=IDTR )	
	if( index.starttime.mode == 'normal' )
	{
		cat(paste('\nUsing index.starttime.mode rnorm(n, 1955, 7)'))
		tmp2		<- rnorm(2*nrow(tmp), 1955, 7)
		tmp2		<- tmp2[ tmp2>1946 & tmp2<1980]
		stopifnot( nrow(tmp)<=length(tmp2) )
		length(tmp2)<- nrow(tmp)			
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	if( grepl('fix',index.starttime.mode) ) 
	{
		tmp2		<- as.numeric(substring(index.starttime.mode, 4))
		cat(paste('\nUsing index.starttime.mode rep=', tmp2))
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', rep(tmp2, nrow(tmp)) )
	}	
	if( grepl('shift',index.starttime.mode))
	{
		tmp2		<- as.numeric(substring(index.starttime.mode, 6))
		cat(paste('\nUsing index.starttime.mode rep=', tmp2))
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', rep(tmp2, nrow(tmp)) )		
	}
	if( index.starttime.mode == 'fix45')
	{
		cat(paste('\nUsing index.starttime.mode runif( n, 1945, 1945.5 )'))
		tmp2		<- runif( nrow(tmp), 1945, 1945.5 )
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	#
	df.trm		<- merge( df.trm, tmp, by='IDTR', all.x=TRUE )
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new) & IDTR_TIME_INFECTED<IDTR_TIME_INFECTED.new)]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED.new', df.trm[tmp2, IDTR_TIME_INFECTED])
	#	
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED', df.trm[tmp2, IDTR_TIME_INFECTED.new])
	df.trm[, IDTR_TIME_INFECTED.new:=NULL]
	#	
	stopifnot( nrow(subset(df.trm, TIME_TR<IDTR_TIME_INFECTED))==0 )
	stopifnot( nrow(subset(df.trm, is.na(TIME_TR)))==0 )
	stopifnot( nrow(subset(df.trm, is.na(IDTR_TIME_INFECTED)))==0 )
	#
	tmp			<- subset( df.trm, IDTR<0, select=c(IDTR, IDTR_TIME_INFECTED) )
	setnames(tmp, c('IDTR','IDTR_TIME_INFECTED'), c('IDPOP','IDTR_TIME_INFECTED.new'))
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	tmp2		<- df.ind[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.ind, tmp2, 'TIME_TR', df.ind[tmp2, IDTR_TIME_INFECTED.new])
	df.ind[, IDTR_TIME_INFECTED.new:=NULL]
	list(df.ind=df.ind, df.trm=df.trm)
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 13-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateStartingTimeOfIndexCase<- function(df.ind, df.trm, index.starttime.mode='normal')
{
	stopifnot(index.starttime.mode%in%c('normal','fix','fix45','shift'))
	#	add transmission time for index case -- this is 40 years back in time so we can sample a starting sequence 
	#	and then generate a long branch to the transmission chain in the population. No hack. :-)
	tmp			<- subset( df.trm, IDTR<0, select=IDTR )	
	if( index.starttime.mode == 'normal' )
	{
		cat(paste('\nUsing index.starttime.mode rnorm(n, 1955, 7)'))
		tmp2		<- rnorm(2*nrow(tmp), 1955, 7)
		tmp2		<- tmp2[ tmp2>1946 & tmp2<1980]
		stopifnot( nrow(tmp)<=length(tmp2) )
		length(tmp2)<- nrow(tmp)			
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	if( index.starttime.mode == 'fix' || index.starttime.mode == 'shift')
	{
		cat(paste('\nUsing index.starttime.mode runif( n, 1954.75, 1955.25 )'))
		tmp2		<- runif( nrow(tmp), 1954.75, 1955.25 )
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	if( index.starttime.mode == 'fix45')
	{
		cat(paste('\nUsing index.starttime.mode runif( n, 1945, 1945.5 )'))
		tmp2		<- runif( nrow(tmp), 1945, 1945.5 )
		set(tmp, NULL, 'IDTR_TIME_INFECTED.new', tmp2 )
	}
	#
	df.trm		<- merge( df.trm, tmp, by='IDTR', all.x=TRUE )
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new) & IDTR_TIME_INFECTED<IDTR_TIME_INFECTED.new)]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED.new', df.trm[tmp2, IDTR_TIME_INFECTED])
	#	
	tmp2		<- df.trm[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.trm, tmp2, 'IDTR_TIME_INFECTED', df.trm[tmp2, IDTR_TIME_INFECTED.new])
	df.trm[, IDTR_TIME_INFECTED.new:=NULL]
	#	
	stopifnot( nrow(subset(df.trm, TIME_TR<IDTR_TIME_INFECTED))==0 )
	stopifnot( nrow(subset(df.trm, is.na(TIME_TR)))==0 )
	stopifnot( nrow(subset(df.trm, is.na(IDTR_TIME_INFECTED)))==0 )
	#
	tmp			<- subset( df.trm, IDTR<0, select=c(IDTR, IDTR_TIME_INFECTED) )
	setnames(tmp, c('IDTR','IDTR_TIME_INFECTED'), c('IDPOP','IDTR_TIME_INFECTED.new'))
	df.ind		<- merge(df.ind, tmp, by='IDPOP', all.x=TRUE)
	tmp2		<- df.ind[, which(!is.na(IDTR_TIME_INFECTED.new))]
	set(df.ind, tmp2, 'TIME_TR', df.ind[tmp2, IDTR_TIME_INFECTED.new])
	df.ind[, IDTR_TIME_INFECTED.new:=NULL]
	list(df.ind=df.ind, df.trm=df.trm)
}
##--------------------------------------------------------------------------------------------------------
#	simulates seroprevalence survey of a proportion of a population	
#	olli originally written 29-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.SeroPrevalenceSurvey<- function(df.inds, epi.adult=13, s.INTERVENTION.start=2015, sp.prop.of.sexactive=0.05, sp.times=c(15, 10, 5, 0), verbose=1)
{		
	df.sp	<- lapply( s.INTERVENTION.start-sp.times, function(yr)
			{
				yr			<- yr+.5
				sxon.all	<- which( (df.inds[['DOB']]+epi.adult)<=yr  &  df.inds[['DOD']]>yr )
				sxon.sp		<- sample(sxon.all, round(length(sxon.all)*sp.prop.of.sexactive) )
				if(verbose)
				{
					cat(paste('\nSero prevalence survey in year', yr))
					cat(paste('\nTotal number of individuals=', length(sxon.all)))
					cat(paste('\nProp to include in survey=', sp.prop.of.sexactive))
					cat(paste('\nTotal number of individuals in survey=', length(sxon.sp)))
				}
				df.sp		<- df.inds[sxon.sp, ]				
				df.sp[, YR:=yr]
				df.sp[, AGE:= df.sp[, cut(YR-DOB, breaks=c(epi.adult-1, 16, 20, 25, 30, 35, 40, 50, 60, Inf))]]
				df.sp
			})
	df.sp	<- do.call('rbind', df.sp)	
	df.sp	<- df.sp[, list( ALIVE_N=length(IDPOP), ALIVE_AND_DIAG_N=length(which(DIAG_T<=YR)), ALIVE_AND_ART_N=length(which(ART1_T<=YR)), ALIVE_AND_SEQ_N=length(which(TIME_SEQ<=YR)) ), by=c('YR', 'GENDER', 'AGE')]
	setkey(df.sp, YR, GENDER, AGE)
	df.sp
}
##--------------------------------------------------------------------------------------------------------
#	simulate imports during the epidemic	
#	olli originally written 11-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.ImportSimulator.SimulateIndexCase<- function(df.ind, df.trm, epi.import)
{	
	#	model imports by re-assigning a fraction of infecteds as 'index cases'
	#	assume there are no imports at the moment, and that we know the time of infection of the imports
	tmp		<- df.trm[, which(IDTR>0)]	  
	cat(paste('\nFound transmissions within population, n=', length(tmp)))
	tmp2	<- round(nrow(df.trm)*epi.import) - ( nrow(df.trm)-length(tmp)-1 )
	tmp2	<- max(tmp2, 0)
	cat(paste('\nRe-setting infecteds as index cases after imports, n=', tmp2))
	stopifnot(length(tmp)>1)
	stopifnot(length(tmp2)>=0)
	tmp2	<- as.integer(sample( tmp, tmp2, replace=FALSE ))	
	#	update df.trm
	setkey(df.trm, TIME_TR)
	set(df.trm, tmp2, 'IDTR', df.trm[, min(IDTR)]-rev(seq_along(tmp2)) )
	if('TR_ACUTE'%in%colnames(df.trm))
		set(df.trm, df.trm[, which(IDTR<0)], 'TR_ACUTE', NA_character_)
	#	update df.ind
	tmp		<- subset(df.trm, select=c(IDTR, IDTR_TIME_INFECTED))
	setnames(tmp, c('IDTR', 'IDTR_TIME_INFECTED'), c('IDPOP', 'TIME_TR'))
	tmp2	<- subset(df.trm, select=c(IDREC, TIME_TR))
	setnames(tmp2, 'IDREC', 'IDPOP')
	tmp		<- rbind( tmp, tmp2 )
	setkey(tmp, IDPOP)
	df.ind	<- merge( df.ind, unique(tmp), by=c('IDPOP','TIME_TR'), all.x=TRUE, all.y=TRUE )
	#	
	cat(paste('\nProportion of imported transmissions, p=', (nrow(subset(df.trm, IDTR<0))-1)/nrow(df.trm) ))
	stopifnot( length(setdiff(df.trm[, IDTR], df.ind[, IDPOP]))==0 )
	stopifnot( length(setdiff(df.trm[, IDREC], df.ind[, IDPOP]))==0 )
	list(df.ind=df.ind, df.trm=df.trm)
}
##--------------------------------------------------------------------------------------------------------
#	simulate guide to sequence sampling times. if not NA, then in every year, individuals in +-6mo to the guide are sampled	
#	olli originally written 26-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2<- function(df.ind, seqtime.mode)
{
	stopifnot(grepl('DUnif|Exp|AtDiag|AtART|AtTrm',seqtime.mode))
	cat(paste('\nUsing seqtime.mode=', seqtime.mode ))
	if(grepl('Exp',seqtime.mode))
	{
		tmp		<- as.numeric(substring(seqtime.mode,4))
		stopifnot(is.finite(tmp))
		cat(paste('\nPANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2: avg time to sequencing since diagnosis=', tmp))		
		df.ind[, T1_SEQ:= rexp(nrow(df.ind), rate=1/tmp) + df.ind[,DIAG_T]]
		#	need T1_SEQ even when no diagnosis for archival samples
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )		
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
	}
	if(grepl('DUnif',seqtime.mode))
	{
		tmp		<- as.numeric(substring(seqtime.mode,6))
		stopifnot(is.finite(tmp))
		cat(paste('\nPANGEA.Seqsampler.SimulateGuideToSamplingTimes.v2: max time to sequencing since diagnosis=', tmp))		
		df.ind[, T1_SEQ:= runif(nrow(df.ind), min=0, max=tmp) + df.ind[,DIAG_T]]
		#	need T1_SEQ even when no diagnosis for archival samples
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )		
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
	}
	if(seqtime.mode=='AtTrm')
	{
		df.ind[, T1_SEQ:= df.ind[, TIME_TR+0.1]]	
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD<T1_SEQ)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, TIME_TR+(TIME_TR+DOD)/2] )
	}	
	if(seqtime.mode=='AtDiag')
	{
		df.ind[, T1_SEQ:= df.ind[, DIAG_T]]	
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
	}
	if(seqtime.mode=='AtART')
	{
		df.ind[, T1_SEQ:= df.ind[, ART1_T]]		
		tmp		<- df.ind[, which(is.na(T1_SEQ) & DOD-TIME_TR>0.5)]
		set( df.ind, tmp, 'T1_SEQ', df.ind[tmp, runif(length(tmp), TIME_TR+0.5, DOD)] )
		tmp		<- df.ind[, which( is.na(DIAG_T) & T1_SEQ>min(DIAG_T, na.rm=1) )]
		set( df.ind, tmp, 'T1_SEQ', NA_real_)
	}
	tmp	<- df.ind[, which(T1_SEQ>ART1_T)]
	set(df.ind, tmp, 'T1_SEQ', df.ind[tmp,ART1_T])
	df.ind
}
##--------------------------------------------------------------------------------------------------------
#	simulate guide to sequence sampling times. if not NA, then in every year, individuals in +-6mo to the guide are sampled	
#	olli originally written 24-10-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.SimulateGuideToSamplingTimes<- function(df.ind, seqtime.mode)
{
	stopifnot(seqtime.mode%in%c('Gamma3','Gamma9','Unif12'))
	cat(paste('\nUsing seqtime.mode=', seqtime.mode ))
	if(seqtime.mode=='Gamma3')
	{
		library(distr) 
		tmp		<- Gammad(shape=3, scale=2)
		tmp		<- Truncate(tmp, lower=0, upper=8)
		df.ind[, T1_SEQ:= df.ind[, r(tmp)(nrow(df.ind)) + TIME_TR]]			
	}
	if(seqtime.mode=='Gamma9')
	{
		df.ind[, T1_SEQ:= df.ind[, rgamma(nrow(df.ind),shape=9,scale=0.25 ) + TIME_TR]]	
	}
	if(seqtime.mode=='Unif12')
	{
		df.ind[, T1_SEQ:= runif(nrow(df.ind), 0, 12) + df.ind[,TIME_TR] ]
	}
	df.ind
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to untreated
#
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.untreated<- function(df.ind, df.epi, pipeline.args)
{	
	s.total					<- round( df.epi[nrow(df.epi), PREV] * pipeline.args['s.PREV.max',][, as.numeric(v)] )	
	cat(paste('\nSample proportional to untreated population'))
	cat(paste('\nSampling sequences, target is n=', s.total))
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )	
	set(df.sample, NULL, 's.nTOTAL', as.numeric(rmultinom(1, s.total, df.sample[, PREV-TREATED]/ df.epi[, sum(PREV-TREATED)])) )	
	cat(paste('\nSampling sequences, scheduled number is n=', df.sample[, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]]))
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))
		stopifnot(length(tmp)>0)
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ] )				
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to diagnoses before and after interventions
#	s% of those newly diagnosed per year until 2015
#	2*s% of those newly diagnosed per year after 2015
#	from before then: 50 (uniform)	
#
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.diagnosis<- function(df.ind, df.epi, pipeline.args)
{
	s.total					<- round( df.epi[nrow(df.epi), PREV] * pipeline.args['s.PREV.max',][, as.numeric(v)] )
	s.archival.yr			<- subset(df.epi, DIAG==0)[, tail(YR,1)]
	s.diagb4intervention.n	<- subset( df.epi, YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )[, tail(DIAG,1)]
	s.diagb4intervention	<- (s.total-pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]) / (s.diagb4intervention.n+pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))
	cat(paste('\nSample proportional to diagnoses up to intervention, and proportional to diagnoses after intervention start'))
	cat(paste('\nSampling sequences, target is n=', s.total))
	cat(paste('\nNo diagnoses up to yr, sampling archival sequences till then. yr=', s.archival.yr))
	cat(paste('\nSampling archival sequences before any diagnoses, n=', pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]))	
	cat(paste('\nSampling sequences before intervention start from diagnosed, p=', s.diagb4intervention))
	cat(paste('\nSampling sequences after intervention start from diagnosed, p=', 2*s.diagb4intervention))
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	tmp			<- df.sample[, which(YR<=s.archival.yr)]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], rep(1/length(tmp), length(tmp)))))
	tmp			<- df.sample[, which( YR>s.archival.yr & YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )]	
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( s.diagb4intervention.n*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	tmp			<- df.sample[, which(YR>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)]) ]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( (df.sample[nrow(df.sample), DIAG]-s.diagb4intervention.n)*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	cat(paste('\nSampling sequences, scheduled number is n=', df.sample[, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]]))
	stopifnot( abs(pipeline.args['s.PREV.max',][, as.numeric(v)]-df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]])<=pipeline.args['s.PREV.max',][, as.numeric(v)]*0.1 )
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))
		stopifnot(length(tmp)>0)
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ] )				
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	sample proportional to diagnoses before and after interventions
#	s% of those newly diagnosed per year until 2015
#	from before any diagnosed: 50 (uniform)
#	after 2015, a fixed amount per year
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.sample.prop.to.diagnosis.b4intervention<- function(df.ind, df.epi, pipeline.args)
{
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	s.total					<- round( df.epi[nrow(df.epi), PREV] * pipeline.args['s.PREV.max',][, as.numeric(v)] )
	s.archival.yr			<- subset(df.epi, DIAG==0)[, tail(YR,1)]
	s.diagb4intervention.n	<- subset( df.epi, YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )[, tail(DIAG,1)]
	s.diagb4intervention	<- (s.total-pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]) / (s.diagb4intervention.n+pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))
	cat(paste('\nSample proportional to diagnoses up to intervention, and then a fixed amount per year'))
	cat(paste('\nSampling sequences, target is n=', s.total))
	cat(paste('\nNo diagnoses up to yr, sampling archival sequences till then. yr=', s.archival.yr))
	cat(paste('\nSampling archival sequences before any diagnoses, n=', pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)]))	
	cat(paste('\nSampling sequences before intervention start from diagnosed, p=', s.diagb4intervention))
	cat(paste('\nSampling sequences after intervention start from diagnosed, n=', round(s.diagb4intervention*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*(df.epi[nrow(df.epi), DIAG]-s.diagb4intervention.n))))
	#	setup number of sequences to be sampled for each year
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	set(df.sample, NULL, 's.nTOTAL', 0)	
	tmp			<- df.sample[, which(YR<=s.archival.yr)]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, pipeline.args['s.ARCHIVAL.n',][, as.numeric(v)], rep(1/length(tmp), length(tmp)))))	
	cat(paste('\nSampling sequences from archival, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which( YR>s.archival.yr & YR<pipeline.args['s.INTERVENTION.start',][, as.numeric(v)] )]	
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( s.diagb4intervention.n*s.diagb4intervention ), df.sample[tmp, NEW_DIAG]/df.sample[tmp, sum(NEW_DIAG)])) )
	cat(paste('\nSampling sequences before intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	tmp			<- df.sample[, which(YR>=pipeline.args['s.INTERVENTION.start',][, as.numeric(v)]) ]
	set(df.sample, tmp, 's.nTOTAL', as.numeric(rmultinom(1, round( (df.sample[nrow(df.sample), DIAG]-s.diagb4intervention.n)*pipeline.args['s.INTERVENTION.mul'][, as.numeric(v)]*s.diagb4intervention ), rep(1/length(tmp),length(tmp))) ) )
	cat(paste('\nSampling sequences after intervention, scheduled number is n=', df.sample[tmp, sum(s.nTOTAL)]))
	stopifnot(df.sample[, all(s.nTOTAL>=0)])
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]]))	
	stopifnot( abs(pipeline.args['s.PREV.max',][, as.numeric(v)]-df.sample[, sum( s.nTOTAL )] / df.sample[, rev(PREV)[1]])<=pipeline.args['s.PREV.max',][, as.numeric(v)]*0.1 )
	#
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample non-incident cases by year
	df.inds		<- copy(df.ind)
	df.inds[, TIME_SEQ:= NA_real_]
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd samples in year',yr,', required=', subset( df.sample, YR==yr )[, s.nTOTAL]))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(TIME_TR)>=0)]
		cat(paste('\navailable non-sampled HIV+ individuals in year=', length(tmp)))
		tmp		<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==0)]
		cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date in year=', length(tmp)))		
		stopifnot(length(tmp)>0)
		k		<- 1
		while(length(tmp)<subset( df.sample, YR==yr )[, s.nTOTAL])
		{			
			cat(paste('\nCannot find samples, fall back to ',k,'th previous year; n=', subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp) ))
			tmp2	<- df.inds[, which(is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-floor(T1_SEQ)==k & (is.na(ART1_T) | ART1_T-yr<=1)) ]
			cat(paste('\navailable non-sampled HIV+ individuals with suggested sampling date one year before=', length(tmp2)))
			tmp2	<- sample(tmp2,  min(subset( df.sample, YR==yr )[, s.nTOTAL]-length(tmp), length(tmp2)) )
			tmp		<- c(tmp, tmp2)
			k		<- k+1
			stopifnot(k<=10)
		}	
		stopifnot(length(tmp)>=subset( df.sample, YR==yr )[, s.nTOTAL])
		tmp		<- sample(tmp, subset( df.sample, YR==yr )[, s.nTOTAL])
		set( df.inds, tmp, 'TIME_SEQ', df.inds[tmp, T1_SEQ+yr-floor(T1_SEQ)] )			
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))
	list(df.inds=df.inds, df.sample=df.sample)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 27-01-2015
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v4<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	stopifnot(grepl('Prop2DiagB4I|Prop2Diag|Prop2Untreated',pipeline.args['s.MODEL',][, v]))
	# 	compute prevalence and incidence by year	
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
				sexactive	<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				diag		<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['DIAG_T']])<=YR )
				diag.new	<- which( floor(df.ind[['DIAG_T']])==YR )
				treated		<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['ART1_T']])<=YR & (is.na(df.ind[['VLS1_TE']]) | floor(df.ind[['VLS1_TE']])>YR) )
				infdead		<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(sexactive), PREV=length(infected), PREVDIED=length(infdead), DIAG=length(diag), NEW_DIAG=length(diag.new), TREATED=length(treated))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )		
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	set(df.epi, NULL, 'ARTcov', df.epi[, TREATED/PREV])
	set(df.epi, NULL, 'UNDIAGp', df.epi[, (PREV-DIAG)/PREV])
	set(df.epi, NULL, 'GROWTHr', c(NA_real_, df.epi[, diff(log(PREV))]))
	
	if(pipeline.args['s.MODEL',][, v]=='Prop2DiagB4I')
		tmp			<- PANGEA.Seqsampler.sample.prop.to.diagnosis.b4intervention(df.ind, df.epi, pipeline.args)
	if(pipeline.args['s.MODEL',][, v]=='Prop2Diag')
		tmp			<- PANGEA.Seqsampler.sample.prop.to.diagnosis(df.ind, df.epi, pipeline.args)
	if(pipeline.args['s.MODEL',][, v]=='Prop2Untreated')
		tmp			<- PANGEA.Seqsampler.sample.prop.to.untreated(df.ind, df.epi, pipeline.args)
	df.inds		<- copy(tmp$df.inds)
	df.sample	<- copy(tmp$df.sample)			
	#	set sampling in df.trm
	tmp		<- subset( df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_SEQ) )
	setnames(tmp, c('IDPOP','TIME_SEQ'), c('IDREC','SAMPLED_REC'))
	df.trms	<- merge(df.trm, tmp, by='IDREC', all.x=TRUE)
	setnames(tmp, c('IDREC','SAMPLED_REC'), c('IDTR','SAMPLED_TR'))
	df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)	
	#
	#	seroprevalence survey
	#
	df.sp	<- PANGEA.SeroPrevalenceSurvey(df.inds, epi.adult=epi.adult, s.INTERVENTION.start=pipeline.args['s.INTERVENTION.start', ][, as.numeric(v)], sp.prop.of.sexactive=pipeline.args['sp.prop.of.sexactive', ][, as.numeric(v)], sp.times=c(15, 10, 5, 0) )
	file	<- gsub('IND.csv','CROSS_SECTIONAL_SURVEY.csv', outfile.ind)
	cat(paste('\nwrite to file', file))
	write.csv(df.sp, file=file)
	#
	#	TRANSMISSION NETWORKS
	#
	require(igraph)
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#	add IDCLU to df.inds
	tmp			<- subset( df.trms, select=c(IDREC, IDTR, IDCLU) )
	tmp			<- subset( melt(tmp, id.var='IDCLU', value.name='IDPOP'), select=c(IDPOP, IDCLU))
	setkey(tmp, IDPOP, IDCLU)
	tmp			<- unique(tmp)
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		tmp		<- names(df.sample)[ which(as.character(sapply(df.sample, class))!='numeric') ]
		for(x in tmp)
			set(df.sample, NULL, x, as.numeric(df.sample[[x]]))
		tmp		<- melt(df.sample, id.vars=c('YR'))
		set(tmp, NULL, 'variable', tmp[, factor(variable, levels=tmp[, unique(variable)], labels=tmp[, unique(variable)])])
		
		ggplot(tmp, aes(x=YR, y=value, group=variable)) + geom_step(with.guide=FALSE) + 
				facet_grid(variable~., scales='free_y')  + 
				theme_bw() + scale_x_continuous(name='year', breaks=seq(1980,2020,2)) +
				theme(strip.text=element_text(size=7))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=16)	
		#	plot CD4 counts at time of diagnosis
		tmp	<- subset(df.inds, !is.na(DIAG_T), select=c(DIAG_T, DIAG_CD4))
		tmp[, YR:= floor(DIAG_T)]
		tmp	<- merge(tmp, tmp[, list(DIAG_CD4_ME=mean(DIAG_CD4)), by='YR'], by='YR')
		setkey(tmp, YR)
		ggplot(tmp) + geom_point(aes(x=DIAG_T, y=DIAG_CD4)) + geom_step(data=unique(tmp), aes(x=YR, y=DIAG_CD4_ME), col='red', size=1.5) + scale_y_continuous(breaks=seq(100,1000,100))		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)					
		ggplot(subset(df.inds, !is.na(DIAG_T)), aes(x=floor(DIAG_T), fill=cut(DIAG_CD4, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent of diagnosed', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP_AMONGDIAG.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)		
		ggplot(subset(df.inds, !is.na(DIAG_T) & !is.na(TIME_SEQ)), aes(x=floor(DIAG_T), fill=cut(DIAG_CD4, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent of sequenced', x='year diagnosed') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP_AMONGSEQ.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)		
		#	plot distribution between transmission time and sequencing time		
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=TIME_SEQ-TIME_TR)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution between transmission time and diagnosis
		ggplot(subset(df.inds, !is.na(DIAG_T)), aes(x=DIAG_T-TIME_TR)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to diagnosis\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Diag.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution of #recipients per infector
		tmp	<- df.trms[, list(N= length(IDREC)), by='IDTR']
		ggplot(tmp, aes(x=N)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='recipients per source case\n(number)', breaks=seq(1,100,1)+0.5, label=seq(1,100,1)) +
				scale_y_continuous(breaks=seq(0,1,0.05))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_RecSource.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot time to death for infected
		tmp	<- subset(df.inds, !is.na(TIME_TR) & IDPOP>0 & DOD<max(DOD, na.rm=1), select=c(TIME_TR, DOD))
		ggplot(tmp, aes(x=DOD-TIME_TR)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='time to death for HIV+\n(years)', breaks=seq(1,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_T2DeathForInf.pdf',sep='')
		ggsave(file=file, w=8, h=8)		
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu & IDPOP>=0, select=c(IDPOP, GENDER, TIME_SEQ))
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
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
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
	if('RISK'%in%colnames(df.inds))
		df.inds[, RISK:=NULL]
	if('INCIDENT_SEQ'%in%colnames(df.inds))
		df.inds[, INCIDENT_SEQ:=NULL]
	if('TIME_SEQYR'%in%colnames(df.inds))
		df.inds[, TIME_SEQYR:=NULL]	
	if('TR_ACUTE'%in%colnames(df.trms))
		df.trms[, TR_ACUTE:=NULL]
	if('YR'%in%colnames(df.trms))
		df.trms[, YR:=NULL]	
	#	add columns that the virus tree simulator needs
	if(!'IDTR_TIME_INFECTED'%in%colnames(df.trms))
	{
		tmp		<- subset( df.inds, !is.na(TIME_TR), c(IDPOP, TIME_TR) )
		setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
		df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)		
	}
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 24-10-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v3<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	#	TODO is IDTR integer?
	# 	compute prevalence and incidence by year
	#	sampling model
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
									sexactive	<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
									infected	<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
									infdead		<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
									list(POP=length(sexactive), PREV=length(infected), PREVDIED=length(infdead))				
								},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	set(df.epi, NULL, 'GROWTHr', c(NA_real_, df.epi[, diff(log(PREV))]))	
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	s.PREV.base	<- exp(1)
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( pipeline.args['s.PREV.max',][, as.numeric(v)]/pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) * pipeline.args['s.PREV.min',][, as.numeric(v)] ]	
	#tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	#tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])
	stopifnot(df.sample[, all(s.n.TOTAL>=0)])
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences to sample=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences to sample=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- subset(df.trm, YR>= pipeline.args['yr.start',][, as.numeric(v)])
	tmp		<- tmp[, {
							z	<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
							z	<- sample(seq_along(IDREC), z)
							list( 	IDPOP=IDREC[z], TIME_TR=TIME_TR[z], 
									TIME_SEQ=TIME_TR[z]+rexp(length(z), rate=1/(3*30))/365, 
									INCIDENT_SEQ=rep('Y',length(z) ) )
						}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)
	tmp		<- df.inds[, which(TIME_SEQ>DOD)]
	set(df.inds, tmp, 'TIME_SEQ', df.inds[tmp, TIME_TR+(DOD-TIME_TR)/2])
	cat(paste('\n total number of incident sequences sampled=', df.inds[, length(which(!is.na(TIME_SEQ)))] ))	
	#	sample non-incident cases by year
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr+1 & yr-TIME_TR>1)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp		<- subset(tmp, yr-TIME_TR<12)
		cat(paste('\navailable with TIME_TR at most 12 years earlier, n=',nrow(tmp)))		
		tmp		<- subset(tmp, is.na(T1_SEQ) | (T1_SEQ-0.5<=yr & T1_SEQ+0.5>yr & T1_SEQ-0.5-TIME_TR>1))
		cat(paste('\navailable with T1_SEQ +- 6mo, n=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		stopifnot(tmp2<=nrow(tmp))
		tmp2	<- sample(seq_len(nrow(tmp)), tmp2)
		#	set variables in df.inds		
		tmp		<- data.table(IDPOP= tmp[tmp2, IDPOP], TIME_SEQ=runif(length(tmp2), min=yr, max=yr+1), INCIDENT_SEQ=rep('N',length(tmp2) ))
		cat(paste('\nsampled non-incident cases in year=',nrow(tmp)))
		tmp2	<- sapply(tmp[,IDPOP], function(x) df.inds[,which(IDPOP==x)])
		set(df.inds, tmp2, 'TIME_SEQ', tmp[,TIME_SEQ])
		set(df.inds, tmp2, 'INCIDENT_SEQ', tmp[,INCIDENT_SEQ])		
	}
	stopifnot( df.inds[, !any(TIME_SEQ>DOD, na.rm=TRUE)] )
	stopifnot( df.inds[, !any(TIME_SEQ<TIME_TR, na.rm=TRUE)] )	
	cat(paste('\n total number of HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR)))))
	cat(paste('\n total number of sampled HIV+ in df.inds=', nrow(subset(df.inds, !is.na(TIME_TR) & !is.na(TIME_SEQ)))))	
	#
	#	interpolate CD4 count at time of sampling
	#
	tmp		<- subset(df.inds, !is.na(TIME_SEQ), select=c(IDPOP, TIME_TR, INCIDENT_SEQ, TIME_SEQ, T1_CD4_500, T1_CD4_350, T1_CD4_200, DOD))
	set(tmp, NULL, 'TIME_TR', tmp[, TIME_TR-0.05])
	set(tmp, NULL, 'DOD', tmp[, DOD+0.05])
	setnames(tmp, c('TIME_TR','DOD'), c('T1_CD4_1000','T1_CD4_100'))
	tmp		<- melt(tmp, id.vars=c('IDPOP','TIME_SEQ','INCIDENT_SEQ'), value.name='T1_CD4', variable.name='CD4')
	set(tmp, NULL, 'CD4', tmp[, as.numeric(sapply(strsplit(as.character(CD4),'_'),'[[',3)) ])
	tmp2	<- tmp[, which(CD4==1000)]
	set(tmp, tmp2, 'CD4', runif(length(tmp2), min=700, max=1000))
	#tmp		<- tmp[, list(TIME_SEQ=TIME_SEQ[1], CD4_SEQ=exp(approx(T1_CD4, log(CD4), xout=TIME_SEQ[1], rule=1)$y)), by='IDPOP']
	tmp		<- tmp[, list(TIME_SEQ=TIME_SEQ[1], CD4_SEQ=round(approx(T1_CD4, CD4, xout=TIME_SEQ[1], rule=1)$y)), by='IDPOP']
	stopifnot( tmp[, !any(is.na(CD4_SEQ))] )
	df.inds	<- merge(df.inds, tmp, all.x=TRUE, by=c('IDPOP','TIME_SEQ'))
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
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#	add IDCLU to df.inds
	tmp			<- subset( df.trms, select=c(IDREC, IDTR, IDCLU) )
	tmp			<- subset( melt(tmp, id.var='IDCLU', value.name='IDPOP'), select=c(IDPOP, IDCLU))
	setkey(tmp, IDPOP, IDCLU)
	tmp			<- unique(tmp)
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 'POP', df.sample[, as.numeric(POP)])
		set(df.sample, NULL, 'PREV', df.sample[, as.numeric(PREV)])
		set(df.sample, NULL, 'INC', df.sample[, as.numeric(INC)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC','GROWTHr'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Proportion incident\nfrom acute in study', 'Proportion incident\nimported', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced','Annual growth rate'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC','GROWTHr'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,2020,2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=16)	
		#	plot CD4 counts at time of sampling
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=TIME_SEQ, y=CD4_SEQ, colour=INCIDENT_SEQ)) + geom_point() + geom_smooth() + scale_y_continuous(breaks=seq(100,1000,100))		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)					
		ggplot(subset(df.inds, !is.na(TIME_SEQ)), aes(x=floor(TIME_SEQ), fill=cut(CD4_SEQ, breaks=c(0,200,350,500,700,2000)))) + geom_bar(binwidth=1, position='fill') + 
				labs(fill='CD4 category', y='percent', x='year sequenced') + scale_y_continuous(breaks=seq(0,1,0.1)) + scale_fill_brewer(palette='Set1')		
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_CD4PROP.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=8, h=8)			
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution of #recipients per infector
		tmp	<- df.trms[, list(N= length(IDREC)), by='IDTR']
		ggplot(tmp, aes(x=N)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='recipients per source case\n(number)', breaks=seq(1,100,1)+0.5, label=seq(1,100,1)) +
				scale_y_continuous(breaks=seq(0,1,0.05))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_RecSource.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot time to death for infected
		tmp	<- subset(df.inds, !is.na(TIME_TR) & IDPOP>0 & DOD<max(DOD, na.rm=1), select=c(TIME_TR, DOD))
		ggplot(tmp, aes(x=DOD-TIME_TR)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='time to death for HIV+\n(years)', breaks=seq(1,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_T2DeathForInf.pdf',sep='')
		ggsave(file=file, w=8, h=8)		
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu & IDPOP>=0, select=c(IDPOP, GENDER, TIME_SEQ))
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
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
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
	if('RISK'%in%colnames(df.inds))
		df.inds[, RISK:=NULL]
	if('INCIDENT_SEQ'%in%colnames(df.inds))
		df.inds[, INCIDENT_SEQ:=NULL]
	if('TIME_SEQYR'%in%colnames(df.inds))
		df.inds[, TIME_SEQYR:=NULL]	
	if('TR_ACUTE'%in%colnames(df.trms))
		df.trms[, TR_ACUTE:=NULL]
	if('YR'%in%colnames(df.trms))
		df.trms[, YR:=NULL]	
	#	add columns that the virus tree simulator needs
	if(!'IDTR_TIME_INFECTED'%in%colnames(df.trms))
	{
		tmp		<- subset( df.inds, !is.na(TIME_TR), c(IDPOP, TIME_TR) )
		setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
		df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)		
	}
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 11-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.Seqsampler.v1<- function(df.ind, df.trm, pipeline.args, outfile.ind, outfile.trm, with.plot=1)
{	
	#	TODO is IDTR integer?
	# 	compute prevalence and incidence by year
	#	sampling model
	epi.adult	<- 13
	suppressWarnings( df.trm[, YR:= df.trm[, floor(TIME_TR)]] )
	df.epi		<- df.trm[, list(INC=length(IDREC), INC_ACUTE=length(which(TR_ACUTE=='Yes')),IMPORT=length(which(IDTR<0))), by='YR']
	tmp			<- df.epi[, 	{
				sexactive	<- which( floor(df.ind[['DOB']]+epi.adult)<=YR  &  ceiling(df.ind[['DOD']])>YR )
				infected	<- which( floor(df.ind[['DOB']])<=YR  &  floor(df.ind[['DOD']])>YR  &  floor(df.ind[['TIME_TR']])<=YR )
				infdead		<- which( floor(df.ind[['DOD']])==YR  &  floor(df.ind[['TIME_TR']])<=YR )
				list(POP=length(sexactive), PREV=length(infected), PREVDIED=length(infdead))				
			},by='YR']
	df.epi		<- merge( tmp, df.epi, by='YR' )	
	set(df.epi, NULL, 'PREVp', df.epi[, PREV/(POP-PREV)])	
	set(df.epi, NULL, 'INCp', df.epi[, INC/(POP-PREV)])
	set(df.epi, NULL, 'IMPORTp', df.epi[, IMPORT/INC])
	set(df.epi, NULL, 'ACUTEp', df.epi[, INC_ACUTE/INC])
	# 	SAMPLING PROBABILITIES and TOTALS PER YEAR
	#
	#	Can we detect a 25% or 50% reduction in HIV incidence in the most recent 2 or 3 years 
	#	with 1%, 5%, 10% of all recent incident cases sampled?
	#
	#	suppose exponentially increasing sampling over time
	#	the number of incident cases sampled is the total sampled in that year * the proportion of incident cases out of all non-sampled cases to date
	#	TODO this needs to be changed to fix the proportion of sequences sampled from incident
	s.PREV.base	<- exp(1)
	df.sample	<- subset( df.epi, YR>= pipeline.args['yr.start',][, as.numeric(v)] & YR<pipeline.args['yr.end',][, as.numeric(v)] )
	#	exponential rate of increasing s.TOTAL (total sampling rate) per year
	tmp			<- log( pipeline.args['s.PREV.max',][, as.numeric(v)]/pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) * pipeline.args['s.PREV.min',][, as.numeric(v)] ]	
	#tmp			<- log( 1+pipeline.args['s.PREV.max',][, as.numeric(v)]-pipeline.args['s.PREV.min',][, as.numeric(v)], base=s.PREV.base ) / df.sample[, diff(range(YR))]
	#tmp			<- df.sample[, s.PREV.base^( tmp*(YR-min(YR)) ) - 1 + pipeline.args['s.PREV.min',][, as.numeric(v)] ]
	set(df.sample, NULL, 's.CUMTOTAL', tmp)		
	set(df.sample, NULL, 's.n.CUMTOTAL', df.sample[, round(PREV*s.CUMTOTAL)])
	set(df.sample, NULL, 's.n.TOTAL', c(df.sample[1, s.n.CUMTOTAL], df.sample[, diff(s.n.CUMTOTAL)]))	
	set(df.sample, NULL, 's.n.INC', df.sample[, round(INC/(PREV-s.n.CUMTOTAL) * s.n.TOTAL)])
	set(df.sample, NULL, 's.n.notINC', df.sample[, round(s.n.TOTAL-s.n.INC)])
	stopifnot(df.sample[, all(s.n.TOTAL>0)])
	cat(paste('\n total number of sequences sampled=', df.sample[, sum( s.n.TOTAL )]))
	cat(paste('\n prop of sequences sampled among HIV+=', df.sample[, sum( s.n.TOTAL )] / df.sample[, rev(PREV)[1]]))		
	cat(paste('\n total number of incident sequences to sample=', df.sample[, sum( s.n.INC )]))
	cat(paste('\n total number of non-incident sequences to sample=', df.sample[, sum( s.n.notINC )]))	
	#	SAMPLE INFECTED INDIVIDUALS BASED ON NUMBERS PER YEAR
	#
	#	sample incident cases by year
	df.inds	<- copy(df.ind)
	tmp		<- subset(df.trm, YR>= pipeline.args['yr.start',][, as.numeric(v)])
	tmp		<- tmp[, {
				z	<- df.sample[['s.n.INC']][ which(df.sample[['YR']]==YR) ]
				z	<- sample(seq_along(IDREC), z)
				list( 	IDPOP=IDREC[z], TIME_TR=TIME_TR[z], 
						TIME_SEQ=TIME_TR[z]+rexp(length(z), rate=1/(3*30))/365, 
						INCIDENT_SEQ=rep('Y',length(z) ) )
			}, by='YR']
	df.inds	<- merge(df.inds, subset(tmp, select=c(IDPOP, TIME_SEQ, INCIDENT_SEQ)), by='IDPOP', all.x=1)
	cat(paste('\n total number of incident sequences sampled=', df.inds[, length(which(!is.na(TIME_SEQ)))] ))	
	#	sample non-incident cases by year
	for(yr in df.sample[, YR])		#TODO took out [-1] because there are s.n.notINC for DSPS in 1980
	{
		#	of all infected and not incident and not yet sampled, sample
		cat(paste('\nadd non-incident samples in year',yr))
		tmp		<- subset(df.inds, is.na(TIME_SEQ) & floor(DOB)<=yr & ceiling(DOD)>yr & floor(TIME_TR)<yr)
		cat(paste('\navailable non-sampled non-incident cases in year=',nrow(tmp)))
		tmp		<- subset(tmp, T1_SEQ-0.5<=yr & T1_SEQ+0.5>yr & T1_SEQ-0.5-TIME_TR>1)
		cat(paste('\navailable with T1_SEQ +- 6mo, n=',nrow(tmp)))
		tmp2	<- df.sample[['s.n.notINC']][ which(df.sample[['YR']]==yr) ]
		stopifnot(tmp2<=nrow(tmp))
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
	#	cluster with index case
	tmp			<- subset(df.trms, select=c(IDTR, IDREC))			
	tmp			<- graph.data.frame(tmp, directed=TRUE, vertices=NULL)
	tmp			<- data.table(IDPOP=as.integer(V(tmp)$name), CLU=clusters(tmp, mode="weak")$membership)
	tmp2		<- tmp[, list(CLU_SIZE=length(IDPOP)), by='CLU']
	setkey(tmp2, CLU_SIZE)
	tmp2[, IDCLU:=rev(seq_len(nrow(tmp2)))]
	tmp			<- subset( merge(tmp, tmp2, by='CLU'), select=c(IDPOP, IDCLU) )
	setnames(tmp, 'IDPOP', 'IDREC')
	df.trms		<- merge( df.trms, tmp, by='IDREC', all.x=TRUE )
	stopifnot( nrow(subset(df.trms, is.na(IDCLU)))==0 )
	cat(paste('\nFound transmission clusters, n=', df.trms[, length(unique(IDCLU))]))
	#	add IDCLU to df.inds
	tmp			<- subset( df.trms, select=c(IDREC, IDTR, IDCLU) )
	tmp			<- subset( melt(tmp, id.var='IDCLU', value.name='IDPOP'), select=c(IDPOP, IDCLU))
	setkey(tmp, IDPOP, IDCLU)
	tmp			<- unique(tmp)
	df.inds		<- merge( df.inds, tmp, by='IDPOP', all.x=TRUE )
	#
	#	PLOTS
	#
	if(with.plot)
	{
		require(ggplot2)
		require(reshape2)
		#	plot numbers sampled, prevalent, incident
		set(df.sample, NULL, 'POP', df.sample[, as.real(POP)])
		set(df.sample, NULL, 'PREV', df.sample[, as.real(PREV)])
		set(df.sample, NULL, 'INC', df.sample[, as.real(INC)])
		tmp	<- data.table(	stat=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC'), 
				stat.long=c('population size','HIV infected', 'HIV incident', 'Proportion incident\nfrom acute in study', 'Proportion incident\nimported', 'Total\nsequenced', 'Total\nincident\nsequenced', 'Total\nnon-incident\nsequenced'))
		tmp	<- merge(	melt(df.sample, id.vars='YR', measure.vars=c('POP','PREV','INC','ACUTEp','IMPORTp','s.n.TOTAL','s.n.INC','s.n.notINC'), variable.name='stat', value.name='v'),
				tmp, by='stat' )
		ggplot(tmp, aes(x=YR, y=v, group=stat.long)) + geom_point() +
				scale_x_continuous(name='year', breaks=seq(1980,2020,2)) + scale_y_continuous(name='total')	+
				facet_grid(stat.long ~ ., scales='free_y', margins=FALSE)
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Totals.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		ggsave(file=file, w=16, h=10)		 
		#	plot distribution between transmission time and sequencing time
		tmp	<- subset(df.inds, !is.na(TIME_SEQ))
		set(tmp, NULL, 'TIME_TO_SEQ', tmp[, TIME_SEQ-TIME_TR])
		ggplot(tmp, aes(x=TIME_TO_SEQ)) + geom_histogram(binwidth=1) + 
				scale_x_continuous(name='time from transmission to sequence sampling\n(years)', breaks=seq(0,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_Time2Seq.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot distribution of #recipients per infector
		tmp	<- df.trms[, list(N= length(IDREC)), by='IDTR']
		ggplot(tmp, aes(x=N)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='recipients per source case\n(number)', breaks=seq(1,100,1)+0.5, label=seq(1,100,1)) +
				scale_y_continuous(breaks=seq(0,1,0.05))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_RecSource.pdf',sep='')
		ggsave(file=file, w=8, h=8)
		#	plot time to death for infected
		tmp	<- subset(df.inds, !is.na(TIME_TR) & IDPOP>0 & DOD<max(DOD, na.rm=1), select=c(TIME_TR, DOD))
		ggplot(tmp, aes(x=DOD-TIME_TR)) + geom_histogram(binwidth=1, aes(y= ..density..)) +
				scale_x_continuous(name='time to death for HIV+\n(years)', breaks=seq(1,100,2))
		file<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_T2DeathForInf.pdf',sep='')
		ggsave(file=file, w=8, h=8)		
		#	plot transmission network
		file		<- paste(substr(outfile.ind, 1, nchar(outfile.ind)-7),'INFO_TrNetworks.pdf',sep='')
		cat(paste('\nPlotting to file',file))
		pdf(file=file, w=20, h=20)
		dummy	<- sapply( df.inds[, sort(na.omit(unique(IDCLU)))], function(clu)
				{
					cat(paste('\nprocess cluster no',clu))
					tmp					<- subset(df.inds, IDCLU==clu & IDPOP>=0, select=c(IDPOP, GENDER, TIME_SEQ))
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
		#ggplot(df.trms, aes(x=IDTR, y=TIME_TR)) + geom_point()
		#ggplot(df.trms, aes(x=IDTR, y=IDCLU)) + geom_point()
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
	if('RISK'%in%colnames(df.inds))
		df.inds[, RISK:=NULL]
	if('INCIDENT_SEQ'%in%colnames(df.inds))
		df.inds[, INCIDENT_SEQ:=NULL]
	if('TIME_SEQYR'%in%colnames(df.inds))
		df.inds[, TIME_SEQYR:=NULL]	
	if('TR_ACUTE'%in%colnames(df.trms))
		df.trms[, TR_ACUTE:=NULL]
	if('YR'%in%colnames(df.trms))
		df.trms[, YR:=NULL]	
	#	add columns that the virus tree simulator needs
	if(!'IDTR_TIME_INFECTED'%in%colnames(df.trms))
	{
		tmp		<- subset( df.inds, !is.na(TIME_TR), c(IDPOP, TIME_TR) )
		setnames(tmp, c('IDPOP','TIME_TR'), c('IDTR','IDTR_TIME_INFECTED') )
		df.trms	<- merge(df.trms, tmp, by='IDTR', all.x=TRUE)		
	}
	cat(paste('\nwrite to file',outfile.ind))
	write.csv(file=outfile.ind, df.inds)
	cat(paste('\nwrite to file',outfile.trm))
	write.csv(file=outfile.trm, df.trms)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 20-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create starting sequence sampler
#' @description Returns a function and function arguments to draw ancestral sequences. 
#' @param root.ctime.grace	Sample a starting sequence with time that matches a query times +- this grace
#' @param sample.grace		Internal parameter to make sure the requested number of samples is obtained. Internally oversample by this multiplier to the sample size, and then check if sequences are unique.
#' @return list of the sampler \code{rANCSEQ} and its arguments \code{rANCSEQ.args}
#' @export
PANGEA.RootSeq.create.sampler<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp				<- lapply(seq_along(root.ctime), function(i)
						{
							tmp	<- subset(rANCSEQ.args$anc.seq.info, 	CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace 	)
							if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
							{
								warning( paste('\nFor root',i,': number of samples is n=',nrow(tmp),'. safe pool size is n=',rANCSEQ.args$sample.grace*100) )
							}					
							data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
						})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 14-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v3<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
					{
						warning( paste('\nFor root',i,': number of samples is n=',nrow(tmp),'. safe pool size is n=',rANCSEQ.args$sample.grace*100) )
					}					
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v2<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	#tree.id.labelsep		<- '|'
	#tree.id.labelidx.ctime	<- 4
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args	<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
					{
						warning( paste('\nFor root',i,': number of samples is n=',nrow(tmp),'. safe pool size is n=',rANCSEQ.args$sample.grace*100) )
					}
					#print( paste('\n',nrow(tmp),'\t',rANCSEQ.args$sample.grace*100) )
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		#	set new calendar time for sequences
		#set(tmp, NULL, 'LABEL_NEW', tmp[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		#set(tmp, NULL, 'LABEL_NEW', tmp[, LABEL_NEW+sample.shift])	
		#tmp			<- tmp[,	{
		#							z							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
		#							z[tree.id.labelidx.ctime]	<- LABEL_NEW
		#							list(LABEL_NEW=paste(z, collapse=tree.id.labelsep,sep=''))
		#						}, by='LABEL']
		#setkey(tmp, LABEL)
		#rownames(anc.seq.draw)		<- tmp[ rownames(anc.seq.draw), ][, LABEL_NEW]		
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return ancestral sequence sampler	
#	olli originally written 22-08-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.RootSeq.create.sampler.v1<- function(root.ctime.grace= 0.5, sample.grace= 3, sample.shift= 40)
{	
	#tree.id.labelsep		<- '|'
	#tree.id.labelidx.ctime	<- 4
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140811_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, sample.shift=sample.shift, 
			anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{				
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME+rANCSEQ.args$sample.shift>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME+rANCSEQ.args$sample.shift<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<=rANCSEQ.args$sample.grace*100)
					{
						print(c(nrow, rANCSEQ.args$sample.grace*100))
						stop()
					}
					data.table( LABEL= tmp[, sample( LABEL, rANCSEQ.args$sample.grace ) ], CALENDAR_TIME=root.ctime[i], DRAW=i )
				})
		tmp		<- do.call('rbind',tmp)
		#	get unique seqs
		setkey(tmp, LABEL)
		tmp				<- unique(tmp)	
		anc.seq.draw	<- do.call( 'cbind', list( rANCSEQ.args$anc.seq.gag[tmp[, LABEL], ], rANCSEQ.args$anc.seq.pol[tmp[, LABEL], ], rANCSEQ.args$anc.seq.env[tmp[, LABEL], ] ) ) 
		anc.seq.draw	<- seq.unique(anc.seq.draw)
		#	check if at least one seq for each draw
		tmp				<- merge( data.table(LABEL=rownames(anc.seq.draw)), tmp, by='LABEL' )	
		stopifnot( !length(setdiff( seq_along(root.ctime), tmp[, unique(DRAW)] )) )
		#	take first seq for each draw	
		tmp				<- tmp[, list(LABEL= LABEL[1]), by='DRAW']
		setkey(tmp, DRAW)
		anc.seq.draw	<- anc.seq.draw[ tmp[, LABEL], ]
		#	set new calendar time for sequences
		#set(tmp, NULL, 'LABEL_NEW', tmp[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		#set(tmp, NULL, 'LABEL_NEW', tmp[, LABEL_NEW+sample.shift])	
		#tmp			<- tmp[,	{
		#							z							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
		#							z[tree.id.labelidx.ctime]	<- LABEL_NEW
		#							list(LABEL_NEW=paste(z, collapse=tree.id.labelsep,sep=''))
		#						}, by='LABEL']
		#setkey(tmp, LABEL)
		#rownames(anc.seq.draw)		<- tmp[ rownames(anc.seq.draw), ][, LABEL_NEW]		
		anc.seq.draw
	}
	list(rANCSEQ=rANCSEQ, rANCSEQ.args=rANCSEQ.args)
}
##--------------------------------------------------------------------------------------------------------
#	return evolutionary rate sampler for transmission edges	
#	olli originally written 16-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create sampler of Evolutionary Rates for transmission edges
#' @description Returns a function to draw evolutionary rates for transmission edges. Currently modelled with a log normal density.
#' @param er.shift		shift to mean of the log normal density
#' @return R function
#' @export
PANGEA.TransmissionEdgeEvolutionaryRate.create.sampler<- function(er.meanlog, er.sdlog)
{
	if(0)
	{
		x		<- seq(0.0005, 0.003, 0.00001)		
		tmp		<- data.table(x=x, y1=dLOGNO(x, mu=-6.714086, sigma=0.15), y2=dLOGNO(x, mu=-6.714086-0.15^2/2+0.075^2/2, sigma=0.075), y4=dLOGNO(x, mu=-6.714086-0.15^2/2+0.0375^2/2, sigma=0.0375))
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.003,0.0005))
		x		<- seq(0.0005, 0.003, 0.00001)
		tmp		<- data.table(x=x, y1=dLOGNO(x, mu=log(0.002239075)-0.05^2/2, sigma=0.05), y5=dLOGNO(x, mu=log(0.002239075-0.0005)-0.065^2/2, sigma=0.065), y10=dLOGNO(x, mu=log(0.002239075-0.001)-0.09^2/2, sigma=0.09))
		
		x		<- seq(0.0005, 0.01, 0.00001)
		tmp		<- data.table(	x=x, 
								y5=dLOGNO(x, mu=log(0.002239075)-0.05^2/2, sigma=0.05), y7=dLOGNO(x, mu=log(0.002239075)-0.07^2/2, sigma=0.07), y10=dLOGNO(x, mu=log(0.002239075)-0.1^2/2, sigma=0.1), y13=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13), 
								W441=dLOGNO(x, mu=log(0.00447743)-0.3^2/2, sigma=0.3), W442=dLOGNO(x, mu=log(0.00447743)-0.4^2/2, sigma=0.4), W443=dLOGNO(x, mu=log(0.00447743)-0.5^2/2, sigma=0.5))
		#	times 2				
		x		<- seq(0.0005, 0.01, 0.00001)
		tmp		<- data.table(	x=x, 
								TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13),
								TransmissionLineage2=dLOGNO(x, mu=log(0.002239075)-0.3^2/2, sigma=0.3),
								TransmissionLineage3=dLOGNO(x, mu=log(0.002239075)-0.2^2/2, sigma=0.2),
								TipBranch=dLOGNO(x, mu=log(0.00447743)-0.5^2/2, sigma=0.5)
								)
		#	times 1.5				
		#tmp		<- data.table(	x=x, 
		#						TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.13^2/2, sigma=0.13), 
		#						TipBranch=dLOGNO(x, mu=log(0.003358613)-0.3^2/2, sigma=0.3))						
		#tmp		<- data.table(	x=x, 
		#						TransmissionLineage=dLOGNO(x, mu=log(0.002239075)-0.01^2/2, sigma=0.01), 
		#						TipBranch=dLOGNO(x, mu=log(0.003)-0.2^2/2, sigma=0.2))
						
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.02,0.001))		
	}
	if(0)
	{
		require(gamlss)		
		file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140907_n390_BEASTlog.R')	
		cat(paste('\nreading GTR parameters from file',file))
		load(file)	# expect log.df
		#	exclude odd BEAST runs
		log.df		<- subset(log.df, !(GENE=='ENV' & FILE=='pool3'))
		#	exclude cols
		log.df[, ucld.mean:=NULL]
		log.df[, ucld.stdev:=NULL]
		log.df[, coefficientOfVariation:=NULL]
		log.df[, treeModel.rootHeight:=NULL]
		#	get posterior distribution of overall meanRate across genes
		tmp		<- log.df[, list(meanRate=mean(meanRate)), by='GENE'][, mean(meanRate)]	
		log.df	<- log.df[, list(meanRate=meanRate/mean(meanRate)*tmp), by='GENE']
		#ggplot(log.df, aes(x=meanRate)) + geom_histogram() + facet_grid(.~GENE, scales='free')
		#	get lognormal density parameters
		sd.log	<- log.df[, sd(log(meanRate))]			# 0.1499262
		mean.log<- log.df[, mean(log(meanRate))]		# -6.113092; 0.002239075
		#	mean is exp( mean.log+sd.log*sd.log/2 ); median is exp(mean.log) 
	}
	rERbw.args	<- list(meanlog=er.meanlog, sdlog=er.sdlog)
	if(er.sdlog>0)
	{
		rERbw		<- function(n, rERbw.args)
				{	
					rlnorm(n, meanlog=rERbw.args$meanlog, sdlog=rERbw.args$sdlog)		
				}		
	}
	if(er.sdlog==0)
	{
		rERbw		<- function(n, rERbw.args)
				{	
					rep( exp(rERbw.args$meanlog), n)		
				}		
	}	
	list(rERbw=rERbw, rERbw.args=rERbw.args)	
}
##--------------------------------------------------------------------------------------------------------
#	return evolutionary rate modifier sampler for transmission edges	
#	olli originally written 21-08-2014
##--------------------------------------------------------------------------------------------------------
PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1<- function(bwerm.mu=1.5, bwerm.sigma=0.12)
{
	require(gamlss)
	if(0)
	{
		#from Vrancken et al:
		#c(3.5/2.5, 7.0/4.2) #1.4, 1.67	draw this from lognormal with mean 1.5 and variance so that tail captures 1.8		
		x		<- seq(1.01, 15, 0.01)
		tmp		<- data.table(x=x, y5=dLOGNO(x, mu=log(1.5), sigma=0.12), y8=dLOGNO(x, mu=log(1.75), sigma=0.103), y2=dLOGNO(x, mu=log(2), sigma=0.09))
		tmp		<- data.table(x=x, y4=dLOGNO(x, mu=log(4), sigma=0.045), y3=dLOGNO(x, mu=log(3), sigma=0.06), y2=dLOGNO(x, mu=log(2), sigma=0.09))
		tmp		<- data.table(x=x, y4=dLOGNO(x, mu=log(4), sigma=0.06), y3=dLOGNO(x, mu=log(3), sigma=0.06), y2=dLOGNO(x, mu=log(2), sigma=0.06))		
		tmp		<- data.table(x=x, y6=dLOGNO(x, mu=log(6)+0.4^2/2, sigma=0.4), y5=dLOGNO(x, mu=log(5)+0.4^2/2, sigma=0.4), y4=dLOGNO(x, mu=log(4)+0.4^2/2, sigma=0.4), y3=dLOGNO(x, mu=log(3)+0.4^2/2, sigma=0.4))
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_log10(breaks=seq(1,20,1))
	}
	rER.bwm<- function(n)
	{
		rLOGNO(n, mu=log(bwerm.mu), sigma=bwerm.sigma)
	}
	rER.bwm
}
##--------------------------------------------------------------------------------------------------------
#	return within host evolutionary rate sampler	
#	olli originally written 21-08-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create sampler of Within Host Evolutionary Rates 
#' @description Returns a function to draw within host evolutionary rates. Currently modelled with a log normal density.
#' @param wher.mu		mean of the log normal density
#' @param wher.sigma 	standard deviation of the log normal density
#' @return R function
#' @export
PANGEA.WithinHostEvolutionaryRate.create.sampler.v1<- function(wher.mu=log(0.005), wher.sigma=0.8)
{
	require(gamlss)
	if(0)
	{	
		#extremely basic model of within host evolutionary rate from HIV-1B pol estimates in the literature
		#median log10 ER of pol in Alizon / Fraser
		df.er	<- data.table(ER= 10^c(-1.85, -2.2, -2.5, -2.7, -2.72, -3.2), GENE='POL')		
		tmp		<- gamlss(ER~1, data=df.er, family=LOGNO)
		x		<- seq(0.0005, 0.02, 0.0001)
		tmp		<- data.table(x=x, y5=dLOGNO(x, mu=log(0.005), sigma=0.8), y4=dLOGNO(x, mu=log(0.004), sigma=0.7), y3=dLOGNO(x, mu=log(0.003), sigma=0.6))
		tmp		<- data.table(x=x, y83=dLOGNO(x, mu=-4.784295+0.045, sigma=0.3), y62=dLOGNO(x, mu=-5.071977+0.08, sigma=0.4), y41=dLOGNO(x, mu=-5.477443+0.18, sigma=0.6))
		#	should have been '-' all along:
		#	mean(rLOGNO(1e4, mu=-5.071977+0.4^2/2, sigma=0.4))
		tmp		<- data.table(x=x, y62=dLOGNO(x, mu=-5.071977+0.08, sigma=0.4), y62b=dLOGNO(x, mu=-5.071977-0.4^2/2, sigma=0.4) )
		tmp		<- data.table(x=x, y62=dLOGNO(x, mu=log(0.006716145)-0.37^2/2, sigma=0.37), y621=dLOGNO(x, mu=-5.071977-0.2^2/2, sigma=0.2), y441=dLOGNO(x, mu=log(0.00447743)-0.3^2/2, sigma=0.3) )		
		tmp		<- melt(tmp, id.var='x')
		ggplot(tmp, aes(x=x, y=value, group=variable, colour=variable)) + geom_line() + scale_x_continuous(breaks=seq(0,0.02,0.002))
		#correlation model: need high multiplier if ER high
		require(compositions)		
		y		<- rlnorm.rplus(1e4, meanlog=c(-5.071977+0.08, log(3)+0.4^2/2),  varlog=diag(c(0.6,0.4))%*%matrix(c(1,0.95,0.95,1),2,2)%*%diag(c(0.6,0.4)) )		
		tmp		<- data.table(wh=y[,1], bm=y[,2])
		ggplot(tmp, aes(x= wh, y=bm)) + geom_point()
		#explore marginal ER along branches
		tmp		<- do.call('rbind',lapply(c(3,4,5,6), function(bm)
				{
					tmp	<- rlnorm.rplus(1e4, meanlog=c(-5.071977+0.08, log(bm)+0.4^2/2),  varlog=diag(c(0.6,0.4))%*%matrix(c(1,0.95,0.95,1),2,2)%*%diag(c(0.6,0.4)) )
					data.table(value=tmp[,1]/tmp[,2], variable=paste('y',bm,sep=''))
				}))
		ggplot(tmp, aes(x=value, fill=variable)) + geom_histogram(binwidth=0.00025) + facet_grid(.~variable, scales='free')
	}	
	if(wher.sigma>0)
	{
		rER.pol<- function(n)
			{		
				ans	<- numeric(0)
				while(length(ans)<n)
				{
					tmp		<- rLOGNO(2*n, mu=wher.mu, sigma=wher.sigma)
					tmp[which(tmp>0.02)]	<- NA
					ans		<- c(ans, na.omit(tmp))						
				}			
				ans[seq_len(n)]
			}
	}
	if(wher.sigma==0)
	{
		rER.pol<- function(n)
			{		
				rep(exp(wher.mu), n)			
			}
	}
	rER.pol
}
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
#	olli originally written 09-09-2014
#	tree 		beast trees in ape format, needed to compute calendar time for each ancestral sequence
#	node.stat	data.table containing meta information in nexus file for nodes
#	return 		list of GAG POL ENV sequences in ape format 
PANGEA.RootSeqSim.get.ancestral.seq.pg<- function(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=5)
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
	set(node.stat, NULL, 'GENE', node.stat[, substr(STAT,1,nchar(STAT)-4)])
	set(node.stat, NULL, 'STAT', node.stat[, substr(STAT,nchar(STAT)-2,nchar(STAT))])
	#	reconstruct genes from codon positions
	node.stat	<- dcast.data.table(node.stat, TREE_ID + NODE_ID + GENE + CALENDAR_TIME ~ STAT, value.var="VALUE")
	node.stat	<- node.stat[, {
									tmp		<- do.call('rbind',sapply(list(CP1,CP2,CP3), strsplit, ''))
									tmp		<- paste(as.vector(tmp), collapse='')
									list(SEQ=tmp, GENE=GENE, CALENDAR_TIME=CALENDAR_TIME)
								}, by=c('TREE_ID','NODE_ID')]	
	set(node.stat, NULL, 'LABEL', node.stat[, paste(TREE_ID, NODE_ID, round(CALENDAR_TIME,d=3), sep='|')])		
	#	remove tree id STATE_xx where xx is smaller than burn-in
	set(node.stat, NULL, 'BEAST_MCMC_IT', node.stat[, as.integer(sapply(strsplit(TREE_ID,tree.id.sep),'[[',tree.id.idx.mcmcit))])
	node.stat			<- subset(node.stat, BEAST_MCMC_IT>tree.id.burnin)
	cat(paste('\nFound ancestral sequences, n=', nrow(node.stat)  ))
	tmp			<- node.stat[, unique(GENE)]
	stopifnot(length(tmp)==1)
	node.stat
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