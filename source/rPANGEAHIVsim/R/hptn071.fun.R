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
	ph<-  seq.read.newick(text=dummy.tree)
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
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create data.table of GTR parameters
#' @description Returns a data.table of GTR parameters. 
#' @return data.table
#' @export
PANGEA.GTR.params<- function()
{		
	file		<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_BEASTlog.R')	
	cat(paste('\nreading GTR parameters from file',file))
	load(file)	# expect log.df
	#	exclude odd BEAST runs
	log.df		<- subset(log.df, !(GENE=='GAG' & FILE=='pool1'))
	log.df		<- subset(log.df, !(GENE=='POL' & FILE=='pool2'))
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
#	return ancestral sequence sampler	
#	olli originally written 10-09-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create starting sequence sampler
#' @description Returns a function and function arguments to draw ancestral sequences. 
#' @param root.ctime.grace	Sample a starting sequence with time that matches a query times +- this grace
#' @param sample.grace		Internal parameter to make sure the requested number of samples is obtained. Internally oversample by this multiplier to the sample size, and then check if sequences are unique.
#' @return list of the sampler \code{rANCSEQ} and its arguments \code{rANCSEQ.args}
#' @export
PANGEA.RootSeq.create.sampler.v1<- function(root.ctime.grace= 0.5, sample.grace= 3)
{	
	#tree.id.labelsep		<- '|'
	#tree.id.labelidx.ctime	<- 4
	file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfgBwhRc-_140902_n390_AncSeq.R')
	cat(paste('\nLoading starting sequences from file', file))
	load(file)		#expect "anc.seq.gag"  "anc.seq.pol"  "anc.seq.env"  "anc.seq.info"
	setkey(anc.seq.info, CALENDAR_TIME)
	rANCSEQ.args<- list(	root.ctime.grace=root.ctime.grace, sample.grace=sample.grace, anc.seq.info=anc.seq.info, anc.seq.gag=anc.seq.gag, anc.seq.pol=anc.seq.pol, anc.seq.env=anc.seq.env)	
	
	rANCSEQ<- function(root.ctime, rANCSEQ.args)
	{		
		tmp		<- lapply(seq_along(root.ctime), function(i)
				{
					tmp	<- subset(rANCSEQ.args$anc.seq.info, CALENDAR_TIME>root.ctime[i]-rANCSEQ.args$root.ctime.grace &  CALENDAR_TIME<=root.ctime[i]+rANCSEQ.args$root.ctime.grace)
					if(nrow(tmp)<rANCSEQ.args$sample.grace*100)
					{
						cat(paste('\n',nrow(tmp),'\t',rANCSEQ.args$sample.grace*100))
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
					stopifnot(nrow(tmp)>rANCSEQ.args$sample.grace*100)
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
#	return within host evolutionary rate sampler	
#	olli originally written 21-08-2014
##--------------------------------------------------------------------------------------------------------
#' @title Create sampler of the Between Host Evolutionary Rate Multipliers
#' @description Returns a function to draw the between host evolutionary rate multipliers. Currently modelled with a log normal density.
#' @param bwerm.mu		mean of the log normal density
#' @param bwerm.sigma	standard deviation of the log normal density
#' @return R function
#' @export
PANGEA.BetweenHostEvolutionaryRateModifier.create.sampler.v1<- function(bwerm.mu=1.5, bwerm.sigma=0.12)
{
	require(gamlss)
	if(0)
	{
		#from Vrancken et al:
		#c(3.5/2.5, 7.0/4.2) #1.4, 1.67	draw this from lognormal with mean 1.5 and variance so that tail captures 1.8
		x		<- seq(1.01, 2, 0.01)
		plot(x, dLOGNO(x, mu=log(1.5), sigma=0.1), type='l')
		lines(x, dLOGNO(x, mu=log(1.5), sigma=0.07), col='blue' )
		lines(x, dLOGNO(x, mu=log(1.5), sigma=0.12), col='red' )		
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
PANGEA.WithinHostEvolutionaryRate.create.sampler.v1<- function(wher.mu=0.005, wher.sigma=0.8)
{
	require(gamlss)
	if(0)
	{	
		#extremely basic model of within host evolutionary rate from HIV-1B pol estimates in the literature
		#median log10 ER of pol in Alizon / Fraser
		df.er	<- data.table(ER= 10^c(-1.85, -2.2, -2.5, -2.7, -2.72, -3.2), GENE='POL')		
		tmp		<- gamlss(ER~1, data=df.er, family=LOGNO)
		x		<- seq(0.0005, 0.02, 0.0001)
		plot(x, dLOGNO(x, mu=coef(tmp, what='mu'), sigma=exp(coef(tmp, what='sigma'))), type='l')
		lines(x, dLOGNO(x, mu=log(0.005), sigma=0.8), col='blue')		
	}	
	rER.pol<- function(n)
				{		
					ans	<- numeric(0)
					while(length(ans)<n)
					{
						tmp		<- rLOGNO(2*n, mu=log(wher.mu), sigma=wher.sigma)
						tmp[which(tmp>0.02)]	<- NA
						ans		<- c(ans, na.omit(tmp))						
					}			
					ans[seq_len(n)]
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