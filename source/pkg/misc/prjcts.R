


prog.hello<- function()	
{
	print('hello')
}
##--------------------------------------------------------------------------------------------------------
##	run BEAST
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.run<- function()
{
	require(hivclust)
	
	DATA		<<- "/work/or105/Gates_2014"
	#indir		<- paste(DATA,'methods_comparison_rootseqsim/140727',sep='/')
	#infile		<- 'ALLv02.n100.rlx.gmrf' 
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140729',sep='/')
	infile		<- 'ALLv04.n97.rlx.gmrf' 
	insignat	<- 'Sun_Jul_27_09-00-00_2014'
	cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median", hpc.tmpdir.prefix="beast", hpc.ncpu=1)
	cat(cmd)
	
	cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=1, hpc.walltime=91, hpc.mem="1800mb")
	outdir		<- indir
	outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
	hivc.cmd.hpccaller(outdir, outfile, cmd)
}
##--------------------------------------------------------------------------------------------------------
##	get sequences from BEAST XML, create concatenated file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.getancestralseq.from.output<- function()
{
	#dir.name	<- "/work/or105/Gates_2014"
	dir.name	<- '/Users/Oliver/duke/2014_Gates'  
	
	if(1)	#	devel
	{		
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140730',sep='/')
		outdir				<- indir
		infile.xml			<- 'working.xml'
		infile.beastparsed	<- 'working.R'
		outfile				<- 'working_ancseq.R'
		
		#	load BEAST PARSER output
		file		<- paste(indir, '/',infile.beastparsed, sep='')
		load(file)	#	expect tree, node.stat		
		#	get original sequences
		file		<- paste(indir, '/',infile.xml, sep='')
		bxml		<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
		bseq		<- hivc.beast.get.sequences(bxml, verbose=1)	
		bseq		<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')				
		#	compute gag pol env		
		tmp			<- PANGEA.RootSeqSim.get.ancestral.seq(tree, node.stat, bseq, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=1, label.sep='|', label.idx.ctime=2)
		ancseq.gag	<- tmp$GAG
		ancseq.env	<- tmp$ENV
		ancseq.pol	<- tmp$POL
		#	save as R
		file		<- paste(outdir, outfile, sep='/')
		save(ancseq.gag, ancseq.env, ancseq.pol, file=file)			
		#	save as FASTA
		file		<- paste(outdir, paste(substr(outfile,1,nchar(outfile)-1),'fasta',sep=''), sep='/')		
		write.dna(cbind(ancseq.gag, ancseq.env, ancseq.pol), file, format = "fasta")		
		#
		#	sample ancestral sequences between 1980-2000 and reconstruct tree with RAxML
		#
		ancseq				<- cbind(ancseq.gag, ancseq.env, ancseq.pol)
		label.sep			<- '|'
		label.idx.tree.id	<- 1
		label.idx.node.id	<- 2
		label.idx.ctime		<- 3
		ancseq.label		<- data.table(	TREE_ID= sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.tree.id),
				NODE_ID= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.node.id)),
				CALENDAR_TIME= as.numeric(sapply(strsplit(rownames(ancseq), label.sep, fixed=1),'[[',label.idx.ctime)))
		hist( ancseq.label[, CALENDAR_TIME], breaks=100 )
	}
	if(0)
	{
		
		tmp	<- subset(node.stat, select=c(TREE_ID, NODE_ID, CALENDAR_TIME))
		setkey(tmp, TREE_ID, NODE_ID)
		hist(tmp[, CALENDAR_TIME], breaks=seq(1930,2011,1))
		#	reconstruct gene sequences and store in ape format
		#	get calendar time for gene sequence
		set(node.stat, NULL, 'VALUE', node.stat[, gsub('\"','',VALUE)])
		#	check seq lengths
		tmp		<- node.stat[, list(NCHAR=nchar(VALUE)), by=c('STAT','NODE_ID','TREE_ID')]
		stopifnot( tmp[, list(CHECK=all(NCHAR[1]==NCHAR)), by='STAT'][, all(CHECK)] )
		
		tmp[, list(CHECK=unique(NCHAR)), by='STAT']
		
		ENV.CP1<- "AGA"
		ENV.CP2<- "XYZ"
		ENV.CP3<- "KLM"
		tmp		<- do.call('rbind',sapply(list(ENV.CP1,ENV.CP2,ENV.CP3), strsplit, ''))
		tmp		<- paste(as.vector(tmp), collapse='')
		
		subset(node.stat, TREE_ID=='STATE_0' & NODE_ID==which(btree[[1]]$tip.label=='C.BW.AF443074|1996'))
		
	}
	if(0)	#devel
	{
		dir.name			<- '/Users/Oliver/duke/2014_Gates'  
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140801',sep='/')
		infile				<- 'ALLv06.n97.rlx.gmrf' 		
		insignat			<- 'Sun_Jul_27_09-00-00_2014'	
		file				<- paste(indir, '/', infile, '_', insignat, '.timetrees', sep='')
		
		indir				<- paste(dir.name,'methods_comparison_rootseqsim/140730',sep='/')		
		infile.timetrees	<- 'working.timetrees'
		outdir				<- indir
		
		file				<- paste(indir, '/', infile.timetrees, sep='')
		tmp					<- hivc.beast2out.read.nexus.and.stats(file, tree.id=NA, method.node.stat='any.node')
		tree				<- tmp$tree
		node.stat			<- tmp$node.stat
		
		file				<- paste(indir,'/',paste(substr(infile.timetrees,1,nchar(infile.timetrees)-9),'R',sep=''),sep='')
		save(tree, node.stat, file=file)
	}	
}
##--------------------------------------------------------------------------------------------------------
##	get sequences from BEAST XML, create concatenated file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.checkDrugResistanceMask<- function()
{
	require(XML)
	require(ape)
	DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140727',sep='/')
	outdir		<- paste(DATA,'methods_comparison_rootseqsim/140728',sep='/')
	infile		<- 'ALLv02.n100.rlx.gmrf' 
	insignat	<- 'Sun_Jul_27_09-00-00_2014'
	
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140727/ALLv01.n100.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140801/ALLv06.n97.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	bxml			<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
	bseq			<- hivc.beast.get.sequences(bxml, verbose=1)	
	bseq			<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')
	#	check all seq of same length per gene
	bseq			<- merge(bseq, bseq[, list(SEQ_N=nchar(SEQ)), by=c('GENE','TAXON_ID')], by=c('GENE','TAXON_ID'))
	stopifnot( bseq[, all(SEQ_N==SEQ_N[1]), by='GENE'][, all(V1)] )
	#	check 3 genes per taxon
	stopifnot( bseq[, length(SEQ)==3, by='TAXON_ID'][, all(V1)] )
	#	check if indeed patterns are compressed
	if(0)
	{
		bseq	<- subset(bseq, GENE=='env')
		tmp		<- tolower(do.call('rbind',strsplit(bseq[, SEQ],'')))
		bseq.CP1<- tmp[, seq.int(1, ncol(tmp), by=3)]
		bseq.CP2<- tmp[, seq.int(2, ncol(tmp), by=3)]
		bseq.CP3<- tmp[, seq.int(3, ncol(tmp), by=3)]
		
		tmp		<- bseq.CP3
		tmp2	<- apply( tmp, 1, function(x) paste(x,sep='',collapse=''))	#identical sequences?	
		cat(paste('\nunique sequences, n=',length(unique(tmp2))))		
		tmp2	<- apply( tmp, 2, function(x) paste(x,sep='',collapse=''))	#identical patterns? 
		cat(paste('\nunique patterns, n=',length(unique(tmp2))))
		tmp2	<- apply( tmp, 2, function(x) all(x==x[1]))					#invariant sites?
		length(which(tmp2))
	}
	#	check for drug resistance mutations
	if(0)
	{
		tmp				<- subset(bseq, GENE=='pol')
		bseq.pol.m		<- do.call('rbind',strsplit(tmp[, SEQ],''))
		
		file			<- '/Users/Oliver/git/hivclust/pkg/data/IAS_primarydrugresistance_201303.rda'
		alignment.start	<- 2085
		load(file)
		IAS_primarydrugresistance_201303		<- as.data.table(IAS_primarydrugresistance_201303)
		set(IAS_primarydrugresistance_201303, NULL, "Alignment.nuc.pos", IAS_primarydrugresistance_201303[,Alignment.nuc.pos]-alignment.start+1)
		#pol not in HXB2 consensus coordinates - TODO would have to align against consensus
		z<- seq.rm.drugresistance(bseq.pol.m, IAS_primarydrugresistance_201303, verbose=1, rtn.DNAbin=0)		
	}
	#	concatenate into single DNAbin matrix and save
	tmp		<- dcast.data.table(bseq, TAXON_ID ~ GENE, value.var="SEQ")	
	tmp[, SEQ_ALL:=paste(gag, pol, env, sep='')]
	tmp2	<- tolower(do.call('rbind',strsplit(tmp[, SEQ_ALL],'')))
	rownames(tmp2)	<- tmp[, TAXON_ID]
	seq		<- as.DNAbin(tmp2)		
	outfile	<- paste(infile,'_conc_',insignat,'.R',sep='')
	save(seq, file=paste(outdir, outfile, sep='/'))
}
##--------------------------------------------------------------------------------------------------------
##	check for recombinants
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.checkRecombinants<- function()
{
	require(XML)
	require(ape)
	require(r3SEQ)
	DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140728',sep='/')	
	infile			<- 'ALLv02.n100.rlx.gmrf_conc'
	outfile			<- 'ALLv03.n97.rlx.gmrf'
	insignat		<- 'Sun_Jul_27_09-00-00_2014'
	
	#file		<- paste(indir, '/', infile, '_', insignat, '.R', sep='')
	#load(file)
	#	run 3SEQ
	infile.3seq		<- paste(infile, '_', insignat, '.R', sep='')
	pipeline.recom.run.3seq(indir, infile.3seq, batch.n=1, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
	#	parse 3SEQ output
	argv			<<-	cmd.recombination.process.3SEQ.output(indir, infile.3seq, '', resume=1, verbose=1) 
	argv			<<- unlist(strsplit(argv,' '))
	df.recomb		<- prog.recom.process.3SEQ.output()	
	#	select potential recombinants with p-value < 1e-3
	df.recomb		<- subset( df.recomb, adjp<1e-3 )
	cat(paste('\nfound potential recombinants, n=',nrow(df.recomb)))
	#
	#	remove potential recombinants from XML file
	#
	file			<- '/Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140727/ALLv02.n100.rlx.gmrf_Sun_Jul_27_09-00-00_2014.xml'
	bxml.template	<- xmlTreeParse(file, useInternalNodes=TRUE, addFinalizer = TRUE)
	bseq			<- hivc.beast.get.sequences(bxml.template, verbose=1)	
	bseq			<- merge(bseq, data.table(ALIGNMENT_ID=paste('alignment',1:3,sep=''), GENE=c('env','gag','pol')), by='ALIGNMENT_ID')
	bseq			<- subset(bseq, !TAXON_ID%in%df.recomb[, child])
	#
	#	write XML file with new sequences
	#		
	bxml			<- newXMLDoc(addFinalizer=T)
	bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
	newXMLCommentNode(text=paste("Generated by HIVCLUST from template",file), parent=bxml.beast, doc=bxml, addFinalizer=T)
	#	add new set of ENV sequences into alignment ID 1	
	set( bseq, NULL, 'SEQ', bseq[, tolower(SEQ)] )
	tmp				<- subset(bseq, GENE=='env')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment1", beast.alignment.dataType= "nucleotide", verbose=1)
	#	add new set of GAG sequences into alignment ID 2
	tmp				<- subset(bseq, GENE=='gag')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment2", beast.alignment.dataType= "nucleotide", verbose=1)
	#	add new set of POL sequences into alignment ID 3
	tmp				<- subset(bseq, GENE=='pol')
	tmp2			<- tmp[, do.call('rbind',strsplit(tmp[, SEQ],''))]
	rownames(tmp2)	<- tmp[, TAXON_ID]
	tmp2			<- as.DNAbin(tmp2)
	bxml			<- hivc.beast.add.seq(bxml, tmp2, df=NULL, beast.label.datepos= 2, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment3", beast.alignment.dataType= "nucleotide", verbose=1)
	#	copy rest from template	
	bt.beast		<- getNodeSet(bxml.template, "//beast")[[1]]
	dummy			<- sapply(seq.int( max(which( xmlSApply(bt.beast, xmlName)=="alignment" ))+1, xmlSize(bt.beast) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	change outfile name 
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(outfile,'_',insignat, '.', tail(x,1), sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	#	write to file
	file		<- paste(indir,'/',outfile,'_',insignat,".xml", sep='')
	if(verbose)	cat(paste("\nwrite xml file to",file))
	saveXML(bxml, file=file)

	if(0)
	{
		#	get RAXML trees
		triplets			<- 1
		#triplets			<- 147:nrow(df.recomb)
		dummy	<- lapply(triplets, function(i)
				{				
					i<- 1
					if(verbose)	cat(paste("\nprocess triplet number",i,"\n"))
					argv				<<- cmd.recombination.check.candidates(indir, infile.3seq, '', i, resume=resume, verbose=1, hpc.walltime=1, hpc.q=NA, hpc.mem='500mb', hpc.nproc=1)
					argv				<<- unlist(strsplit(argv,' '))
					prog.recom.get.incongruence()		#this starts ExaML for the ith triplet			
				})	
	}		
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
	
	
	setup.df<- data.table(stat= c('yr.start','yr.end','s.INC.recent','s.INC.recent.len', 's.PREV.min', 's.PREV.max', 'epi.dt'), v=c(1980, 2020, 0.1, 5, 0.01, 0.25, 1/48) )
	setkey(setup.df, stat)	
	
	df.trm	<- as.data.table(read.csv(fin.trm, stringsAsFactors=FALSE, sep=' ', dec='.'))
	setnames(df.trm, c("IdInfector","IdInfected","TimeOfInfection","IsInfectorAcute"), c('IDTR','IDREC','TIME_TR','TR_ACUTE'))		
	#	transmissions happen either at baseline, or at unique times.
	#	the epi simulation allocates transmissions in 1/48 of a year, so draw a uniform number if there are more transmission per TIME_TR
	df.trm	<- df.trm[, {
							z<- TIME_TR
							if(TIME_TR>setup.df['yr.start',][,v] & length(IDTR)>1)
								z<- sort(runif(length(IDTR), z, max(setup.df['yr.end',][,v], z+setup.df['epi.dt',][,v])))
							list(IDTR=IDTR, IDREC=IDREC, TIME_TR.new=z, TR_ACUTE=TR_ACUTE, l=length(IDTR))
						}, by='TIME_TR']
	df.trm[, TIME_TR:=NULL]
	setnames(df.trm, 'TIME_TR.new', 'TIME_TR')
	set(df.trm, NULL, 'YR', df.trm[, floor(TIME_TR)])
	#	check that all transmission times except baseline are unique
	tmp		<- subset(df.trm, TIME_TR>setup.df['yr.start',][,v])
	stopifnot( nrow(tmp)==tmp[,length(unique(TIME_TR))] )
	
	
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

