


prog.hello<- function()	
{
	print('hello')
}
##--------------------------------------------------------------------------------------------------------
##	select between host sequences
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.SSAfg.select.between.host<- function()
{
	label.sep				<- '|'
	label.idx.country.id	<- 2
	label.idx.label			<- 3
	label.idx.ctime			<- 4		
	#read hand aligned sequences
	DATA		<<- '~/git/HPTN071sim/raw_rootseq'
	infile.gag	<- 'PANGEA_SSAfg_gag_140806_n556_final.fasta'
	infile.pol	<- 'PANGEA_SSAfg_pol_140806_n556_final.fasta'
	infile.env	<- 'PANGEA_SSAfg_env_140806_n556_final.fasta'
	seq.gag		<- read.dna( paste(DATA,'/',infile.gag, sep=''), format='fasta' )
	seq.pol		<- read.dna( paste(DATA,'/',infile.pol, sep=''), format='fasta' )
	seq.env		<- read.dna( paste(DATA,'/',infile.env, sep=''), format='fasta' )
	#check seq names are the same
	stopifnot( all(rownames(seq.gag)==rownames(seq.pol)), all(rownames(seq.gag)==rownames(seq.env)) )	
	#ensure unique seq names
	tmp			<- paste( rownames(seq.gag), seq_len(nrow(seq.gag)), sep=label.sep )
	rownames(seq.gag)	<- rownames(seq.pol)	<- rownames(seq.env)	<- tmp	
	tmp			<- strsplit( rownames(seq.gag), label.sep, fixed=1 )
	label.info	<- data.table( SEQ_NAME=sapply(tmp, function(x) paste(x, collapse=label.sep)), COUNTRY_ID= sapply(tmp, '[[', label.idx.country.id), LABEL= sapply(tmp, '[[', label.idx.label), CALENDAR_TIME= sapply(tmp, '[[', label.idx.ctime))
	label.info	<- label.info[, {
							list(SEQ_NAME=SEQ_NAME, COUNTRY_ID=COUNTRY_ID, LABEL_UNIQUE=paste(LABEL,seq_along(SEQ_NAME),sep=label.sep), CALENDAR_TIME=CALENDAR_TIME)
						}, by=LABEL]
	setkey(label.info, SEQ_NAME)
	tmp			<- label.info[rownames(seq.gag), ][, paste('C',COUNTRY_ID,LABEL_UNIQUE,CALENDAR_TIME,sep=label.sep) ]
	rownames(seq.gag)	<- rownames(seq.pol)	<- rownames(seq.env)	<- tmp
	set(label.info, NULL, 'SEQ_NAME', label.info[, paste('C',COUNTRY_ID,LABEL_UNIQUE,CALENDAR_TIME,sep=label.sep) ])	
	#concatenate gag pol env
	seq.c		<- do.call('cbind', list(seq.gag, seq.pol, seq.env))
	#keep only first sequence with same label	
	tmp			<- label.info[, as.integer(sapply( strsplit( LABEL_UNIQUE, label.sep, fixed=1 ), '[[', 2 ))]
	label.info[, LABEL_NO:= tmp]
	label.info[, SELECT:= FALSE]	
	#search for further within host seqs: break up label by '.' or '_'
	setkey(label.info, LABEL)
	tmp			<- label.info[, which( grepl('\\.', LABEL) | grepl('_', LABEL) ) ]	
	#select manually not within host seqs
	include		<- c( "00BW3876_9", "00BW3886_8", "00BW3891_6", "00BW3970_2", "00BW5031_1","02ET_288","4403bmLwk4_fl7",                   
		"702010141_CH141.w12","702010293_CH293.w8a","702010432_CH432.w4","702010440_CH440.w4","703010085_CH085.w4a",
		"703010131_CH131_TF", "703010167_CH167.w8", "703010200_CH200_TFa", "703010228_CH228_TFa", "703010256_CH256.w96", 
		"703010269_CH269.w24", "704010042_CH042.mo6", "705010067_CH067_TF", "705010162_CH162.mo6", "705010185_CH185.mo6", 
		"705010198_CH198_TF", "705010534_CH534.w12", "706010164_CH164.mo6", "707010457_CH457.w8", "89SM_145", "90SE_364", 
		"93MW_965", "96BWMO1_5", "BD16_10", "BD22_11", "BD39_8", "BD9_11", "C.703010159.S.0dps.fl", "C.704010042.S.0dps.fl", 
		"C.705010162.e.wg2", "C.CAP210.w02.0dps.1_00_F4","C.CAP239.w02.0dps.1_02_F32", "C.CAP45.w02.0dps.1_05_T1", "CAP174_4w", "CAP206_8w_F1",	              
		"CAP210_5w", "CAP228_8w_F2", "CAP229_7w", "CAP239_5w_F1", "CAP244_8w_F1", "CAP248_9w", "CAP255_8w_F1", "CAP256_6w",                    
		"CAP257_7w_F1", "CAP30_5w_F4", "CAP45_5w_F1", "CAP61_8w_F3", "CAP63_5w_F4", "CAP65_6w", "CAP84_3w_F2", "CAP85_5w_F1",                  
		"CAP88_5w_F2", "CAP8_3w_F2", "C_ZA_1069MB", "C_ZA_1184MB", "C_ZA_1189MB", "C_ZA_J112MA", "TV001_patent", "TV002_patent",
		"ZM246F_flA1", "ZM247F_flA1", "ZM249M_flC1", "pZAC_R3714")
	include		<- c(include, label.info[ !grepl('\\.', LABEL) & !grepl('_', LABEL),   ][, LABEL])
	#these are excluded based on the '.' and '_' search:
	#c("4403bmLwk4_fl11","702010293_CH293.w8b","703010200_CH200_TFb","703010200_CH200_TFc", "703010228_CH228_TFb",
	#	"703010256_CH256_TF", "704010042_CH042_TF", "705010162_CH162_TF", "705010185_CH185_TF", "706010164_CH164_TF",
	#	"C.CAP210.w02.0dps.1_00_T11", "C.CAP210.w02.0dps.1_00_T2B", "C.CAP210.w02.0dps.1_00_T3", "C.CAP210.w02.0dps.1_00_T36",
	#	"C.CAP210.w02.0dps.1_00_T3C", "C.CAP210.w02.0dps.1_00_T4", "C.CAP210.w02.0dps.1_00_T43", "C.CAP210.w02.0dps.1_00_T5",
	#	"C.CAP210.w02.0dps.1_00_T6", "C.CAP210.w05.21dps.2_00", "C.CAP210.w12.70dps.2_05_T13", "C.CAP210.w12.70dps.2_05_T13C",
	#	"C.CAP210.w12.70dps.2_05_T2", "C.CAP210.w12.70dps.2_05_T39w", "C.CAP210.w12.70dps.2_05_T42", "C.CAP210.w12.70dps.2_05_T5",
	#	"C.CAP210.w12.70dps.2_05_T8", "C.CAP210.w26.168dps.3_10_T13B", "C.CAP210.w26.168dps.3_10_T20w", "C.CAP210.w26.168dps.3_10_T23B", 
	#	"C.CAP210.w26.168dps.3_10_T24B", "C.CAP210.w26.168dps.3_10_T24C", "C.CAP210.w26.168dps.3_10_T28", "C.CAP210.w26.168dps.3_10_T40",
	#	"C.CAP210.w26.168dps.3_10_T42", "C.CAP210.w26.168dps.3_10_T43B", "C.CAP210.w26.168dps.3_10_T47", "C.CAP210.w26.168dps.3_10_T49B",
	#	"C.CAP239.w02.0dps.1_02_T8", "C.CAP239.w05.21dps.2_00", "C.CAP239.w05.21dps.2_00_T11", "C.CAP239.w05.21dps.2_00_T17", 
	#	"C.CAP239.w05.21dps.2_00_T18", "C.CAP239.w05.21dps.2_00_T19", "C.CAP239.w05.21dps.2_00_T21", "C.CAP239.w05.21dps.2_00_T3",
	#	"C.CAP239.w05.21dps.2_00_T49", "C.CAP239.w05.21dps.2_00_T8", "C.CAP239.w11.63dps.2_05_T37", "C.CAP239.w11.63dps.2_05_T47",
	#	"C.CAP239.w117.805dps.4_21_T44", "C.CAP239.w117.805dps.4_21_T46", "C.CAP239.w117.805dps.4_21_T50", "C.CAP239.w22.140dps.3_09_F1", 
	#	"C.CAP239.w22.140dps.3_09_T10", "C.CAP239.w22.140dps.3_09_T17", "C.CAP239.w22.140dps.3_09_T20", "C.CAP239.w22.140dps.3_09_T36", 
	#	"C.CAP239.w22.140dps.3_09_T39", "C.CAP239.w22.140dps.3_09_W16", "C.CAP45.w02.0dps.1_05_F3", "C.CAP45.w02.0dps.1_05_T2", "C.CAP45.w05.21dps.2_00", 
	#	"C.CAP45.w05.21dps.2_00_T11", "C.CAP45.w05.21dps.2_00_T12", "C.CAP45.w05.21dps.2_00_T14", "C.CAP45.w05.21dps.2_00_T5", "C.CAP45.w05.21dps.2_00_T9", 
	#	"C.CAP45.w12.70dps.2_05_F1", "C.CAP45.w12.70dps.2_05_T11", "C.CAP45.w12.70dps.2_05_T13b", "C.CAP45.w12.70dps.2_05_T18",
	#	"C.CAP45.w65.455dps.4_17_T14", "C.CAP45.w65.455dps.4_17_T14B", "ZM246F_flA10", "ZM246F_flA2", "ZM246F_flA6", "ZM246F_flB1", 
	#	"ZM246F_flC12", "ZM246F_flC3" , "ZM246F_flC5", "ZM246F_flC7", "ZM246F_flD5", "ZM247F_flA12", "ZM247F_flA2", "ZM247F_flB8", 
	#	"ZM247F_flB9", "ZM247F_flE10", "ZM247F_flE11", "ZM247F_flE3", "ZM247F_flF10", "ZM247F_flF7", "ZM247F_flG11", "ZM247F_flH1",
	#	"ZM249M_flC5", "ZM249M_flE10", "ZM249M_flE8", "ZM249M_flF1", "TV001_patent", "TV001_patent", "TV002_patent", "chimeric_MJ4")
	exclude	<- c(	'C|ZA|C_ZA_1184MB|1|2000','C|ZA|C_ZA_1189MB|1|2000','C|ZA|C_ZA_J112MA|1|2000','C|ZA|C_ZA_1069MB|1|2000',
					'C|ZA|TV001_patent|1|1998','C|ZA|TV002_patent|1|1998','C|ZA|03ZASK212B1|1|2003',
					'C|ZA|C.CAP239.w02.0dps.1_02_F32|1|2005',"C|ZA|C.CAP210.w02.0dps.1_00_F4|1|2005","C|ZA|C.CAP45.w02.0dps.1_05_T1|1|2005",
					'C|MW|C.703010159.S.0dps.fl|1|2007','C|ZA|C.704010042.S.0dps.fl|1|2007','C|ZA|705010162_CH162.mo6|1|2007',
					'C|BW|96BW15C05|1|1996','C|BW|96BW15B03|1|1996',
					'C|BW|96BW01B22|1|1996','C|BW|96BW01B03|1|1996',"C|BW|96BW11B01|1|1996","C|BW|96BW1104|1|1996",
					'C|BW|96BW16B01|1|1996','C|BW|96BW1626|1|1996', "C|BW|96BW0502|1|1996", "C|BW|96BW06H51|1|1996","C|BW|96BW06|1|1996",
					"C|BW|96BW0407|1|1996","C|BW|96BW0402|1|1996","C|BW|96BW0410|1|1996","C|BW|96BW0408|1|1996")
	set(label.info, label.info[, which(  LABEL_NO==1L & LABEL%in%include )], 'SELECT', TRUE)
	set(label.info, label.info[, which(  SEQ_NAME%in%exclude )], 'SELECT', FALSE)
	#select
	tmp			<- subset(label.info, SELECT)[, SEQ_NAME]
	seq.c		<- seq.c[tmp, ]
	#save the whole lot
	tmp			<- rownames(seq.c)
	seq.gag		<- seq.gag[tmp, ]
	seq.pol		<- seq.pol[tmp, ]
	seq.env		<- seq.env[tmp, ]
	seq			<- seq.c			#need 'seq' because expected for 3SEQ
	outdir		<- '~/duke/2014_Gates/methods_comparison_rootseqsim/140811'
	outfile		<- 'PANGEA_SSAfgBwh_140811_n415_final.R'
	file		<- paste(outdir, '/', outfile, sep='')
	save(seq, seq.gag, seq.pol, seq.env, file=file)
	#check for recombinants
}
##--------------------------------------------------------------------------------------------------------
##	run 3SEQ
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.3SEQ.SSAfg.rm.recombinants<- function()
{
	require(XML)
	require(ape)
	require(r3SEQ)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140811',sep='/')	
	infile			<- 'PANGEA_SSAfgBwh_140811_n415_final.R'
	outfile			<- 'PANGEA_SSAfgBwhRc-_140811_n390.R'
	#
	#	run 3SEQ
	r3seq.pipe.run.3seq(indir, infile, batch.n=5, hpc.walltime=1, hpc.q=NA, hpc.mem="500mb", hpc.nproc=1)
	#	parse 3SEQ output
	argv			<<-	r3seq.cmd.process.3SEQ.output(indir, infile, '', resume=1, verbose=1) 
	argv			<<- unlist(strsplit(argv,' '))
	df.recomb		<- r3seq.prog.process.3SEQ.output()	
	#	subset( df.recomb, adjp<1e-4 & min_rec_length>500)[, hist(log10(adjp), breaks=100)]
	df.recomb		<- subset( df.recomb, adjp<1e-7 & min_rec_length>500)
	cat(paste('\nfound potential recombinants, n=',nrow(df.recomb)))
	#
	file		<- paste(indir, '/', infile, sep='')
	load(file)
	tmp			<- setdiff(rownames(seq), df.recomb[, child])
	seq.gag		<- seq.gag[tmp,]
	seq.pol		<- seq.pol[tmp,]
	seq.env		<- seq.env[tmp,]
	seq			<- seq[tmp,]
	file		<- paste(indir, '/', outfile, sep='')
	save(seq, seq.gag, seq.pol, seq.env, file=file)
}
##--------------------------------------------------------------------------------------------------------
##	run BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.run<- function()
{
	#DATA		<<- "/work/or105/Gates_2014"
	DATA		<<- '/Users/Oliver/duke/2014_Gates'	
	indir		<- paste(DATA,'methods_comparison_rootseqsim/140811',sep='/')
	#search for XML files in indir
	infiles		<- list.files(indir, pattern=paste(".xml$",sep=''))
	insignat	<- ''	
	hpc.ncpu	<- 8
	
	for(infile in infiles)
	{
		infile		<- substr(infile, 1, nchar(infile)-4) 		
		cmd			<- hivc.cmd.beast.runxml(indir, infile, insignat, prog.beast=PR.BEAST, prog.beast.opt=" -beagle -working", hpc.tmpdir.prefix="beast", hpc.ncpu=hpc.ncpu)
		tmp			<- paste(infile,'.timetrees',sep='')	
		cmd			<- paste(cmd, hivc.cmd.beast.read.nexus(indir, tmp, indir, tree.id=NA, method.node.stat='any.node'), sep='\n')
		cmd			<- paste(cmd, hivc.cmd.beast.run.treeannotator(indir, infile, insignat, prog.beastmcc=PR.BEASTMCC, beastmcc.burnin=500, beastmcc.heights="median"), sep='\n')
		cat(cmd)	
		cmd			<- hivc.cmd.hpcwrapper(cmd, hpc.q="pqeph", hpc.nproc=hpc.ncpu, hpc.walltime=91, hpc.mem="3700mb")		
		outdir		<- indir
		outfile		<- paste("b2m.",paste(strsplit(date(),split=' ')[[1]],collapse='_',sep=''),sep='')					
		hivc.cmd.hpccaller(outdir, outfile, cmd)		
	}
}
##--------------------------------------------------------------------------------------------------------
##	create BEAST XML file
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.createXML<- function()
{
	require(hivclust)
	require(XML)
	require(ape)
	require(r3SEQ)
	#DATA			<<- "/work/or105/Gates_2014"
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	
	if(1)	
	{		
		#
		#	to define sequences for each BEAST run
		#	compute NJ tree and define clusters
		#
		infile.beast	<- '/Users/Oliver/git/HPTN071sim/raw_rootseq/BEAST_template_v08.xml'
		indir			<- paste(DATA,'methods_comparison_rootseqsim/140813',sep='/')	
		infile			<- 'PANGEA_SSAfgBwhRc-_140811_n390.R'
		file			<- paste(indir, '/', infile, sep='')
		load(file)		
		#	remove sequences without calendar time 		
		label.sep				<- '|'
		label.idx.ctime			<- 5		
		tmp						<- sapply( strsplit( rownames(seq.gag), label.sep, fixed=1 ), '[[', label.idx.ctime )
		tmp						<- rownames(seq.gag)[ which(is.na(as.numeric(tmp))) ]
		cat(paste('\nExclude sequences with no calendar date, ', paste(tmp, collapse=' ')))
		tmp						<- setdiff(rownames(seq.gag), tmp)		
		seq.gag					<- seq.gag[tmp,]
		seq.pol					<- seq.pol[tmp,]
		seq.env					<- seq.env[tmp,]
		seq						<- seq[tmp,]
		#	get NJ tree and plot
		tmp				<- dist.dna( seq )
		seq.ph			<- nj(tmp)				
		file			<- paste( indir, '/', substr(infile,1,nchar(infile)-2), '_njtree.pdf', sep='' )	
		pdf(file=file, w=10, h=80)
		plot(seq.ph, show.tip=TRUE)
		dev.off()			
		#
		#	get 3 sequence pools of equal size
		#
		pool.n			<- 3
		tmp				<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)
		thresh.brl		<- 0.055
		clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=tmp, retval="all")
		#	allocate clustering tips into 3 distinct clusters
		seq.clumem		<- data.table( PH_NODE_ID=seq_len(Ntip(seq.ph)), CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph)) ] )
		setkey(seq.clumem, CLU_ID)		
		tmp				<- which(!is.na(seq.clumem[, CLU_ID]))
		tmp				<- seq.clumem[tmp,][, list(CLU_N=-length(PH_NODE_ID)), by='CLU_ID']
		setkey(tmp, CLU_N)
		set(tmp, NULL, 'POOL_ID', tmp[, cumsum(-CLU_N)]) 		
		set(tmp, NULL, 'POOL_ID', tmp[, ceiling( POOL_ID / max(POOL_ID) * pool.n ) ] )
		seq.clumem		<- merge(seq.clumem, subset(tmp, select=c(CLU_ID, POOL_ID)), by='CLU_ID', all.x=TRUE)
		#	allocate non-clustering tips into 3 distinct clusters
		tmp				<- subset(seq.clumem,!is.na(POOL_ID))[, list(NOCLU_N= ceiling( nrow(seq.clumem) / pool.n ) - length(PH_NODE_ID)), by='POOL_ID']		
		set(tmp, 1L, 'NOCLU_N', tmp[1,NOCLU_N] - ( tmp[, sum(NOCLU_N)] - ( Ntip(seq.ph) - nrow(subset(seq.clumem,!is.na(POOL_ID))) )) )		
		set(seq.clumem, seq.clumem[, which(is.na(POOL_ID))], 'POOL_ID',  rep(tmp[,POOL_ID], tmp[,NOCLU_N]) )
		seq.clumem[, table(POOL_ID)]	
		#
		#	for each sequence pool, set up BEAST run
		#
		verbose			<- 1
		bxml.template	<- xmlTreeParse(infile.beast, useInternalNodes=TRUE, addFinalizer = TRUE)
		for(pool.id in seq_len(pool.n))
		{
			pool.infile		<- paste(  substr(infile,1,nchar(infile)-2),'_pool',pool.id, sep='' )
			pool.seqnames	<- seq.ph$tip.label[ subset(seq.clumem, POOL_ID==pool.id)[, PH_NODE_ID] ]
			cat(paste('\ncreate BEAST XML file for seqs=',paste(pool.seqnames, collapse=' ')))
			#	write XML file with new sequences
			bxml			<- newXMLDoc(addFinalizer=T)
			bxml.beast		<- newXMLNode("beast", doc=bxml, addFinalizer=T)
			newXMLCommentNode(text=paste("Generated by HIVCLUST from template",infile.beast), parent=bxml.beast, doc=bxml, addFinalizer=T)
			#	add new set of ENV sequences into alignment ID 1
			tmp				<- seq.env[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment1", beast.alignment.dataType= "nucleotide", verbose=1)
			#	add new set of GAG sequences into alignment ID 2
			tmp				<- seq.gag[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment2", beast.alignment.dataType= "nucleotide", verbose=1)
			#	add new set of POL sequences into alignment ID 3
			tmp				<- seq.pol[pool.seqnames,]
			bxml			<- hivc.beast.add.seq(bxml, tmp, df=NULL, beast.label.datepos= 5, beast.label.sep= '|', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment3", beast.alignment.dataType= "nucleotide", verbose=1)
			#	copy from template	
			bt.beast		<- getNodeSet(bxml.template, "//beast")[[1]]
			dummy			<- sapply(seq.int( 1, xmlSize(bt.beast) ), function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					})
			#	change gmrf dimensions	
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1  
			tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
			if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
			tmp			<- tmp[[1]]
			xmlAttrs(tmp)["dimension"]	<-	length(pool.seqnames)-1			
			#	change outfile name 
			bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
			tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
			tmp			<- gsub("(time).","time",tmp,fixed=1)
			tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
			tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(pool.infile, '.', tail(x,1), sep=''))
			dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
			#	write to file
			file		<- paste(indir,'/',pool.infile,".xml", sep='')
			if(verbose)	cat(paste("\nwrite xml file to",file))
			saveXML(bxml, file=file)
		}		
	}
	#
}
##--------------------------------------------------------------------------------------------------------
##	get anecestral sequences from BEAST XML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.SSAfg.getancestralseq.from.output<- function()
{
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	outdir				<- indir
	#	search for BEAST output
	files				<- list.files(indir)
	files				<- files[ sapply(files, function(x) grepl('pool[0-9].R$',x) ) ]	
	if(!length(files))	stop('cannot find files matching criteria')
	
	#	load and process BEAST PARSER output
	anc.seq				<- lapply(files, function(file)
			{
				cat(paste('\nProcess file=', file  ))
				load( paste(indir, file, sep='/') )	#	expect tree, node.stat
				#	compute gag pol env ancestral sequences		
				anc.seq	<- PANGEA.RootSeqSim.get.ancestral.seq(tree, node.stat, tree.id.sep='_', tree.id.idx.mcmcit=2, tree.id.burnin=tree.id.burnin, label.sep=tree.id.labelsep, label.idx.ctime=5)				
				set(anc.seq, NULL, 'LABEL', anc.seq[, paste( substr(file,1,nchar(file)-2), LABEL, sep=tree.id.labelsep )] )				
				set(anc.seq, NULL, 'TREE_ID', NULL )
				set(anc.seq, NULL, 'NODE_ID', NULL )
				set(anc.seq, NULL, 'BEAST_MCMC_IT', NULL )
				anc.seq
			})
	anc.seq				<- do.call('rbind',anc.seq)
	#
	#	return DNAbin
	#
	anc.seq.gag				<- tolower(do.call('rbind',strsplit(anc.seq[, GAG],'')))
	rownames(anc.seq.gag)	<- anc.seq[, LABEL]
	anc.seq.gag				<- as.DNAbin(anc.seq.gag)		
	anc.seq.pol				<- tolower(do.call('rbind',strsplit(anc.seq[, POL],'')))
	rownames(anc.seq.pol)	<- anc.seq[, LABEL]
	anc.seq.pol				<- as.DNAbin(anc.seq.pol)		
	anc.seq.env				<- tolower(do.call('rbind',strsplit(anc.seq[, ENV],'')))
	rownames(anc.seq.env)	<- anc.seq[, LABEL]
	anc.seq.env				<- as.DNAbin(anc.seq.env)	
	
	set( anc.seq, NULL, 'GAG', NULL )
	set( anc.seq, NULL, 'POL', NULL )
	set( anc.seq, NULL, 'ENV', NULL )
	anc.seq.info			<- anc.seq
	#anc.seq					<- cbind(anc.seq.gag, anc.seq.pol, anc.seq.env)
	#
	outfile				<- paste( substr(files[1],1,nchar(files[1])-7), 'AncSeq.R',sep='' )
	file				<- paste(outdir, outfile, sep='/')
	cat(paste('\nwrite Ancestral Sequences to ',file))
	save(file=file, anc.seq.gag, anc.seq.pol, anc.seq.env, anc.seq.info)
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, create random draw to check
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.createdataset<- function()
{
	tree.id.burnin		<- 2e7
	tree.id.labelsep	<- '|'
	dir.name			<- '/Users/Oliver/duke/2014_Gates'  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')	
	infile				<- "PANGEA_SSAfgBwhRc-_140811_n390_AncSeq.R"	
	outdir				<- indir
	outfile.fg			<- infile
	outfile.partialpol	<- "PANGEA_SSApolBwhRc-_140811_n390_AncSeq.R"
	outsignat			<- "Mon_Aug_17_17:05:23_2014"	
	#	load ancestral sequences
	load( paste(indir, infile, sep='/') ) #	expect anc.seq.gag, anc.seq.pol, anc.seq.env, anc.seq.info
	#
	#	basic checks
	#	
	if(0)
	{
		tmp					<- anc.seq.info[, 	sapply( strsplit(LABEL, tree.id.labelsep, fixed=TRUE),'[[', 1) ]
		anc.seq.info[, POOL:= regmatches(tmp, regexpr('pool[0-9]+', tmp))]
		set(anc.seq.info, NULL, 'POOL', anc.seq.info[, as.numeric(substr(POOL, 5, nchar(POOL)))])
		ggplot(anc.seq.info, aes(x=CALENDAR_TIME)) + geom_histogram(binwidth=2) + facet_grid(.~POOL, margins=0)
		#	most sequences between 1940 - 1980
		subset(anc.seq.info, floor(CALENDAR_TIME)==1940)		
	}
	#	sample data set to run ExaML
	anc.seq.info[, CALENDAR_YR:=anc.seq.info[, floor(CALENDAR_TIME)]]	
	#	~ 1400 sequences from 1940
	#	sample 10 times 1e3 sequences randomly from exp increasing prevalence and estimate tree to calculate root to tip divergence + clustering on fg and partial pol
	s.seed						<- 42
	s.PREV.MAX					<- 0.25
	s.PREV.MIN					<- 0.01
	s.RANGE						<- 40
	s.size						<- 1e3
	s.baseline.ancseq.time		<- 1940
	s.baseline.calendar.time	<- 1980
	s.LENGTH.PARTIAL.POL		<- 1500
	tree.id.labelidx.ctime		<- 4		
	tmp							<- log( 1+s.PREV.MAX-s.PREV.MIN ) / (s.RANGE-1)
	tmp							<- exp( tmp*seq.int(0,s.RANGE-1) ) - 1 + s.PREV.MIN
	seq.s						<- c(s.PREV.MIN, diff(tmp))
	seq.s						<- data.table( SEQ_N= round( seq.s/s.PREV.MAX*s.size ), ANCSEQ_YR=seq_along(seq.s)+s.baseline.ancseq.time-1 )
	anc.seq.info				<- subset(anc.seq.info,  CALENDAR_YR>=s.baseline.ancseq.time & CALENDAR_YR<=(s.baseline.ancseq.time+s.RANGE-1))
	
	#	draw partial pol genome sequences - length is first 1500 sites
	#	and draw full genome sequences
	set.seed(s.seed)	
	for(check.draw in 1:5)
	{
		#	draw a large enough number of ancseq labels from the 3 pools for each year to accommodate SEQ_N/3 anc seqs from each pool
		anc.seq.infodraw			<- anc.seq.info[, {
					tmp	<- seq.s$SEQ_N[ which(seq.s$ANCSEQ_YR==CALENDAR_YR) ]
					list(LABEL= sample(LABEL, 2*tmp, replace=FALSE), SEQ_N=tmp, SEQ_N_GRACE=2*tmp)
				}, by='CALENDAR_YR']	
		anc.seq.draw				<- anc.seq.pol[anc.seq.infodraw[, LABEL], seq_len(s.LENGTH.PARTIAL.POL)]
		#	make sure the drawn sequences are unique on partial POL
		anc.seq.draw				<- seq.unique(anc.seq.draw)
		anc.seq.infodraw			<- merge( data.table(LABEL=rownames(anc.seq.draw)), anc.seq.infodraw, by='LABEL' )	
		stopifnot( anc.seq.infodraw[, list(SEQ_N_GRACE=length(LABEL), SEQ_N=SEQ_N[1]), by='CALENDAR_YR'][, all(SEQ_N_GRACE>=SEQ_N)] )
		anc.seq.infodraw			<- anc.seq.infodraw[, list(LABEL= LABEL[seq_len(SEQ_N[1])]), by='CALENDAR_YR']
		#	set new calendar time for sequences
		set(anc.seq.infodraw, NULL, 'LABEL_NEW', anc.seq.infodraw[, as.numeric( sapply( strsplit(LABEL,tree.id.labelsep,fixed=TRUE), '[[', tree.id.labelidx.ctime) ) ])
		set(anc.seq.infodraw, NULL, 'LABEL_NEW', anc.seq.infodraw[, LABEL_NEW-s.baseline.ancseq.time+s.baseline.calendar.time])	
		anc.seq.infodraw			<- anc.seq.infodraw[,	{
					tmp							<- strsplit(LABEL,tree.id.labelsep,fixed=TRUE)[[1]]
					tmp[tree.id.labelidx.ctime]	<- LABEL_NEW
					list(LABEL_NEW=paste(tmp, collapse=tree.id.labelsep,sep=''))
				}, by='LABEL']
		setkey(anc.seq.infodraw, LABEL)
		#	select partial POL seqs and save to file
		anc.seq.draw				<- anc.seq.pol[anc.seq.infodraw[, LABEL], seq_len(s.LENGTH.PARTIAL.POL)]		
		rownames(anc.seq.draw)		<- anc.seq.infodraw[ rownames(anc.seq.draw), ][, LABEL_NEW]
		file						<- paste( outdir, '/', substr(outfile.partialpol,1,nchar(outfile.partialpol)-2),'_checkdraw', check.draw,'_', insignat, '.R', sep='' )
		cat(paste('\nsave to file', file))
		save(anc.seq.draw, file=file)
		#	select the same full genome seqs and save to file
		anc.seq.draw				<- do.call( 'cbind', list( anc.seq.gag[anc.seq.infodraw[, LABEL], ], anc.seq.pol[anc.seq.infodraw[, LABEL], ], anc.seq.env[anc.seq.infodraw[, LABEL], ] ) ) 
		rownames(anc.seq.draw)		<- anc.seq.infodraw[ rownames(anc.seq.draw), ][, LABEL_NEW]
		file						<- paste( outdir, '/', substr(outfile.fg,1,nchar(outfile.fg)-2),'_checkdraw', check.draw,'_', insignat, '.R', sep='' )
		cat(paste('\nsave to file', file))
		save(anc.seq.draw, file=file)
	}
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, run ExaML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.runExaML<- function()
{
	#DATA			<<- "/work/or105/Gates_2014"
	DATA				<<- '/Users/Oliver/duke/2014_Gates'
	dir.name			<- DATA  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	#	ExaML bootstrap args
	bs.from		<- 0
	bs.to		<- 1
	bs.n		<- 100
	
	#	search for 'checkdraw' files
	infiles		<- list.files(indir)
	infiles		<- infiles[ sapply(infiles, function(x) grepl('.*checkdraw[0-9]+.*R$',x) ) ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	
	outdir		<- indir
	for(infile in infiles)
	{
		#infile		<- files[1]
		infile		<- substr(infile, 1, nchar(infile)-2)
		insignat	<- regmatches(infile, regexpr('checkdraw[0-9]+_.*', infile))
		insignat	<- regmatches(insignat,regexpr('_.*',insignat))
		insignat	<- substr(insignat,2,nchar(insignat))
		infile		<- regmatches(infile, regexpr('.*checkdraw[0-9]+', infile))
		
		
		cmd			<- hivc.cmd.examl.bootstrap(indir, infile, insignat, insignat, bs.from=bs.from, bs.to=bs.to,bs.n=bs.n,outdir=outdir, resume=1, verbose=1)
		dummy		<- lapply(cmd, function(x)
				{				
					x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q= NA, hpc.mem="450mb", hpc.nproc=1)
					#x		<- hivc.cmd.hpcwrapper(x, hpc.walltime=24, hpc.q="pqeph", hpc.mem="3850mb", hpc.nproc=8)
					signat	<- paste(strsplit(date(),split=' ')[[1]],collapse='_',sep='')
					outfile	<- paste("exa",signat,sep='.')
					#cat(x)
					hivc.cmd.hpccaller(outdir, outfile, x)
					Sys.sleep(1)
				})
		stop()
	}
}
##--------------------------------------------------------------------------------------------------------
##	check ancestral sequences from BEAST XML, run ExaML
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.SIMU.SSAfg.checkancestralseq.evalExaML<- function()
{
	#DATA			<<- "/work/or105/Gates_2014"
	DATA				<<- '/Users/Oliver/duke/2014_Gates'
	dir.name			<- DATA  	
	indir				<- paste(dir.name,'methods_comparison_rootseqsim/140813',sep='/')
	label.sep			<- '|'
	label.idx.ctime		<- 4
	clu.thresh.bs		<- 0.9
	clu.thresh.brl		<- 0.045	#0.035
	bs.n				<- 100
	#	search for 'checkdraw examl' files
	infiles				<- list.files(indir)
	infiles				<- infiles[ grepl('.*checkdraw[0-9]+_examl.*newick$',infiles)  ]	
	if(!length(infiles))	stop('cannot find files matching criteria')
	
	#	read tree
	infile	<- infiles[1]
	file	<- paste(indir, infile, sep='/')
	ph		<- read.tree(file)
	ph		<- ladderize(ph)	
	#file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Tree.pdf', sep='' )
	#pdf(file=file, h=150, w=10)
	#plot(ph, cex=0.7)
	#dev.off()
	#
	#	check root to tip divergence
	#
	tmp		<- node.depth.edgelength(ph)
	ph.info	<- data.table(LABEL=ph$tip.label, ROOT2TIP=tmp[seq_len(Ntip(ph))] )
	set(ph.info, NULL, 'CALENDAR_TIME', ph.info[, as.numeric(sapply(strsplit(LABEL, label.sep, fixed=TRUE),'[[',label.idx.ctime))] )
	tmp		<- lm(ROOT2TIP~CALENDAR_TIME, data=ph.info)		 
	set( ph.info, NULL, 'ROOT2TIP_LM', predict(tmp, type='response') ) 	
	tmp2	<- c( R2=round(summary(tmp)$r.squared,d=3), SLOPE= as.numeric(round(coef(tmp)['CALENDAR_TIME'],d=4)), TMRCA=as.numeric(round( -coef(tmp)['(Intercept)']/coef(tmp)['CALENDAR_TIME'], d=1 )) )
	ggplot(ph.info, aes(x=CALENDAR_TIME, y=ROOT2TIP)) + geom_point(alpha=0.5) + geom_line(aes(y=ROOT2TIP_LM)) +
			#scale_x_continuous(breaks=seq(1980,2020,2)) +						
			labs(x='Sequence sampling date', y='root-to-tip divergence') +
			annotate("text", x=ph.info[, min(CALENDAR_TIME)], y=ph.info[, 0.9*max(ROOT2TIP)], label=paste("R2=", tmp2['R2'],'\nSlope=',tmp2['SLOPE'],'\nTMRCA=',tmp2['TMRCA'], sep=''), hjust = 0, size = 4) +
			theme(legend.position=c(0,1), legend.justification=c(0,1))
	file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Root2Tip.pdf', sep='' )
	ggsave(file=file, w=10, h=6)
	#	check brl divergence
	#brl.tips	<- distTips(ph , method='patristic')
	#
	#	check clustering among root seqs - this should be minimal at cut-off 0.045
	#
	dist.brl		<- hivc.clu.brdist.stats(ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade )
	quantile(dist.brl,p=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5))
	hist(dist.brl, breaks=30)
	ph.node.bs		<- as.numeric( ph$node.label )		
	ph.node.bs[is.na(ph.node.bs)]	<- 0
	ph.node.bs		<- ph.node.bs/bs.n
	ph$node.label	<- ph.node.bs	
	ph.clu			<- hivc.clu.clusterbythresh(ph, thresh.nodesupport=clu.thresh.bs, thresh.brl=clu.thresh.brl, dist.brl=dist.brl, nodesupport=ph.node.bs, retval="all")
	#	plot clustering tree
	file	<- paste( indir, '/', substr(infile, 1, nchar(infile)-7), '_Clutree.pdf', sep='' )
	pdf(file=file, h=150, w=10)
	hivc.clu.plot(ph, ph.clu[["clu.mem"]], show.tip.label=TRUE, cex.edge.incluster=1.5)
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	get anecestral sequences from BEAST XML - developper version
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.BEAST.devel.getancestralseq.from.output<- function()
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
##	prepare LANL download
##	downloaded 556 full genome sequences from SSA on 06-08-2014
##--------------------------------------------------------------------------------------------------------
project.PANGEA.RootSeqSim.DATA.SSAfg.process.LosAlamos<- function()
{
	DATA			<<- '/Users/Oliver/duke/2014_Gates'
	indir			<- paste(DATA,'methods_comparison_rootseqsim/140806',sep='/')
	infile.fasta	<- 'PANGEA_SSAfg_140806_n557.fasta'
	infile.fasta.gag<- 'PANGEA_SSAfg_gag_140806_n557.fasta'
	infile.fasta.pol<- 'PANGEA_SSAfg_pol_140806_n557.fasta'
	infile.fasta.env<- 'PANGEA_SSAfg_env_140806_n557.fasta'
	file	<- paste(indir, '/', infile.fasta,sep='')
	seq		<- read.dna(file, format='fasta')	
	
	label.sep				<- '|'
	label.idx.country.id	<- 2
	label.idx.label			<- 3
	label.idx.ctime			<- 4
	tmp						<- strsplit(names(seq), label.sep, fixed=1)
	
	seq.label				<- data.table(	LABEL= names(seq),
			COUNTRY_ID= sapply(tmp,'[[',label.idx.country.id),
			NAME= sapply(strsplit(names(seq), label.sep, fixed=1),'[[',label.idx.label),
			CALENDAR_TIME= sapply(strsplit(names(seq), label.sep, fixed=1),'[[',label.idx.ctime))
	#	remove sequences without a sampling date								
	seq.label	<- subset( seq.label, !is.na(as.numeric(CALENDAR_TIME)) )							
	set(seq.label, NULL, 'CALENDAR_TIME', seq.label[, as.numeric(CALENDAR_TIME)])
	#	histogram of sample dates and sample location	
	ggplot(seq.label, aes(x=CALENDAR_TIME)) + geom_histogram(binwidth=1)
	ggplot(seq.label, aes(x=COUNTRY_ID)) + geom_histogram()
	#
	#	some seqs are crap and have duplicate runs, need to edit manually. ISOLATING GAG POL ENV
	#
	data(refseq_hiv1_hxb2)
	seq.gag.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_gag_DNA.fasta', format='fasta'))
	seq.pol.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_pol_DNA.fasta', format='fasta'))
	seq.env.c	<- as.character(read.dna(file='~/git/hivclust/pkg/data/LosAlamos_HIV1B_CONSENSUS_2004_env_DNA.fasta', format='fasta'))
	#	located gag start at 1009 manually, rm everything before
	seq			<- as.character(seq)
	for(i in seq_along(seq))	seq[[i]][1:1008]<- "-"
	tmp		<- c(1041:1043, 1218:1220, 1321:1322, 1326, 1341:1352, 1366:1389, 1395:1397, 1402:1407, 1411:1413, 1417:1440, 1443:1445, 1447:1449, 1896:1898, 2271:2291, 2455:2466, 2507:2509, 
			2511:2513, 2521:2523, 2551:2577, 2596, 2640, 2644:2655, 2668:2670 )	
	for(i in seq_along(seq))	seq[[i]][tmp]	<- '-'	
	#	add '-' to get all seqs of same length					
	tmp			<- max(sapply(seq, length))
	seq			<- lapply(seq, function(x) c(x, rep('-',tmp-length(x))) 	)
	#	get into matrix DNAbin	
	seq			<- as.DNAbin(do.call('rbind', seq))
	#	rm '-' columns so far
	seq			<- seq.rmallchar(seq, rm.char='-', verbose=1)
	file		<- paste(indir, '/fixup1_', infile.fasta,sep='')
	write.dna(seq, file=file, format='fasta')
	#	located pol/prot start at 1657 manually, cut and rm gaps and re-align gag
	seq				<- as.character(seq)
	seq.gag			<- seq[, 1:1656]
	tmp				<- rownames(seq.gag)
	seq.gag			<- lapply(seq_len(nrow(seq.gag)), function(i){	seq.gag[i, seq.gag[i,]!="-" ]	})
	seq.gag[[length(seq.gag)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.gag[[length(seq.gag)+1]]	<- seq.gag.c['CONSENSUS_C',]
	names(seq.gag)	<- c(tmp,'HXB2','CONSENSUS_C')	
	seq.gag			<- as.DNAbin(seq.gag)
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')
	#	align GAG
	tmp			<- hivc.cmd.clustalo(indir, infile.fasta.gag, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#	continue fixup for POL	
	file			<- paste(indir, '/fixup1_', infile.fasta,sep='')
	seq				<- read.dna(file, format='fasta')
	seq				<- as.character(seq)
	seq[, 1:1656]	<- "-"
	#	located end of POL at 4597	- found length issues, align with own HXB2 
	seq.pol			<- seq[, 1:4700]	
	tmp				<- rownames(seq.pol)
	seq.pol			<- lapply(seq_len(nrow(seq.pol)), function(i){	seq.pol[i, seq.pol[i,]!="-" ]	})	
	seq.pol[[length(seq.pol)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.pol[[length(seq.pol)+1]]	<- seq.pol.c['CONSENSUS_C',]	
	names(seq.pol)	<- c(tmp,'HXB2','CONSENSUS_C')	
	seq.pol			<- as.DNAbin( seq.pol )
	file			<- paste(indir, '/', infile.fasta.pol,sep='')
	write.dna(seq.pol, file=file, format='fasta')
	#	align POL
	#	/Users/Oliver/git/hivclust/pkg/inst/mafft-mac/mafft.bat --op 1.8 --ep 0.4 --maxiterate 15 --thread 4 /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_pol_140806_n557.fasta > /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_pol_140806_n557.fasta.mafft
	#tmp				<- hivc.cmd.clustalo(indir, infile.fasta.pol, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#
	#	located start of ENV around 5918			which( seq[nrow(seq), 4597:ncol(seq)]=='-' )		which( seq[nrow(seq), 5918:ncol(seq)]=='-' )	-> 660 -> 9058
	#	located end of ENV around 9176 
	file			<- paste(indir, '/fixup1_', infile.fasta,sep='')
	seq				<- read.dna(file, format='fasta')
	seq				<- as.character(seq)
	seq.env			<- seq[, 5850:9250]
	tmp				<- rownames(seq.env)
	seq.env			<- lapply(seq_len(nrow(seq.env)), function(i){	seq.env[i, seq.env[i,]!="-" ]	})	
	seq.env[[length(seq.env)+1]]	<- hxb2[, as.character(HXB2.K03455)]
	seq.env[[length(seq.env)+1]]	<- seq.env.c['CONSENSUS_C',]	
	names(seq.env)	<- c(tmp,'HXB2','CONSENSUS_C')
	seq.env			<- as.DNAbin( seq.env )
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta')
	#	align ENV
	tmp			<- hivc.cmd.clustalo(indir, infile.fasta.env, signat='', outdir=indir, prog= PR.CLUSTALO, hmm=PR.CLUSTALO.HMM, nproc=1, verbose=1)
	#
	#	process GAG
	#	
	file		<- paste(indir, '/', infile.fasta.gag,'.mafft',sep='')
	seq.gag		<- read.dna(file, format='fasta')
	seq.gag		<- as.character(seq.gag)
	seq.gag		<- seq.gag[,790:2698]
	tmp			<- rownames(seq.gag)
	seq.gag		<- lapply(seq_len(nrow(seq.gag)), function(i){	seq.gag[i, seq.gag[i,]!="-" ]	})
	names(seq.gag)	<- tmp
	seq.gag			<- as.DNAbin(seq.gag)
	infile.fasta.gag<- 'PANGEA_SSAfg_gag2_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')	
	#/Users/Oliver/git/hivclust/pkg/inst/mafft-mac/mafft.bat --op 3.0 --ep 1.5 --maxiterate 1000 --thread 4 /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_gag2_140806_n557.fasta > /Users/Oliver/duke/2014_Gates/methods_comparison_rootseqsim/140806/PANGEA_SSAfg_gag2_140806_n557.fasta.mafft
	#	manual edits
	file			<- paste(indir, '/', infile.fasta.gag,'.mafft',sep='')
	seq.gag			<- read.dna(file, format='fasta')
	seq.gag			<- as.character(seq.gag)
	seq.gag[, c(1361:1404, 1406:1413, 1421:1458)]	<- '-'
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=1)
	infile.fasta.gag<- 'PANGEA_SSAfg_gag3_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')
	#	final without HXB2 / CONSENSUS	
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=1)	#len 1466 -- 3 more than expected; there is an AA insertion for the SA seqs at pos 367
	seq.gag			<- seq.gag[ -c( which(grepl('HXB2',rownames(seq.gag))), which(grepl('CONSENSUS',rownames(seq.gag))) ),  ]
	seq.gag			<- seq.rmgaps(seq.gag, rm.only.col.gaps=1, verbose=1)
	infile.fasta.gag<- 'PANGEA_SSAfg_gag_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.gag,sep='')
	write.dna(seq.gag, file=file, format='fasta')	
	#
	#	process POL
	#	
	file			<- paste(indir, '/', infile.fasta.pol,'.mafft',sep='')
	seq.pol			<- read.dna(file, format='fasta')
	#	start: 2253		end: 5096	
	seq.pol			<- seq.pol[,2253:5096]
	seq.pol			<- seq.pol[ -c( which(grepl('HXB2',rownames(seq.pol))), which(grepl('CONSENSUS',rownames(seq.pol))) ),  ]
	infile.fasta.pol<- 'PANGEA_SSAfg_pol_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.pol,sep='')
	write.dna(seq.pol, file=file, format='fasta')
	#
	#	process ENV
	#
	#	start 6216	end 9734
	file		<- paste(indir, '/', infile.fasta.env,'.mafft',sep='')
	seq.env		<- read.dna(file, format='fasta')
	seq.env		<- as.character(seq.env)
	tmp			<- names(seq.env)
	seq.env		<- t(sapply( seq_along(seq.env), function(i){	seq.env[[i]][ 6216:9734 ]	}))
	seq.env		<- lapply(seq_len(nrow(seq.env)), function(i){	seq.env[i, seq.env[i,]!="-" ]	})	
	names(seq.env)	<- tmp
	seq.env			<- as.DNAbin(seq.env)
	infile.fasta.env<- 'PANGEA_SSAfg_env2_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta')	
	#	align on LosAlamos + manual edits
	file			<- paste(indir, '/', infile.fasta.env,'.hivalign',sep='')
	seq.env			<- read.dna(file, format='fasta')	
	seq.env			<- as.character(seq.env)
	seq.env[,c(1613:1645,2852:2932) ]	<- '-'  
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=1)	
	infile.fasta.env<- 'PANGEA_SSAfg_env3_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
	#	align on LosAlamos + manual edits
	file			<- paste(indir, '/', infile.fasta.env,'.hivalignnt',sep='')
	seq.env			<- read.dna(file, format='fasta')	
	seq.env			<- as.character(seq.env)
	tmp				<- max(sapply(seq.env, length))
	seq.env			<- t(sapply(seq.env, function(x) c(x, rep('-',tmp-length(x))) 	))
	seq.env[,c(430:543) ]	<- '-'		#rm non-consensus jitter
	seq.env[,c(1294:1344) ]	<- '-'		#rm non-consensus jitter
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=1)
	infile.fasta.env<- 'PANGEA_SSAfg_env4_140806_n557.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
	#seq.env[,c(1:9) ]	<- '-'			#C consensus start
	#!! There are two extra AA that are not in consensus C TAGTAG at pos 561
	#!! There are two extra AA that are not in consensus C at pos 2310
	#fixed some final issues manually
	seq.env			<- read.dna(file, format='fasta')
	seq.env 		<- seq.env[ -c( which(grepl('HXB2',rownames(seq.env))), which(grepl('CONSENSUS',rownames(seq.env))) ),  ]
	seq.env			<- seq.rmgaps(seq.env, rm.only.col.gaps=1, verbose=1)
	infile.fasta.env<- 'PANGEA_SSAfg_env_140806_n556_final.fasta'
	file			<- paste(indir, '/', infile.fasta.env,sep='')
	write.dna(seq.env, file=file, format='fasta', colsep='', nbcol=-1)
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

