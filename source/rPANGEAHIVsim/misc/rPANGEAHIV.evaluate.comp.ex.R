##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate<- function()
{
	require(RColorBrewer)
	dfa		<- project.PANGEA.TEST.pipeline.Aug2015.evaluate.read()
	#	check for updated submissions, and keep 
	dfa		<- project.PANGEA.TEST.pipeline.Aug2015.keep.most.recent.submission(dfa, format='%d.%m.%Y')
	#	save submissions
	outdir	<- '~/Dropbox (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_08_results/results'
	save(dfa, file=paste(outdir,'/submissions.R',sep=''))
	load(paste(outdir,'/submissions.R',sep=''))
	#	add objective legend
	dfa		<- merge(dfa, data.table(USED_GENES=c('pol','all'), USED_GENES_L=c('pol gene','pol+gag+env\ngenome') ), by='USED_GENES')
	set(dfa, NULL, 'TEAM', dfa[, factor(TEAM)])
	tmp		<- data.table( 	OBJ=	c('OBJ_i','OBJ_ii','OBJ_iii','OBJ_iv','OBJ_v','OBJ_vi'),
			OBJ_L=	c('Incidence\nTrend', '%Incidence', 'Incidence\nreduction', '%Acute Ctgr\n(baseline)', '%Acute\n(baseline)', '%Acute\n(endpoint)'))
	set(tmp, NULL, 'OBJ_L2', tmp[, factor(OBJ_L, levels=OBJ_L, labels=OBJ_L)])
	set(tmp, NULL, 'OBJ_L', tmp[, factor(OBJ_L, levels=rev(OBJ_L), labels=rev(OBJ_L))])
	dfa		<- merge(dfa, tmp, by='OBJ')
	#	add data legend
	dfa[, DATA_T2:='NA_character_']
	set(dfa, dfa[, which(DATA_T=='seq')], 'DATA_T2', 'using\nsequences')
	set(dfa, dfa[, which(DATA_T=='phy')], 'DATA_T2', 'using\ntrue tree')
	set(dfa, NULL, 'DATA_T2', dfa[, factor(DATA_T2, levels=rev(c('using\nsequences','using\ntrue tree')), labels=rev(c('using\nsequences','using\ntrue tree')))])		
	#	add scenario type
	set(dfa, NULL, 'DATA_T', dfa[, factor(DATA_T, levels=c('seq','phy'), labels=c('seq','phy'))])
	set(dfa, NULL, 'INT_T', dfa[, factor(INT_T, levels=c('fast','slow','none'), labels=c('fast','slow','none'))])
	set(dfa, NULL, 'AC_T', dfa[, factor(AC_T, levels=c('low','high'), labels=c('low','high'))])
	set(dfa, NULL, 'IMPRT', dfa[, factor(IMPRT*100, levels=c(5,20,2,0), labels=paste(c(5,20,2,0),'%',sep=''))])
	set(dfa, NULL, 'SMPL_C', dfa[, factor(SMPL_C*100, levels=c(8, 16, 30, 60), labels=paste(c(8, 16, 30, 60),'%',sep=''))])
	set(dfa, NULL, 'SMPL_D', dfa[, factor(SMPL_D, levels=c(5,3), labels=c(5,3))])	
	set(dfa, dfa[, which(SMPL_M=='overs')], 'SMPL_M', 'much')
	set(dfa, dfa[, which(SMPL_M=='extrs')], 'SMPL_M', 'extreme')
	set(dfa, dfa[, which(is.na(SMPL_M))], 'SMPL_M', 'extreme')
	set(dfa, NULL, 'SMPL_M', dfa[, factor(SMPL_M, levels=c('much','extreme'), labels=c('much','extreme'))])	
	tmp		<- unique(subset( dfa, select=c(DATAT_L, SC_RND, DATA_T, SC, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D) ))
	setkey(tmp, DATAT_L, AC_T, INT_T, DATA_T, IMPRT, SMPL_C, SMPL_D, SMPL_M)
	tmp[, SCENARIO_L:= paste('%AC=',AC_T,' ARTup=',INT_T,' EXT=',IMPRT,'\n',DATA_T,' ',SMPL_N,' ',SMPL_C,' ',SMPL_D,' ',SMPL_M, ' (',SC_RND,')',sep='')]
	dfa		<- merge(dfa, subset(tmp, select=c(SC_RND, SCENARIO_L)), by='SC_RND')
	#	add intervention legend
	dfa[, INT_L:= dfa[, paste('ART scale up\n',as.character(INT_T),sep='')]]
	setkey(dfa, INT_T)
	set(dfa, NULL, 'INT_L', dfa[, factor(INT_L, levels=dfa[, unique(INT_L)], labels=dfa[, unique(INT_L)])])
	#	add %Acute legend
	dfa[, AC_L:= dfa[, paste('%Acute\n',as.character(AC_T),sep='')]]
	setkey(dfa, AC_T)
	set(dfa, NULL, 'AC_L', dfa[, factor(AC_L, levels=dfa[, unique(AC_L)], labels=dfa[, unique(AC_L)])])
	#	get data.table of data sets ~ all primary and secondary objectives
	dfd		<- subset(dfa, select=c(SC_RND, DATA_T, DATAT_L, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C,  SMPL_M, SMPL_D))
	setkey(dfd, DATAT_L, INT_T, AC_T, DATA_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)
	dfd		<- unique(dfd)
	#	Primary Objectives, on sequences 
	tmp		<- data.table(expand.grid(ANA='Pr_Seq', SC_RND=subset(dfd, DATA_T=='seq')[, SC_RND], stringsAsFactors=FALSE))
	dfr		<- copy(tmp)
	#	Primary Objectives, on trees
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & (DATAT_L=='Regional' & IMPRT=='5%' & SMPL_M=='much' & SMPL_N==1600 | DATAT_L=='Village' & SMPL_C=='30%')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Pr_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: sequence coverage
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & (DATAT_L=='Regional' & IMPRT=='5%' & SMPL_M=='much' | DATAT_L=='Village')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SeqCoverage_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: imports
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & AC_T=='high' & (DATAT_L=='Regional' & SMPL_M=='much' & SMPL_C=='8%')  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_Imports_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: focussed sampling
	tmp		<- subset(dfd, DATA_T=='phy' & SMPL_D==5 & INT_T!='none' & AC_T=='low' & SMPL_C=='8%' & IMPRT=='5%' & DATAT_L=='Regional'  )
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SmplFc_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	Secondary: sampling duration
	tmp		<- subset(dfd, DATA_T=='phy' & INT_T!='none' & DATAT_L=='Regional' & SMPL_M=='much' & SMPL_C=='8%' & IMPRT=='5%' & INT_T=='fast')
	dfr		<- rbind(dfr, data.table(expand.grid(ANA='Sc_SmplD_Phy', SC_RND=tmp[, SC_RND], stringsAsFactors=FALSE)))
	#	merge with dfa
	dfr		<- dcast.data.table(dfr, SC_RND~ANA, value.var='SC_RND')
	set(dfr, NULL, 'Pr_Phy', dfr[, as.numeric(!is.na(Pr_Phy))])	
	set(dfr, NULL, 'Pr_Seq', dfr[, as.numeric(!is.na(Pr_Seq))])	
	set(dfr, NULL, 'Sc_Imports_Phy', dfr[, as.numeric(!is.na(Sc_Imports_Phy))])
	set(dfr, NULL, 'Sc_SeqCoverage_Phy', dfr[, as.numeric(!is.na(Sc_SeqCoverage_Phy))])
	set(dfr, NULL, 'Sc_SmplD_Phy', dfr[, as.numeric(!is.na(Sc_SmplD_Phy))])
	set(dfr, NULL, 'Sc_SmplFc_Phy', dfr[, as.numeric(!is.na(Sc_SmplFc_Phy))])	
	dfr		<- merge(unique(subset(dfa, select=c(DATAT_L,SC_RND))), dfr, by='SC_RND')
	dfa		<- merge(dfa, dfr, by=c('SC_RND','DATAT_L'))	
	#
	#	set team color
	#
	tmp				<- c('Cambridge','Cambridge/Imperial','ETH Zurich','Imperial','Vancouver','True','Cambridge/Imperial (chronos)','Cambridge/Imperial (lsd)','Cambridge/Imperial (mh15)','Cambridge/Imperial (mh30)','Cambridge/Imperial (merged)')
	TEAM_CL			<- c( brewer.pal(5, 'Set1'), 'black', brewer.pal(5, 'Set2') )
	names(TEAM_CL)	<- tmp	
	#	write info on scenario IDs to table
	dfi<- subset( dfa, select=c(SC_RND, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D) )
	setkey(dfi, DATAT_L, DATA_T, AC_T, INT_T, IMPRT, SMPL_N, SMPL_C, SMPL_M, SMPL_D)
	dfi	<- unique(dfi)
	file<- paste(outdir,'/SC_RND_info.csv',sep='')
	write.csv(dfi, file=file, row.names=FALSE)
	#
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.overallnumbers(dfa, outdir)
	#
	#	no results on pol submitted, focus on full genome only. get sample sizes per objective.
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.samplesize(dfr, dfa, outdir)
	#
	tmp	<- subset(dfa, select=c(TEAM, DATA_T, DATAT_L, SIM_SCENARIO))
	setkey(tmp, TEAM, DATA_T, DATAT_L, SIM_SCENARIO)
	tmp	<- unique(tmp)
	tmp[, table(TEAM, DATA_T, DATAT_L )]
	#	for each primary objective
	#	compare results across teams
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence(dfa, outdir, onSeq=1)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence(dfa, outdir, onSeq=0)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute(dfa, outdir, onSeq=1)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute(dfa, outdir, onSeq=0)
	#	for each secondary objective
	#	compare results across teams	
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.imports(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.focussedsampling(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.sduration(dfa, outdir)
	project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.seqcoverage(dfa, outdir)
	#	for each team
	#	all results
	invisible(sapply(setdiff(dfa[, unique(TEAM)],'True'), function(x)
					{		
						#x	<- 'Imperial'
						df		<- subset(dfa, (TEAM=='True' | TEAM==x) & USED_GENES=='all')
						set(df, df[, which(TEAM==x)], 'TEAM', 'estimate')
						set(df, df[, which(TEAM=='True')], 'TEAM', 'true value')
						set(df, NULL, 'TEAM', df[, factor(TEAM, levels=c('estimate','true value'), labels=c('estimate','true value'))])
						ggplot(df, aes(y=SCENARIO_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM, pch=TEAM)) + 
								geom_errorbarh(height=0.3) + geom_point(size=3) + 
								scale_colour_manual(values = c("red","black")) +
								scale_shape_manual(values = c(13,18), guide = FALSE) +
								labs(x='', y='', title= paste('TEAM',x,'\n'), colour='')  +
								facet_grid(DATAT_L~OBJ_L2, scales='free', space='free_y') +
								theme_bw() + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/res_obj_TEAM_',gsub(' ','_',gsub('\\/|\\(|\\)','',x)),'.pdf',sep=''), w=14, h=0.5*df[, length(unique(SCENARIO_L))])
						#	results using seq data
						df		<- subset(dfa, (TEAM=='True' | TEAM==x) & USED_GENES=='all' & DATA_T=='seq')
						set(df, df[, which(TEAM==x)], 'TEAM', 'estimate')
						set(df, df[, which(TEAM=='True')], 'TEAM', 'true value')
						set(df, NULL, 'TEAM', df[, factor(TEAM, levels=c('estimate','true value'), labels=c('estimate','true value'))])
						ggplot(df, aes(y=SCENARIO_L, x=central, xmin=lower95, xmax=upper95, colour=TEAM, pch=TEAM)) + 
								geom_errorbarh(height=0.3) + geom_point(size=3) + 
								scale_colour_manual(values = c("red","black")) +
								scale_shape_manual(values = c(13,18), guide = FALSE) +
								labs(x='', y='', title= paste('TEAM',x,'\n'), colour='')  +
								facet_grid(DATAT_L~OBJ_L2, scales='free', space='free_y') +
								theme_bw() + theme(legend.position='bottom')
						ggsave(file=paste(outdir,'/res_objonseq_TEAM_',gsub(' ','_',gsub('\\/|\\(|\\)','',x)),'.pdf',sep=''), w=14, h=0.7*df[, length(unique(SCENARIO_L))])	
					}))
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.overallnumbers<- function(dfa, outdir)
{
	#	count total submissions primary vs secondary
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1))
	tmp		<- tmp[, list(	Village=length(which(grepl('Vill',SIM_SCENARIO))), Regional=length(which(grepl('Regional',SIM_SCENARIO)))), by=c('TEAM','OBJ_L','USED_GENES_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(USED_GENES_L~variable) +			
			guides(fill=guide_legend(ncol=2)) +
			scale_fill_manual(values=TEAM_CL) +
			labs(x='', y='submissions\n(#)', title='Total scenarios submitted\n(using sequence data or true trees)\n') +
			theme_bw()+ theme(legend.position='bottom') + coord_flip()	
	ggsave(file=paste(outdir,'/res_scenarios_total.pdf',sep=''), w=10, h=8)
	
	#	count all submissions for primary objectives
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1) & DATA_T=='seq')
	tmp		<- tmp[, list(	Village=length(which(grepl('Vill',SIM_SCENARIO))), Regional=length(which(grepl('Regional',SIM_SCENARIO)))), by=c('TEAM','OBJ_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(~variable) +
			labs(x='', y='submissions\n(#)', title='Total scenarios submitted\n(using sequence data)\n') +
			scale_fill_manual(values=TEAM_CL) +
			guides(fill=guide_legend(ncol=2)) +
			theme_bw() + theme(legend.position='bottom') + coord_flip()
	ggsave(file=paste(outdir,'/res_scenarios_total_seqonly.pdf',sep=''), w=10, h=5)
	
	#	count complete submissions for primary objectives
	tmp		<- subset(dfa, TEAM!='True' & !grepl('(', TEAM, fixed=1) & DATA_T=='seq')
	tmp		<- tmp[, list(	Village=as.numeric(length(setdiff(c('01','02','03','04'),SC_RND))==0), Regional=as.numeric(length(setdiff(c('A','B','C','D'),SC))==0)), by=c('TEAM','OBJ_L','USED_GENES_L')]	
	tmp		<- melt(tmp, measure.vars=c('Village','Regional'))	
	ggplot(tmp, aes(x=OBJ_L, y=value, fill=TEAM)) + geom_bar(stat='identity') +
			facet_grid(USED_GENES_L~variable) +
			scale_y_continuous(breaks=seq(1,10,1), minor_breaks=NULL) +
			scale_fill_manual(values=TEAM_CL) +
			labs(x='', y='complete set of 4 submissions\n(#)', title='Complete submissions to evalute primary objectives\n(either village or regional)') +
			guides(fill=guide_legend(ncol=2)) +
			theme_bw() + theme(legend.position='bottom') + coord_flip()
	ggsave(file=paste(outdir,'/res_scenarios_total_seqonlycomplete.pdf',sep=''), w=10, h=7)
	
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.keep.most.recent.submission<- function(dfa, format='%d.%m.%Y')
{
	
	tmp		<- dfa[, list(SUB_N=length(unique(SUBMISSION_DATE)), SUB_DATE=unique(SUBMISSION_DATE)), by=c('TEAM','DATAT_L','OBJ')]
	tmp		<- subset(tmp, SUB_N>1)
	set(tmp, NULL, 'SUB_DATE', tmp[,as.Date(SUB_DATE, format=format)])
	#	for each objective, determine submissions that are to be discarded
	tmp		<- tmp[, list(SUB_DATE=SUB_DATE[SUB_DATE!=max(SUB_DATE)]), by=c('TEAM','DATAT_L','OBJ')]
	for(i in seq_len(nrow(tmp)))
		set(dfa, dfa[, which(TEAM==tmp$TEAM[i] & DATAT_L==tmp$DATAT_L[i] & OBJ==tmp$OBJ[i] & SUBMISSION_DATE==tmp[, as.character(tmp$SUB_DATE[i], format=format)])], 'central', NA_real_)
	subset(dfa, !is.na(central))	
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.incidence<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	#,'OBJ_v','OBJ_vi'
	#	compare objectives with / without seq data, village + regional	
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\nIncidence from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nPrimary objective\nIncidence when phylogeny known\n'
	}			
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p1		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 10)) +
			scale_x_continuous(breaks=seq(1,9,1)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n%Incidence', y='')
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p2		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			coord_cartesian(xlim=c(0, 2)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.2,1.8,0.4), minor_breaks=seq(0.2,1.8,0.2)) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='')
	pdf(file=paste(outdir,'/res_acrossTEAM_PrimaryIncidence_onSeq',onSeq,'.pdf',sep=''), width = 15, height = 8)
	print(grid.arrange(p1, p2, nrow=1, main=title))
	dev.off() 
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute2<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\n%Acute from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nSecondary objective\n%Acute when phylogeny known\n'
	}			
	tmp		<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))	
	tmp3	<- subset(tmp, TEAM=='True', select=c(SC_RND, central))
	setnames(tmp3, 'central', 'central_true')
	tmp		<- subset(merge(tmp, tmp3, by='SC_RND', all.x=TRUE), TEAM!='True')
	
	ggplot(tmp, aes(x=central_true, y=central, colour=TEAM, pch= gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')))) +
			geom_point() +
			geom_abline(slope=1, intercept=0) +
			scale_colour_manual(values=TEAM_CL) +
			facet_wrap(~DATAT_L, scales='free') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3), pch=guide_legend(ncol=3)) +
			labs(x= '\ntrue % Acute at baseline', y='\nestimated % Acute at baseline', pch='')	
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.primary.acute<- function(dfa, outdir, onSeq=1)
{
	require(gridExtra)
	if(onSeq)
	{
		df		<- subset(dfa, Pr_Seq==1)
		title	<- '\nPrimary objective\n%Acute from sequence data\n'
	}		
	if(!onSeq)
	{
		df	<- subset(dfa, Pr_Phy==1)
		title	<- '\nSecondary objective\n%Acute when phylogeny known\n'
	}			
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% transmissions from individuals in their first 3 months of infection\nat baseline', y='')
	tmp	<- subset(df, !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')		
	p2		<- ggplot(subset(tmp, TEAM!='True'), aes(y=gsub('\n',':',paste(INT_L,'  ',AC_L,sep='')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +			
			scale_colour_manual(values=TEAM_CL) +			
			facet_grid(TEAM~DATAT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% transmissions from individuals in their first 3 months of infection\nat end of intervention', y='')
	pdf(file=paste(outdir,'/res_acrossTEAM_PrimaryAcute_onSeq',onSeq,'.pdf',sep=''), width = 15, height = 8)
	print(grid.arrange(p1, p2, nrow=1, main=title))
	dev.off() 
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.seqcoverage<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of sequence coverage\n'
	
	tmp		<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp2	<- as.data.table(expand.grid(SMPL_C=tmp[, unique(SMPL_C)], central=1, DATAT_L=tmp[, unique(DATAT_L)], INT_L=tmp[, unique(INT_L)], AC_L=tmp[,unique(AC_L)]))	
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p1		<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='sequence coverage\n')
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp2	<- as.data.table(expand.grid(SMPL_C=tmp[, unique(SMPL_C)], central=1, DATAT_L=tmp[, unique(DATAT_L)], INT_L=tmp[, unique(INT_L)], AC_L=tmp[,unique(AC_L)]))
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 2)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,1.5,0.5), minor_breaks=seq(0.25,1.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='sequence coverage\n')
	
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='sequence coverage\n')
	
	
	tmp	<- subset(dfa, Sc_SeqCoverage_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=SMPL_C, x=central)) +
			geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='sequence coverage\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SeqCoverage.pdf',sep=''), width = 15, height = 30)
	print(grid.arrange(p1, p2, p3, p4, nrow=4, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.imports<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of 20% / year transmissions from outside\n'
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='transmissions from outside\n')
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='transmissions from outside\n')
	
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='transmissions from outside\n')
	
	
	tmp	<- subset(dfa, Sc_Imports_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(IMPRT,' / year',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='transmissions from outside\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_Imports.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.focussedsampling<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of 50% vs 85% of samples obtained after intervention start\n'
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='% samples obtained after\nintervention start\n')
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='% samples obtained after\nintervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='% samples obtained after\nintervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplFc_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'INT_L',tmp[,gsub('\n',': ',INT_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, INT_L, SMPL_M)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=factor(as.character(SMPL_M),levels=c('much','extreme'),labels=c('50%','85%')), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+INT_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='% samples obtained after\nintervention start\n')
	
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SFocus.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.secondary.sduration<- function(dfa, outdir)
{
	title	<- '\nSecondary objective\nImpact of sampling duration 3 yrs vs 5 yrs\n'
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_ii'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p1	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 8)) +
			scale_x_continuous(breaks=seq(1,7,1)) +			
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% Incidence', y='duration of intensified sampling\nafter intervention start\n')
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_iii'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p2	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			geom_vline(xintercept=1, colour='grey50', size=0.8) +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) +
			coord_cartesian(xlim=c(0, 4)) +
			scale_colour_manual(values=TEAM_CL) +
			scale_x_continuous(breaks=seq(0.5,3.5,0.5), minor_breaks=seq(0.25,3.75,0.25)) +
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +			
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\nIncidence reduction', y='duration of intensified sampling\nafter intervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_v'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p3	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat baseline', y='duration of intensified sampling\nafter intervention start\n')
	
	
	tmp	<- subset(dfa, Sc_SmplD_Phy==1 & !grepl('(',TEAM,fixed=1) & USED_GENES=='all' & OBJ%in%c('OBJ_vi'))
	set(tmp, NULL, 'AC_L',tmp[,gsub('\n',': ',AC_L)])
	tmp3	<- subset(tmp, TEAM=='True')
	setnames(tmp3, 'TEAM','team')	
	setkey(tmp, TEAM, AC_L, IMPRT)	
	p4	<- ggplot(subset(tmp, TEAM!='True'), aes(y=paste(SMPL_D,' yrs',sep=''), x=central)) +
			#geom_point(data=tmp2, size=1, colour='transparent') +
			geom_errorbarh(aes(xmin=lower95, xmax=upper95, colour=TEAM), height=0.3) + 
			geom_point(size=4, aes(colour=TEAM), pch=13) +			
			geom_point(data=tmp3, size=3, colour='black', pch=18) +
			coord_cartesian(xlim=c(0, 50)) +
			scale_x_continuous(breaks=seq(10,40,10), minor_breaks=seq(0,50,2)) +
			scale_colour_manual(values=TEAM_CL) +
			facet_grid(TEAM~DATAT_L+AC_L, scales='free', space='free_y') +
			theme_bw() + theme(panel.margin.x= unit(0.1, "lines"), legend.position='bottom') +
			guides(colour=guide_legend(ncol=3)) +
			labs(x= '\n% trms from individuals in their first 3 months of infection\nat end of intervention', y='duration of intensified sampling\nafter intervention start\n')
	
	pdf(file=paste(outdir,'/res_acrossTEAM_Secondary_SDuration.pdf',sep=''), width = 12, height = 12)
	print(grid.arrange(p1, p2, p3, p4, nrow=2, main=title))
	dev.off()
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 08.05.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.read<- function()
{
	#	read truth for regional simus	
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_05_results'	
	file	<- paste(indir, '/answers_Regional_Feb2015_rFormat.csv', sep='')
	df		<- read.submission.Feb2015(file, verbose=0, reset.OBJiv.conservative=1)
	#	read truth for village simus
	file	<- paste(indir, '/answers_Village_Feb2015-yr43_rFormat.csv', sep='')
	tmp		<- read.submission.Feb2015(file, verbose=0, reset.OBJiv.conservative=1)
	set(tmp, NULL, 'TEAM', 'True')
	df		<- rbind(df, tmp)
	#	read submissions from May 2015
	tmp		<- list.files(indir, pattern='csv$')
	tmp		<- tmp[!grepl('answers',tmp)]
	#	read Eriks multiple submissions from May 2015
	tmp2	<- data.table(FILE=tmp[grepl('cambImp',tmp)])
	tmp2[, RUN:= tmp2[,  sapply( strsplit(FILE,'_'), function(x) rev(x)[1] )]]
	set(tmp2, NULL, 'RUN', tmp2[, substr(RUN, 1, nchar(RUN)-4)])
	set(tmp2, NULL, 'RUN', tmp2[, gsub('results0','',RUN)])
	dfs		<- do.call('rbind',lapply(seq_len(nrow(tmp2)), function(i)
					{
						z	<- read.submission.Feb2015( paste(indir, '/', tmp2[i, FILE], sep=''), verbose=0, reset.OBJiv.conservative=1 )
						set(z, NULL, 'TEAM', z[, paste(TEAM, ' (', tmp2[i, RUN], ')', sep='')])
						z
					}))
	tmp		<- tmp[!grepl('cambImp',tmp)]
	tmp		<- do.call('rbind',lapply(tmp, function(x) read.submission.Feb2015(paste(indir,'/',x,sep=''), verbose=0, reset.OBJiv.conservative=1)))
	dfs		<- rbind(dfs, tmp)
	#	read submissions from August 2015
	indir	<- '/Users/Oliver/Dropbox\ (Infectious Disease)/PANGEAHIVsim_internal/documents/external/2015_08_results'
	tmp		<- list.files(indir, pattern='csv$')
	tmp		<- tmp[!grepl('answers',tmp)]
	tmp2	<- tmp[grepl('Vancouver',tmp)]
	stopifnot(length(tmp2)==1)
	tmp2	<- read.submission.May2015(paste(indir,'/',tmp2,sep=''), verbose=0)
	dfs		<- rbind(dfs, tmp2)
	tmp2	<- tmp[!grepl('Vancouver',tmp)]
	tmp2	<- do.call('rbind',lapply(tmp2, function(x) read.submission.Aug2015(paste(indir,'/',x,sep=''), verbose=0, reset.OBJiv.conservative=1)))
	dfs		<- rbind(dfs, tmp2)
	# 	change team name
	set(dfs, dfs[, which(TEAM=='Colijn')],'TEAM','Imperial')
	#	construct Erik's gold submission
	#	for regional tree, use mergedTab
	tmp		<- subset(dfs, grepl('merged', TEAM) & grepl('Regional',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']	
	#tmp		<- subset(dfs, grepl('mh30', TEAM) & grepl('Regional',SIM_SCENARIO))	
	#tmp[, TEAM:='Cambridge/Imperial']
	#tmp2	<- subset(dfs, grepl('mh15', TEAM) & grepl('Regional',SIM_SCENARIO))	
	#tmp2[, TEAM:='Cambridge/Imperial']
	#tmp		<- merge(tmp, tmp2, by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES','OBJ'), all=1)
	#tmp2	<- tmp[, which(is.na(central.x))]
	#set(tmp, tmp2, 'central.x', tmp[tmp2, central.y])
	#set(tmp, tmp2, 'lower95.x', tmp[tmp2, lower95.y])
	#set(tmp, tmp2, 'upper95.x', tmp[tmp2, upper95.y])
	#setnames(tmp, c('central.x', 'lower95.x', 'upper95.x'), c('central', 'lower95', 'upper95'))
	#set(tmp, NULL, c('central.y', 'lower95.y', 'upper95.y'), NULL)
	dfs		<- rbind(dfs, tmp)
	#	for village tree, use mh30 where available and mh15 where mh30 not available
	tmp		<- subset(dfs, grepl('mh30', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']
	tmp2	<- subset(dfs, grepl('mh15', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp2[, TEAM:='Cambridge/Imperial']
	tmp		<- merge(tmp, tmp2, by=c('TEAM','SUBMISSION_DATE','SIM_SCENARIO','USED_GENES','OBJ'), all=1)
	tmp2	<- tmp[, which(is.na(central.x))]
	set(tmp, tmp2, 'central.x', tmp[tmp2, central.y])
	set(tmp, tmp2, 'lower95.x', tmp[tmp2, lower95.y])
	set(tmp, tmp2, 'upper95.x', tmp[tmp2, upper95.y])
	setnames(tmp, c('central.x', 'lower95.x', 'upper95.x'), c('central', 'lower95', 'upper95'))
	set(tmp, NULL, c('central.y', 'lower95.y', 'upper95.y'), NULL)
	dfs		<- rbind(dfs, tmp)
	#	for village seq, use LSD
	tmp		<- subset(dfs, grepl('lsd', TEAM) & grepl('Vill',SIM_SCENARIO))	
	tmp[, TEAM:='Cambridge/Imperial']
	dfs		<- rbind(dfs, tmp)
	#	define data types (seq or phylo)
	dfa		<- rbind(dfs, df)
	dfa[, DATA_T:=NA_character_]
	set(dfa, dfa[, which(grepl('Vill_0[1-4]', SIM_SCENARIO))], 'DATA_T', 'seq')
	set(dfa, dfa[, which(!grepl('Vill_0[1-4]', SIM_SCENARIO))], 'DATA_T', 'phy')	
	set(dfa, dfa[, which(grepl('FirstObj', SIM_SCENARIO))], 'DATA_T', 'seq')
	set(dfa, dfa[, which(grepl('SecondObj', SIM_SCENARIO))], 'DATA_T', 'phy')
	stopifnot(!any(is.na(dfa[, DATA_T])))
	#	define randomized scenario IDs
	dfa[, SC_RND:=NA_character_]
	tmp		<- dfa[, which(grepl('Regional',SIM_SCENARIO))]
	set(dfa, tmp, 'SC_RND', dfa[tmp, substring(regmatches(SIM_SCENARIO,regexpr('sc[A-Z]',SIM_SCENARIO)),3)])
	tmp		<- dfa[, which(grepl('Vill',SIM_SCENARIO))]
	set(dfa, tmp, 'SC_RND', dfa[tmp, substring(regmatches(SIM_SCENARIO,regexpr('Vill_[0-9]+',SIM_SCENARIO)),6)])
	stopifnot(!any(is.na(dfa[, SC_RND])))
	
	#	describe regional simulations in terms of fast/low intervention high/low acute	
	set.seed(42)
	dfi			<- data.table(FILE=list.files('/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/FINAL', '.*zip$', full.names=FALSE))	
	dfi[, SC:= sapply(strsplit(FILE, '_'),'[[',3)]
	dfi[, CONFIG:= sapply(strsplit(SC, '-'),'[[',2)]
	set(dfi, NULL, 'SC', dfi[, sapply(strsplit(SC, '-'),'[[',1)])
	dfi[, DATAT:= sapply(strsplit(FILE, '_'),'[[',5)]
	set(dfi, NULL, 'DATAT', dfi[, gsub('.zip','',DATAT,fixed=T)])	
	set(dfi, NULL, 'OBJECTIVE', 'SecondObj')
	set(dfi, dfi[,which(CONFIG=='sq')],'OBJECTIVE', 'FirstObj')
	dfi			<- merge(dfi,dfi[, list(FILE=FILE, DUMMY=sample(length(FILE),length(FILE))), by='OBJECTIVE'],by=c('OBJECTIVE','FILE'))
	tmp			<- dfi[, which(OBJECTIVE=='SecondObj')]
	set(dfi, tmp, 'DUMMY', dfi[tmp, DUMMY] + dfi[OBJECTIVE=='FirstObj', max(DUMMY)])	
	setkey(dfi, DUMMY)
	dfi[, SC_RND:= toupper(letters[seq_len(nrow(dfi))])]
	dfi			<- subset(dfi, select=c(SC, SC_RND, CONFIG))
	set(dfi, NULL, 'SC', dfi[, substring(SC, 3)])
	dfi			<- merge( dfi, data.table(SC= c('A','B','C','D','E','F'), AC_T=c('low','high','low','high','low','high'), INT_T=c('fast','fast','slow','slow','none','none')), by='SC' )
	tmp			<- data.table(	CONFIG=	c('sq','s2x','y3','mFP85','ph','tr20'),
			IMPRT=	c(.05, .05, .05, .05, .05, .2),
			SMPL_N=	c(1600, 3200, 1280, 1600, 1600, 1600),
			SMPL_C= c(0.08, 0.16, 0.08, 0.08, 0.08, 0.08),
			SMPL_M=	c('overs', 'overs', 'overs', 'extrs', 'overs', 'overs'),
			SMPL_D= c(5, 5, 3, 5, 5, 5))
	dfi			<- merge( dfi, tmp, by='CONFIG')					
	set(dfi, NULL, c('CONFIG'), NULL)
	#	add info for village
	tmp			<- data.table(	SC_RND= c('03','02','01','04','05','08','06','07','11','09','12','10','00'),
			AC_T=	c('low','low','high','high','low','low','high','high','low','low','high','high','low'),
			INT_T=	c('fast','slow','fast','slow','fast','slow','fast','slow','fast','slow','fast','slow','none'),
			#SMPL_C=	c(0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25, 0.25),
			SMPL_C=	c(0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 0.3, 0.3, 0.3, 0.3, 0.3),
			SMPL_D= 5,
			SMPL_N= c(777, 857, 957, 1040, 1469, 1630, 1831, 1996, 638, 686, 956, 1012, 872),
			IMPRT=	c(0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0, 0, 0, 0, 0.02))	
	dfi			<- rbind(dfi, tmp, fill=TRUE,use.names=TRUE)
	#	merge to dfa
	cat(paste('\nnumber of rows before merge with dfi, n=', nrow(dfa)))
	dfa			<- merge(dfa, dfi, by='SC_RND')
	cat(paste('\nnumber of rows before merge with dfi, n=', nrow(dfa)))
	#	set answers to numerical
	set(dfa, dfa[, which(OBJ%in%c('OBJ_i','OBJ_iv'))], c('lower95','upper95'), NA_character_)
	set(dfa, dfa[, which(central=='decreasing')], c('central'), '-1')
	set(dfa, dfa[, which(central=='stable')], c('central'), '0')
	set(dfa, dfa[, which(central=='increasing')], c('central'), '1')	
	set(dfa, dfa[, which(central=='<15%')], c('central'), '-1')
	set(dfa, dfa[, which(central=='15%-30%')], c('central'), '0')
	set(dfa, dfa[, which(central=='>30%')], c('central'), '1')
	set(dfa, NULL, 'central', dfa[, as.numeric(central)])
	set(dfa, NULL, 'lower95', dfa[, as.numeric(lower95)])
	set(dfa, NULL, 'upper95', dfa[, as.numeric(upper95)])	
	#	add simulation type
	dfa[, DATAT_L:='NA_character_']
	set(dfa, dfa[, which(grepl('Vill',SIM_SCENARIO))], 'DATAT_L','Village')
	set(dfa, dfa[, which(grepl('Regional',SIM_SCENARIO))], 'DATAT_L','Regional')
	#	transform %incidence 0.01 to 1%
	tmp		<- dfa[, which(OBJ=='OBJ_ii')]
	set(dfa, tmp, 'central', dfa[tmp, 100*central])
	set(dfa, tmp, 'lower95', dfa[tmp, 100*lower95])
	set(dfa, tmp, 'upper95', dfa[tmp, 100*upper95])	
	#	transform %Acute 0.01 to 1%
	tmp		<- dfa[, which(OBJ%in%c('OBJ_v','OBJ_vi'))]
	set(dfa, tmp, 'central', dfa[tmp, 100*central])
	set(dfa, tmp, 'lower95', dfa[tmp, 100*lower95])
	set(dfa, tmp, 'upper95', dfa[tmp, 100*upper95])	
	#	fix submission dates
	set(dfa, NULL, 'SUBMISSION_DATE', dfa[, gsub('\\.15','\\.2015',SUBMISSION_DATE)])
	dfa
}
##--------------------------------------------------------------------------------------------------------
##	evaluate results
##	olli 12.08.15
##--------------------------------------------------------------------------------------------------------
project.PANGEA.TEST.pipeline.Aug2015.evaluate.samplesize<- function(dfr, dfa, outdir)
{
	tmp		<- dfr[, lapply(.SD, sum ), .SDcol=c('Pr_Phy','Pr_Seq','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy'), by='DATAT_L']
	tmp2	<- dfa[, unique(TEAM[TEAM!='True' & !grepl('(',TEAM,fixed=T)])]
	tmp		<- melt(tmp, id.vars='DATAT_L')
	set(tmp, NULL, 'value', tmp[,value]*length(tmp2))
	tmp		<- dcast.data.table(tmp, DATAT_L~variable, value.var='value')	
	tmp[, N:='TOTAL']
	tmp2	<- unique(subset(dfa, select=c(DATAT_L,OBJ)))
	tmp2[, N:='TOTAL']
	tmp		<- merge(tmp, tmp2, by=c('DATAT_L','N'))	
	tmp2	<- subset(dfa, USED_GENES=='all' & TEAM!='True' & !grepl('(',TEAM,fixed=T))[, lapply(.SD, sum ), by=c('DATAT_L','OBJ','TEAM'), .SDcol=c('Pr_Phy','Pr_Seq','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy')]
	tmp2[, N:='SUBMITTED']	
	tmp		<- rbind(tmp, tmp2, use.names=T, fill=T)	
	tmp2	<- dcast.data.table(melt(tmp, id.vars=c('DATAT_L','OBJ','N','TEAM')), DATAT_L+OBJ~variable+N,fun.aggregate=sum,value.var='value')
	file	<- paste(outdir,'/SampleSizesByAnalysis_PolGagEnv','.csv',sep='')
	cat('\nWrite to',file)
	write.csv(tmp2, file=file, row.names=FALSE)
	
	set(tmp, tmp[, which(is.na(TEAM))], 'TEAM', 'True')
	z	<- melt(subset(tmp, OBJ=='OBJ_i'), id.vars=c('DATAT_L','N','OBJ','TEAM'))
	set(z, NULL, 'variable', z[, factor(variable, 	levels=c('Pr_Seq','Pr_Phy','Sc_Imports_Phy','Sc_SeqCoverage_Phy','Sc_SmplD_Phy','Sc_SmplFc_Phy'), 
							labels=c('Primary Objective\non sequences','Primary Objective\non true tree','Secondary\nimports','Secondary\nsequence coverage','Secondary\nsampling duration','Secondary\nfocussed sampling'))])
	ggplot(z, aes(x=factor(N, levels=c('SUBMITTED','TOTAL'), labels=c('submitted','if submissions\nhad been\ncomplete')), fill=TEAM, y=value)) + geom_bar(stat='identity') + facet_grid(DATAT_L~variable) +
			scale_fill_manual(values=TEAM_CL) +			
			theme_bw() + theme(legend.position='bottom') + labs(x='', y='total')
	file	<- paste(outdir,'/SampleSizesByAnalysis_PolGagEnv.pdf',sep='')
	cat('\nPlot to',file)
	ggsave(file=file, w=12, h=4)
	NULL	
}