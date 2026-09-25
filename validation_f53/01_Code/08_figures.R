invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
dir.create("02_Results/Figures",showWarnings=FALSE)
theme_set(theme_classic(base_size=10,base_family="Arial")+theme(axis.text=element_text(colour="black"),plot.title=element_text(size=11,face="bold"),legend.position="none",plot.margin=margin(10,15,10,10)))
save_figure<-function(p,name,height=4.5) {
 stem<-file.path("02_Results/Figures",name)
 ggsave(paste0(stem,".pdf"),p,width=7.2,height=height,device=cairo_pdf)
 ggsave(paste0(stem,".png"),p,width=7.2,height=height,dpi=600,bg="white")
 ggsave(paste0(stem,".tiff"),p,width=7.2,height=height,dpi=600,compression="lzw",bg="white")
 if(requireNamespace("svglite",quietly=TRUE))ggsave(paste0(stem,".svg"),p,width=7.2,height=height,device=svglite::svglite)
}
forest<-function(d,xlab,title) {
 d<-copy(d);d[,label:=paste0(label,"  (n=",format(n,big.mark=",",trim=TRUE),")")]
 d[,label:=factor(label,levels=rev(unique(label)))]
 ggplot(d,aes(x=effect,y=label))+geom_vline(xintercept=1,linetype="dashed",colour="#777777",linewidth=.4)+
 geom_errorbar(aes(xmin=lower,xmax=upper),width=.15,orientation="y",linewidth=.55,colour="#286D8E")+
 geom_point(size=2.4,colour="#286D8E")+scale_x_log10(breaks=c(.5,1,2,4))+
 labs(x=xlab,y=NULL,title=title)+theme(axis.text.y=element_text(size=9))
}
t<-fread("02_Results/V04_temporal_observed_models.csv")[term=="map_qQ4" & model %in% c("Earlier <=2016 Primary","Later >=2017 Primary","Earlier <=2016 Phenotype","Later >=2017 Phenotype")]
t[,`:=`(effect=HR,label=fcase(model=="Earlier <=2016 Primary","Earlier period: primary",model=="Later >=2017 Primary","Later period: primary",model=="Earlier <=2016 Phenotype","Earlier period: phenotype",default="Later period: phenotype"))]
save_figure(forest(t,"Q4 versus Q2 hazard ratio (95% CI)","Internal temporal stability: observed-hour MAP"),"Figure_V1_internal_temporal",3.4)
p<-fread("02_Results/P08_full_cohort_Q4_changes.csv")[model %in% c("Primary model","Plus ICU type","Plus admission type","Plus hospital service","F52 phenotype model","Extended clinical context","SOFA domain replacement","Noncardiovascular SOFA","Non-CNS SOFA")]
p[,`:=`(effect=HR,label=model)]
save_figure(forest(p,"Q4 versus Q2 hazard ratio (95% CI)","Clinical context and SOFA adjustment"),"Figure_P1_phenotype_models",4.8)
s<-fread("02_Results/V02_MI_subgroup_models.csv")[term=="map_qQ4" & model %in% c("Medical ICU","Cardiac ICU","Surgical/Trauma ICU","Neuro ICU/Stepdown","No MV","MV","No RRT","RRT")]
s[,`:=`(effect=HR,label=model)]
save_figure(forest(s,"Q4 versus Q2 hazard ratio (95% CI)","Clinical subgroup associations: m=40"),"Figure_V2_clinical_subgroups",4.6)
e<-fread("02_Results/E03_harmonized_hospital_models.csv")[term=="map_qQ4" & model %in% c("eICU shared primary","eICU plus clinical context","eICU complete24","eICU antibiotic culture SOFA proxy","MIMIC-IV aligned early-sepsis hospital endpoint","MIMIC-IV aligned complete24")]
e[,`:=`(effect=OR,label=fcase(model=="MIMIC-IV aligned early-sepsis hospital endpoint","MIMIC-IV: aligned hospital endpoint",model=="MIMIC-IV aligned complete24","MIMIC-IV: complete 24 hours",model=="eICU antibiotic culture SOFA proxy","eICU: antibiotic/culture SOFA proxy",default=model))]
save_figure(forest(e,"Fixed high versus reference MAP category odds ratio (95% CI)","Exploratory in-hospital mortality transportability"),"Figure_E1_hospital_transportability",4.1)
c<-fread("04_QC/mice_chain_means.csv")[hour %in% c("h00","h01")]
av<-c[,.(mean=mean(mean)),by=.(hour,iteration)]
q<-ggplot(c,aes(iteration,mean,group=chain))+geom_line(alpha=.3,linewidth=.3,colour="#888888")+
 geom_line(data=av,aes(iteration,mean,group=1),linewidth=.8,colour="#286D8E")+facet_wrap(~hour,ncol=1,scales="free_y")+
 scale_x_continuous(breaks=1:10)+labs(x="MICE iteration",y="Mean of imputed hourly MAP (mmHg)",title="Primary m=40 imputation traces")+
 theme(strip.background=element_blank(),strip.text=element_text(face="bold"))
save_figure(q,"Figure_QC1_imputation_traces",4.2)
fwrite(rbindlist(list(t[,.(figure="V1",label,effect,lower,upper,n)],p[,.(figure="P1",label,effect,lower,upper,n)],s[,.(figure="V2",label,effect,lower,upper,n)],e[,.(figure="E1",label,effect,lower,upper,n)])),"02_Results/Figures/figure_source_data.csv")
