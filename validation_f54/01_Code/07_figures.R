invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
dir.create("04_Figures",showWarnings=FALSE)
theme_set(theme_classic(base_size=10,base_family="Arial")+theme(legend.position="bottom",legend.title=element_blank(),plot.margin=margin(8,10,8,8),axis.text=element_text(color="black")))
savefig<-function(p,n,w=7,h=4.4){
 base<-file.path("04_Figures",n)
 ggsave(paste0(base,".pdf"),p,width=w,height=h,device=cairo_pdf)
 ggsave(paste0(base,".svg"),p,width=w,height=h,device=svglite::svglite)
 ggsave(paste0(base,".png"),p,width=w,height=h,dpi=600,device=ragg::agg_png,bg="white")
 ggsave(paste0(base,".tiff"),p,width=w,height=h,dpi=600,device=ragg::agg_tiff,compression="lzw",bg="white")
 ggsave(paste0(base,"_preview.png"),p,width=w,height=h,dpi=150,device=ragg::agg_png,bg="white")
}
curves<-fread("02_Results/spline_curves.csv")
curve<-function(d,xlab){
 ggplot(d,aes(x,HR,color=adjustment,fill=adjustment))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,colour=NA)+geom_line(linewidth=.8)+geom_hline(yintercept=1,linetype=2,color="gray40",linewidth=.4)+geom_vline(xintercept=unique(d$reference),linetype=3,color="gray40",linewidth=.4)+scale_color_manual(values=c(Primary="#236B72",Extended="#9A4664"))+scale_fill_manual(values=c(Primary="#236B72",Extended="#9A4664"))+labs(x=xlab,y="Hazard ratio (95% CI)")
}
savefig(curve(curves[exposure=="mean_map"],"Early 24-hour TWA-MAP (mmHg)"),"Figure2_MAP_spline",7,4.4)
savefig(curve(curves[exposure=="ttr65"],"Fraction of hourly MAP values below 65 mmHg"),"FigureS4_hypotension_spline",7,4.4)
forest<-function(d){
 d<-copy(d)
 d[,model:=gsub("No vasoactive infusion","No NE-equivalent vasopressor",model,fixed=TRUE)]
 d[,model:=gsub("Vasoactive infusion","NE-equivalent vasopressor",model,fixed=TRUE)]
 d[,model:=gsub("Plus vascular history","Plus vascular codes",model,fixed=TRUE)]
 d<-copy(d);d[,label:=factor(model,levels=rev(unique(model)))];d[,effect:=sprintf("%.2f (%.2f-%.2f)",HR,lower,upper)]
 ggplot(d,aes(HR,label))+geom_vline(xintercept=1,linetype=2,linewidth=.4,color="gray40")+geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=.13,linewidth=.55,color="#236B72")+geom_point(size=2.3,color="#236B72")+scale_x_log10()+labs(x="Q4 versus Q2 hazard ratio (95% CI)",y=NULL)+theme(axis.text.y=element_text(size=9))
}
m<-fread("02_Results/main_models.csv")[term=="map_qQ4"]
whichmodels<-c("Primary","Plus vascular history","Plus ICU type","Plus admission type","Plus service","Plus procedures/devices","Plus NE dose","Phenotype","Extended context","SOFA components")
savefig(forest(m[match(whichmodels,model)]),"Figure3_clinical_context",7,4.7)
t<-fread("02_Results/temporal_models.csv")[term=="map_qQ4" & status=="estimated"]
t[,model:=gsub("Earlier <=2016","Earlier period",model,fixed=TRUE)]
t[,model:=gsub("Later >=2017","Later period",model,fixed=TRUE)]
savefig(forest(t),"FigureS1_temporal",7,4.2)
s<-fread("02_Results/subgroup_models.csv")[term=="map_qQ4" & status=="estimated"]
savefig(forest(s),"FigureS2_clinical_strata",7,5)
cm<-fread("04_QC/chain_means.csv")[hour %in% c("h00","h01")]
p<-ggplot(cm,aes(iteration,mean,group=chain))+geom_line(alpha=.22,color="gray40",linewidth=.3)+stat_summary(aes(group=1),fun=mean,geom="line",linewidth=.8,color="#236B72")+facet_wrap(~hour,ncol=1,scales="free_y")+scale_x_continuous(breaks=c(1,10,20,30))+labs(x="Iteration",y="Mean imputed MAP (mmHg)")+theme(strip.background=element_blank())
savefig(p,"FigureS3_imputation_traces",7,5)
fwrite(rbindlist(list(m[,.(figure="Figure3",model,HR,lower,upper)],t[,.(figure="FigureS1",model,HR,lower,upper)],s[,.(figure="FigureS2",model,HR,lower,upper)])),"04_Figures/forest_source_data.csv")
fwrite(curves,"04_Figures/spline_source_data.csv")
cat("Six figures exported in PDF, SVG, 600 dpi PNG and TIFF\n")
