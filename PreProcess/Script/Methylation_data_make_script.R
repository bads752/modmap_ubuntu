source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')
source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_function.R')

#cl <- makeCluster(7)
#registerDoParallel(cl)

projects<-tcga_projects()
projects<-str_split(projects,'-')
projects<-lapply(projects,function(x){
  x[[2]]
})
projects<-unlist(projects)
#f<-list.files('//home/giordano/Methylation/Pre_Process/')
#f<-str_split(f,'_')
#f<-lapply(f,function(x){
#  x[[1]]
##})
#f<-unlist(f)
#projects<-projects[!projects%in%f]


CancerSite <- projects[as.numeric(commandArgs(TRUE)[1])]
CancerSite <- 'CHOL'
print(CancerSite)
TargetDirectory <- paste0("C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/Input_Data/Methylation/")

directories <- Download_DNAmethylation(CancerSite, TargetDirectory)



assign(paste0(CancerSite,'_Data'),Preprocess_DNAmethylation(CancerSite, 
directories,lg=F, mvt = 0.2, bct='limma',normal=T,svtmp='C:/Users/giordano/')) 

save(list = c(paste0(CancerSite,'_Data')),
    file=paste0('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/PreProcessed_Data/Methylation/',
    paste0(CancerSite,'_Data'),'.RData'))
  




#stopCluster(cl)





#METdirectories<-list()
#METdirectories$METdirectory27k<-"/home/giordano/Methylation/gdac_20160128/"
#METdirectories$METdirectory450k<-"/home/giordano/Methylation/gdac_20160128/"


#rm(list = ls())

#.rs.restartR()































