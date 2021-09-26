source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')







########## illumina preprocess #################################











Methylation_data_make<-function(CancerSite,directories,lg=F,dtype='Meth_450k',
                               mvt=0.2,bct='limma',normal=T){
  
  
  dtype<-str_split(dtype,'_')
  if(dtype[[1]][1]=='Meth'){
    data.category = "DNA Methylation"
    data.type = "Methylation Beta Value"
    experimental.strategy = "Methylation Array"
    workflow.type = 'Liftover'
  }
  
  if(dtype[[1]][1]=='Meth'&dtype[[1]][2]=='27k'){
    platform = 'Illumina Human Methylation 27'
    f=getFunction('Preprocess_CancerSite_Methylation27k')
    n<-names(directories)[grepl(pattern = '27',names(directories))]
    dir<-eval(parse(text=paste0('directories$',n)))
    dtn='_data_meth_27k'
  }
  
  if(dtype[[1]][1]=='Meth'&dtype[[1]][2]=='450k'){
    platform = 'Illumina Human Methylation 450'
    f=getFunction('Preprocess_CancerSite_Methylation450k')
    n<-names(directories)[grepl(pattern = '450',names(directories))]
    dir<-eval(parse(text=paste0('directories$',n)))
    dtn='_data_meth_450k'
  }
  
  print('Downloading Tumor Metadata from TCGA')
  projects<-tcga_projects()
  projects_d<-as.character(lapply(str_split(projects,'-'),function(x){x[2]}))
  tumor<-projects[grepl(CancerSite,projects)]
  
  
  tumor_query.har<-GDCquery(project = tumor,
                            legacy = lg,
                            data.category = data.category,
                            data.type = data.type,
                            experimental.strategy = experimental.strategy,
                            workflow.type = workflow.type,
                            platform = platform)
  
  
  metainfo=tumor_query.har
  
  if(dtype[[1]][1]=='Meth'){
    print('Removing NAs, SNPs and sex chromossomens from methylation data')
  }
  
  ProcessedData=f(data=NULL,
                  CancerSite,directory = dir,metadata=tumor_query.har, 
                  MissingValueThreshold=mvt, na=TRUE, batchcor=F, batchtp=bct,
                  BatchData=NULL, sex=F,SNP=T)
  
  if(normal){
    
    MET_Data_Cancer<-ProcessedData$MET_Data_Cancer
    MET_Data_Cancer<-MET_Data_Cancer[,order(colnames(MET_Data_Cancer))]  
    if(length(ProcessedData$MET_Data_Normal)>0){
      MET_Data_Normal<-ProcessedData$MET_Data_Normal
      MET_Data_Normal<-MET_Data_Normal[,order(colnames(MET_Data_Normal))]
    }else{
      MET_Data_Normal<-NULL
      print(paste0('No Control data wwas found for ',tumor,'  ',platform,
                   '. Using tumor data only!'))
    }
    
    
    
  }else{
    MET_Data_Cancer<-ProcessedData$MET_Data_Cancer
    MET_Data_Cancer<-MET_Data_Cancer[,order(colnames(MET_Data_Cancer))]  
  }
  
  rm(ProcessedData)
  gc()
  
  batch_info<-function(x){
    a<-x[,2]
    names(a)<-x[,1]
    return(list(a))
  }
  
  

  
  
  desc_data<-function(CancerSite,data,metainfo){
    print('Measuring DSC for all batch types')
    batch_00<-batch_detect(data=data,metainfo=metainfo)
    
    
    
    
    
    
    batch_lst<-list()
    batch_lst[length(batch_lst)+1]<-list('Origin')
    batch_lst[[length(batch_lst)]][2]<-batch_info(batch_00)
    
    bsel<-batch_00[batch_00$DSC>=0.3,]
    bsel<-bsel[order(bsel$DSC,decreasing = T),]
    
    
    #2.2.3 Density plot
    print('Creating density plot for original data')
    plot_ls<-list()
    
    den(data = data,main=paste0(
      "Methylation Data Density Distribution (",CancerSite,")"))
    dist0 <- recordPlot()
    plot_ls[length(plot_ls)+1]<-'dist0'
    plot_ls[[length(plot_ls)]][2]<-list(dist0)
    
    
    # 2.2.3.2 Quantile normalization
    print('Creating density plot for normalized data')
    N_Data<-normalizeBetweenArrays(data,method = 'quantile')
    den(data=data,main=paste0(
      "Methylation Data Density Distribution (",CancerSite,") Quantile Norm"))
    distQ <- recordPlot()
    plot_ls[[length(plot_ls)+1]]<-'distQ'
    plot_ls[[length(plot_ls)]][2]<-list(distQ)
    graphics.off()
    
    
    return(list('Normalized_data'=N_Data,'Batch_list'=batch_lst,
                'Plot_list'=plot_ls,'bsel'=bsel))
    rm(list = c('batch_00','dist0','distQ','bsel','N_Data','batch_lst',
                'plot_ls'))
    gc()
    
  } 
  
  if(normal){
    a=desc_data(CancerSite,data=MET_Data_Cancer,metainfo = tumor_query.har)
    N_Cancer_Data<-a$Normalized_data
    batch_Cancer_lst=a$Batch_list
    plot_Cancer_lst=a$Plot_list
    bsel_Cancer=a$bsel
    
    if(!is.null(MET_Data_Normal)){
      a=desc_data(CancerSite,MET_Data_Normal,metainfo = tumor_query.har)
      N_Normal_Data<-a$Normalized_data
      batch_Normal_lst=a$Batch_list
      plot_Normal_lst=a$Plot_list
      bsel_Normal=a$bsel
    }else{
      N_Normal_Data<-NULL
      batch_Normal_lst=NULL
      plot_Normal_lst=NULL
      bsel_Normal=NULL
      #print(paste0('No Control data wwas found for ',tumor,'  ',platform,
       #            '. Using tumor data only!'))
    }
    
    
    
    
  }else{
    a=desc_data(CancerSite,MET_Data_Cancer,metainfo = tumor_query.har)
    N_Cancer_Data<-a$Normalized_data
    batch_Cancer_lst=a$Batch_list
    plot_Cancer_lst=a$Plot_list
    bsel_Cancer=a$bsel
  }
  
  
  
  
  # 2.2.4 PCA plot
  
  
  
  ppcess_data<-function(x,vname='_Cancer',bsel,metainfo,dir,mvt,bct,
                        plot_ls=plot_Cancer_lst,batch_lst=batch_Cancer_lst){
    print('Starting batch correction')
    
    i=0
    
    while (length(bsel$Covariate)>0 & i<=10) {
      
      i=i+1
      print(paste0('Batch correction iteration ',i))
      if(i==1){
        
        batch=bsel$Covariate[i]
        
        
        a=PCA(x,batch=batch,main=paste0(
          'PCA plot Origin (',batch,')'),meta = metainfo,outrm = F)
        assign(paste0('PCA_Ori_',CancerSite),a)
        plot_ls[[length(plot_ls)+1]]<-paste0('PCA_Ori_',CancerSite)
        plot_ls[[length(plot_ls)]][2]<-list(get(paste0('PCA_Ori_',CancerSite)))
        rm(list = c('a',paste0('PCA_Ori_',CancerSite)))
        gc()
        
        #removing outliers with PCA
        PCAtmp<-PCA(x,batch=batch,main=paste0(
          'PCA plot Non-Corrected (',batch,')'),meta = metainfo,outrm = T)
        if(PCAtmp$Outlier_rm>0){
          x<-PCAtmp$Data
        }
        rm(PCAtmp)
        gc()
        
        assign(paste0('BT_0',1),
               batch_detect(data = x, metainfo = metainfo))
        batch_lst[[length(batch_lst)+1]]<-paste0('BT_0',1)
        batch_lst[[length(batch_lst)]][2]<-batch_info(get(paste0('BT_0',1)))
        btmp<-get(paste0('BT_0',1))
        rm(list = c(paste0('BT_0',1)))
        gc()
        
        a=PCA(x,batch=batch,main=paste0(
          'PCA plot Non-Corrected (',batch,')'),meta = metainfo,outrm = F)
        assign(paste0('PCA_NC_',CancerSite),a)
        plot_ls[[length(plot_ls)+1]]<-paste0('PCA_NC_',CancerSite)
        plot_ls[[length(plot_ls)]][2]<-list(get(paste0('PCA_NC_',CancerSite)))
        rm(list = c('a',paste0('PCA_NC_',CancerSite)))
        gc()
        
        
        
        # 2.3.0 Removing batch effects
        BatchData<-snames(x = metainfo,colnames(x),batch)
        
        
        bcdt=f(data=x,CancerSite,dir, 
               metadata = metainfo,MissingValueThreshold=mvt, 
               na=F, batchcor=T, batchtp=bct, BatchData=BatchData,
               sex=F,SNP=F)
        
        bcdt_f<-bcdt$GEN_Data_Corrected
        
        #2.3.0 Checking batch effects removal
        
        assign(paste0('BT_',batch,'_0',i),
               batch_detect(data = bcdt_f, metainfo = metainfo))
        batch_lst[[length(batch_lst)+1]]<-paste0('BT_',batch,'_0',i)
        batch_lst[[length(batch_lst)]][2]<-batch_info(
          get(paste0('BT_',batch,'_0',i)))
        btmp<-get(paste0('BT_',batch,'_0',i))
        rm(list = c(paste0('BT_',batch,'_0',i)))
        gc()
        
        a=PCA(bcdt_f,batch=batch,main=paste0(
          'PCA plot Corrected (',batch,'_0',i,')'),meta = metainfo,outrm = F)
        assign(paste0('PCA_Co_',batch,'_0',i,'_',CancerSite),a)
        plot_ls[[length(plot_ls)+1]]<-paste0(
          'PCA_Co_',batch,'_0',i,'_',CancerSite)
        plot_ls[[length(plot_ls)]][2]<-list(get(paste0(
          'PCA_Co_',batch,'_0',i,'_',CancerSite)))
        rm(list = c('a',paste0('PCA_Co_',batch,'_0',i,'_',CancerSite)))
        gc()
        
        bsel<-btmp[btmp$DSC>=0.4,]
        bsel<-bsel[order(bsel$DSC,decreasing = T),]
        if(length(colnames(bcdt_f))<10){
          bsel<-data.frame('Covariate'=c(),'DSC'=c())
        }
        
        
      }
      
      
      if(i>1){
        batch=bsel$Covariate[1]
        
        # 2.3.0 Removing batch effects
        BatchData<-snames(metainfo,colnames(bcdt_f),batch)
        
        
        bcdt=f(data=bcdt_f,CancerSite,dir, 
               metadata = metainfo,MissingValueThreshold=mvt, na=F, batchcor=T,
               batchtp=bct, BatchData=BatchData, sex=F,SNP=F)
        
        bcdt_f<-bcdt$GEN_Data_Corrected
        
        #2.3.0 Checking batch effects removal
        
        assign(paste0('BT_',batch,'_0',i),
               batch_detect(data = bcdt_f, metainfo = metainfo))
        batch_lst[[length(batch_lst)+1]]<-paste0('BT_',batch,'_0',i)
        batch_lst[[length(batch_lst)]][2]<-batch_info(get(
          paste0('BT_',batch,'_0',i)))
        btmp<-get(paste0('BT_',batch,'_0',i))
        rm(list = c(paste0('BT_',batch,'_0',i)))
        gc()
        
        a=PCA(bcdt_f,batch=batch,main=paste0(
          'PCA plot Corrected (',batch,'_0',i,')'),meta = metainfo,outrm = F)
        assign(paste0('PCA_Co_',batch,'_0',i,'_',CancerSite),a)
        plot_ls[[length(plot_ls)+1]]<-paste0(
          'PCA_Co_',batch,'_0',i,'_',CancerSite)
        plot_ls[[length(plot_ls)]][2]<-list(get(paste0(
          'PCA_Co_',batch,'_0',i,'_',CancerSite)))
        rm(list=c('a',paste0('PCA_Co_',batch,'_0',i,'_',CancerSite)))
        gc()
        bsel<-btmp[btmp$DSC>=0.4,]
      }
    }
    
    if(length(bsel_Normal$Covariate)==0&is.null(bcdt_f)){
      bcdt_f<-x
    }
    
    
    assign(paste0(CancerSite,vname,'_data'),bcdt_f)
    
    batch_dt<-function(x){
      
      dt<-data.frame('Iteration'=unlist(lapply(x, function(x){x[1]})),
                     'batch'=unlist(lapply(x, function(x){x[[2]][1]})),
                     'center'=unlist(lapply(x, function(x){x[[2]][2]})),
                     'gender'=unlist(lapply(x, function(x){x[[2]][3]})),
                     'TSS'=unlist(lapply(x, function(x){x[[2]][4]})))
      
      
      return(dt)
    }
    
    
    
    bdt<-batch_dt(batch_lst)
    
    plots<-function(x){
      
      p<-lapply(plot_ls, function(x){x[[2]]})
      pn<-lapply(plot_ls, function(x){x[[1]]})
      names(p)<-unlist(pn)
      
      return(p)
      rm(list = c('pn','p'))
      gc()
    }
    p<-plots(plot_ls)
    
    a<-list()
    a[[1]]<-get(paste0(CancerSite,vname,'_data'))
    a[[2]]<-bdt
    a[[3]]<-p
    
    names(a)=c(paste0(CancerSite,vname,dtn),'batch_info','plot_list')
    assign(paste0(CancerSite,vname,dtn),a)
    #keep(list = c(paste0(CancerSite,dtn)),sure = T)
    
    
    return(get(paste0(CancerSite,vname,dtn)))
    
    rm(list = c('a','bdt','p',paste0(CancerSite,vname,'_data'),
                paste0(CancerSite,vname,dtn),'bcdt','bcdt_f','batch_lst',
                'plot_ls','bsel','btmp','metainfo','x','BatchData'))
    gc()
  }
  
  if(normal){
    assign(paste0(CancerSite,'_Cancer',dtn),ppcess_data(x=N_Cancer_Data,
    vname='_Cancer',bsel = bsel_Cancer,metainfo = tumor_query.har, mvt = mvt,
    dir = dir,bct = bct, plot_ls=plot_Cancer_lst, batch_lst=batch_Cancer_lst))
    n1<-paste0(CancerSite,'_Cancer',dtn)
    
    if(!is.null(MET_Data_Normal)){
      assign(paste0(CancerSite,'_Normal',dtn),ppcess_data(x=N_Normal_Data,
      vname='_Normal', bsel = bsel_Normal,metainfo = tumor_query.har, mvt = mvt, 
      dir = dir, bct = bct, plot_ls=plot_Normal_lst, batch_lst=batch_Normal_lst))
    }else{
      a<-list()
      a[[1]]<-NA
      a[[2]]<-NA
      a[[3]]<-NA
      names(a)=c(paste0(CancerSite,'_Normal',dtn),'batch_info','plot_list')
      assign(paste0(CancerSite,'_Normal',dtn),a)
    }
    n2<-paste0(CancerSite,'_Normal',dtn)
    
    return(list(n1=get(
      paste0(CancerSite,'_Cancer',dtn)),
      n2=get(
        paste0(CancerSite,'_Normal',dtn))
      
      ))
    
    rm(list = c('n1','n2',paste0(CancerSite,'_Cancer',dtn),
                paste0(CancerSite,'_Normal',dtn)))
    gc()
  }else{
    n1<-paste0(CancerSite,'_Cancer',dtn)
    assign(n1,ppcess_data(x=N_Cancer_Data,
    vname='_Cancer',bsel = bsel_Cancer,metainfo = tumor_query.har, mvt = mvt,
   dir = dir,bct = bct, plot_ls=plot_Cancer_lst, batch_lst=batch_Cancer_lst))
    
    return(list(n1=get(
      paste0(CancerSite,'_Cancer',dtn))))
    
    rm(list = c('n1',paste0(CancerSite,'_Cancer',dtn)))
    gc()
  }
  
  
}




Preprocess_DNAmethylation <- function(CancerSite, directories,lg=F, mvt = 0.2,
                                      bct='limma',normal=T,svtmp='C:/Users/giordano/') {    
  save(list = c('CancerSite','directories','lg','mvt','bct','normal','svtmp'),
       file = paste0(svtmp,'PP_params.RData'))
  cat("\tProcessing data for",CancerSite,"\n")
  
  
  if (!is.na(directories$METdirectory27k)) {
    cat('\tLoading data for 27k.\n')
    cm=paste0('Rscript ', svtmp,'Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/PP_27k.R')
    system(cm)
    print('Finished Illumina 27k Preprocessing')
  }else{
    print(paste0('No Illumina Methylation 27k data',' was found for ',CancerSite,'!'))
  }
  
  if (!is.na(directories$METdirectory450k)) {
    cat('\tLoading data for 450k.\n')
    cm=paste0('Rscript ', svtmp,'Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/PP_450k.R')
    system(cm)
    print('Finished Illumina 450k Preprocessing')
  }else{
    print(paste0('No Illumina Methylation 450k data',' was found for ',CancerSite,'!'))
  }
  
  if(!is.na(directories$METdirectory27k)){
    ld=paste0(svtmp,'ProcessedData27k.Rdata')
    load(ld)
  }else{
    ProcessedData27k<-c()
  }
  if(!is.na(directories$METdirectory450k)){
    ld=paste0(svtmp,'ProcessedData450k.Rdata')
    load(ld)
  }else{
    ProcessedData450k<-c()
  }
  if(is.na(directories$METdirectory27k)&is.na(directories$METdirectory450k)){
    stop(message('No data could be found!!'))
  }
  
  # check if we want to combine 27k and 450k
  if (length(ProcessedData27k)!=0 & length(ProcessedData450k)!=0 && CancerSite !="LAML") {
    # only do it when enough samples are in both
    if (ncol(ProcessedData27k$n1[[1]])>20 & ncol(ProcessedData450k$n1[[1]])>20) {               
      Mode='450kon27k'               
      cat("\tCombining 450k and 27k by mapping 450k probes to 27k array.\n")
      # Check if there are duplicate samples, remove the 27k ones. 
      OverlapSamplesCancer=intersect(colnames(ProcessedData27k$n1[[1]]),colnames(ProcessedData450k$n1[[1]]))
      if (length(OverlapSamplesCancer)>0) {
        cat("\tCancer sample overlap is not empty: ",length(OverlapSamplesCancer))
        ProcessedData27k$n1$MET_Data_Cancer=ProcessedData27k$n1[[1]][,-OverlapSamplesCancer,drop=FALSE]
      }
      OverlapSamplesNormal=intersect(colnames(ProcessedData27k$n2[[1]]),colnames(ProcessedData450k$n2[[1]]))
      if (length(OverlapSamplesNormal)>0) {
        cat("\tNormal sample overlap is not empty: ",length(OverlapSamplesNormal))
        ProcessedData27k$n2$MET_Data_Normal=ProcessedData27k$MET_Data_Normal[,-OverlapSamplesNormal,drop=FALSE]
      }
      
      # Overlap the probes
      ProcessedData=list(MET_Data_Cancer=c(),MET_Data_Normal=c())
      
      OverlapProbesCancer=intersect(rownames(ProcessedData27k$n1[[1]]),rownames(ProcessedData450k$n1[[1]]))               
      
      ProcessedData$MET_Data_Cancer=cbind(ProcessedData27k$n1[[1]][OverlapProbesCancer,,drop=FALSE],ProcessedData450k$n1[[1]][OverlapProbesCancer,,drop=FALSE])
      gc()
      if ( length(colnames(ProcessedData27k$n2[[1]]))>0 & length(colnames(ProcessedData450k$n2[[1]]))>0 ) {
        OverlapProbesNormal=intersect(rownames(ProcessedData27k$n2[[1]]),rownames(ProcessedData450k$n2[[1]]))               
        ProcessedData$MET_Data_Normal=cbind(ProcessedData27k$n2[[1]][OverlapProbesNormal,,drop=FALSE],ProcessedData450k$n2[[1]][OverlapProbesNormal,,drop=FALSE])
      } else if ( length(colnames(ProcessedData27k$n2[[1]]))>0 ) {
        ProcessedData$MET_Data_Normal=ProcessedData27k$n2[[1]]                    
      } else if ( length(colnames(ProcessedData450k$n2[[1]]))>0 ) {
        ProcessedData$MET_Data_Normal=ProcessedData450k$n2[[1]] 
      } 
      
      
      gc()
    } else if (ncol(ProcessedData27k$n1[[1]])>ncol(ProcessedData450k$n1[[1]])) { 
      cat("\tNot enough 450k samples, only using the 27k (need min 20 samples).\n")
      Mode='27k'
      ProcessedData=ProcessedData27k$n1[[1]]          
    } else {
      cat("\tNot enough 27k samples, only using the 450k (need min 50 samples).\n")
      Mode='450k'
      ProcessedData=ProcessedData450k$n1[[1]]     
    }
  } else if (CancerSite == "LAML") {
    cat("\tLAML is a special case, only using 450k data.\n")
    OverlapSamplesCancer=intersect(colnames(ProcessedData27k$n1[[1]]),colnames(ProcessedData450k$n1[[1]]))
    cat("\tOverlap length is:",length(OverlapSamplesCancer),".\n")
    Mode='450k'          
    ProcessedData=ProcessedData450k          
  } else if (length(ProcessedData27k)!=0) {
    cat("\tOnly 27k samples.\n")
    Mode='27k'
    ProcessedData=ProcessedData27k$n1[[1]]         
  } else {       
    cat("\tOnly 450k samples.\n")
    Mode='450k'
    ProcessedData=ProcessedData450k$n1[[1]]           
  } 
  
  PreProcessed_Data=list('ProcessedData27k'=ProcessedData27k,
                         'ProcessedData450k'=ProcessedData450k,
                         'ProcessedData'= ProcessedData)
  return('PreProcessed_Data'=PreProcessed_Data)
  rm(list=c('PreProcessed_Data','ProcessedData27k','ProcessedData450k',
            'ProcessedData','OverlapSamplesCancer','OverlapProbesCancer',
            'OverlapProbesNormal','OverlapProbesNormal','mode'))
  gc()
}
