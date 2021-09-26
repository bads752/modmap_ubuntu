source('/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Script/aux_scripts/source_lib.R')



####################### Combat Functions #################






#dat=GEN_Data
#saminfo=BatchDataSelected
ComBat_NoFiles <- function(dat, saminfo, type='txt', write=F, covariates='all', par.prior=F, filter=F, skip=0, prior.plots=T){
  #debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
  
  # 'expression_xls' is the expression index file (e.g. outputted by dChip). I think it was replaced by dat, a matrix
  # 'sample_info_file' is a tab-delimited text file containing the colums: Array  name, sample name, Batch, and any other covariates to be included in the modeling. Also I think it was replaced by the data as R objects.
  
  cat('Reading Sample Information File\n')
  #saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
  if(sum(colnames(saminfo)=="Name")!=1){return('ERROR: Sample Information File does not have a Batch column!')}
  
  cat('Reading Expression Data File\n')
  #      if(type=='csv'){
  #           dat <- read.csv(expression_xls,header=T,row.names=1,as.is=T)
  #           #print(dat[1:2,])
  #           #	dat <- dat[,trim.dat(dat)]  
  #           #print(colnames(dat))
  #           #colnames(dat)=scan(expression_xls,what='character',nlines=1,sep=',',quiet=T)[1:ncol(dat)]
  #           #print(colnames(dat))
  #      }  else {
  #dat <- read.table(expression_xls,header=T,comment.char='',fill=T,sep='\t', as.is=T)
  dat <- dat[,trim.dat(dat)]
  #colnames(dat)=scan(expression_xls,what='character',nlines=1,sep='\t',quiet=T)[1:ncol(dat)]
  #      }
  
  
  if (skip>0){
    geneinfo <- as.matrix(dat[,1:skip])
    dat <- dat[,-c(1:skip)]
  } else {
    geneinfo=NULL
  }
  
  if(filter){
    ngenes <- nrow(dat)
    col <- ncol(dat)/2
    present <- apply(dat, 1, filter.absent, filter)
    dat <- dat[present, -(2*(1:col))]
    if (skip>0){geneinfo <- geneinfo[present,]}
    cat('Filtered genes absent in more than',filter,'of samples. Genes remaining:',nrow(dat),'; Genes filtered:',ngenes-nrow(dat),'\n')
  }
  
  if(any(apply(dat,2,mode)!='numeric')){return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
  
  tmp <- match(colnames(dat),saminfo[,1])
  if(any(is.na(tmp))){return('ERROR: Sample Information File and Data Array Names are not the same!')}
  tmp1 <- match(saminfo[,1],colnames(dat))
  saminfo <- saminfo[tmp1[!is.na(tmp1)],]		
  
  if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
  design <- design.mat(saminfo)	
  
  
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  ## Check for missing values
  NAs = any(is.na(dat))
  if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  #print(dat[1:2,])
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}
  
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
  
  ##Get regression batch effect parameters
  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[,1:n.batch]
  if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
  }
  
  ##Find Priors
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  
  
  ##Plot empirical and parametric priors
  
  if (prior.plots & par.prior){
    pdf(file='prior_plots.pdf')
    par(mfrow=c(2,2))
    tmp <- density(gamma.hat[1,])
    plot(tmp,  type='l', main="Density Plot")
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
    qqnorm(gamma.hat[1,])	
    qqline(gamma.hat[1,], col=2)	
    
    tmp <- density(delta.hat[1,])
    invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
    tmp1 <- density(invgam)
    plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
    lines(tmp1, col=2)
    qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
    lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
    title('Q-Q Plot')
    dev.off()
  }
  
  ##Find EB batch adjustments
  
  gamma.star <- delta.star <- NULL
  if(par.prior){
    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch){
      temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }else{
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch){
      temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }
  
  
  ### Normalize the Data ###
  cat("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  if(write) {
    # output_file <- paste(expression_xls,'Adjusted','.txt',sep='_')
    output_file <- 'Adjusted.txt'
    #print(geneinfo[1:2])
    #print(bayesdata[1:2,1:4])
    #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
    #suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
    outdata <- cbind(ProbeID=rownames(dat), bayesdata)
    write.table(outdata, file=output_file, sep="\t")
    cat("Adjusted data saved in file:",output_file,"\n")
  } else {
    return(cbind(rownames(dat),bayesdata))
  }
  
}

# filters data based on presence/absence call
filter.absent <- function(x,pct){
  present <- T
  col <- length(x)/2
  pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
  if(pct.absent > pct){present <- F}
  present
}

# Next two functions make the design matrix (X) from the sample info file 
build.design <- function(vec, des=NULL, start=2){
  tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
  for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
  cbind(des,tmp)
}

design.mat <- function(saminfo){
  tmp <- which(colnames(saminfo) == 'Name')
  tmp1 <- as.factor(saminfo[,tmp])
  cat("Found",nlevels(tmp1),'batches\n')
  design <- build.design(tmp1,start=1)
  ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
  cat("Found",ncov,'covariate(s)\n')
  if(ncov>0){
    for (j in 1:ncov){
      tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
      design <- build.design(tmp1,des=design)
    }
  }
  design
}

# Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(saminfo){
  tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Name')])
  batches <- NULL
  for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
  batches
}

# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
  tmp <- strsplit(colnames(dat),'\\.')
  tr <- NULL
  for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
  tr
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments
it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration function to find the nonparametric adjustments
#sdat=as.matrix(s.data[,batches[[i]]])
#g.hat=gamma.hat[i,]
#d.hat=delta.hat[i,]

int.eprior <- function(sdat,g.hat,d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    print(paste0('### ',i,' ###'))
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x),length(g),n,byrow=T)
    resid2 <- (dat-g)^2
    sum2 <- resid2%*%j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star,sum(g*LH)/sum(LH))
    d.star <- c(d.star,sum(d*LH)/sum(LH))
    #if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust	
} 

#fits the L/S model in the presence of missing data values
Beta.NA = function(y,X){
  des=X[!is.na(y),]
  y1=y[!is.na(y)]
  B <- solve(t(des)%*%des)%*%t(des)%*%y1
  B
}






####################### Dowwnload and Preprocessing functions #################

TCGA_GENERIC_LoadIlluminaMethylationData <- function(Filename) {
  
  # read in an illumina methylation file with the following format: 
  # header row with sample labels
  # 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate
  # The first column has the probe names. 
  MET_Data<-data.table::fread(Filename)
  MET_Data=as.matrix(MET_Data)
  Probes=MET_Data[,1]
  rownames(MET_Data)=Probes
  MET_Data=MET_Data[,-1]
  MET_Data=MET_Data[-1,]
  MET_Data=MET_Data[,seq(1,ncol(MET_Data),4)]     
  class(MET_Data)='numeric'
  
  return(MET_Data)
}



TCGA_GENERIC_LoadIlluminaRNAseqData <- function(Filename,metadata) {
  
  # read in an illumina methylation file with the following format: 
  # header row with sample labels
  # 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate
  # The first column has the probe names. 
  RNA_Data=data.table::fread(Filename)
  RNA_Data=as.matrix(RNA_Data)
  Probes=RNA_Data[,1]
  Probes=str_split(Probes,pattern = '\\|')
  Probes=unlist(lapply(Probes, function(x){
    x[[1]][1]
  }))
  id=Probes=='?'
  RNA_Data=RNA_Data[!id,]
  Probes=Probes[!id]
  
  rownames(RNA_Data)=Probes
  RNA_Data=RNA_Data[,-1]
  class(RNA_Data)='numeric'
  
  x<-metadata
  n<-(x$cases)
  n<-n[order(n)]
  cn<-colnames(RNA_Data)
  
  n<-lapply(cn, function(x){
    n[grepl(x,n)]
  })
  trm<-c()
  for(i  in seq_along(n)){
    if(length(n[[i]])==0|length(n[[i]])>1){
      trm<-append(trm,i)
    }
    
  }
  
  if(!is.null(trm)|length(trm)>0){
    n<-n[-trm]
  }
  
  
  n<-(sort(unlist(n)))
  grepl(paste0(cn,collapse = '|'),n)
  
  tkp<-c()
  for (i in seq_along(cn)) {
    
    a<-grepl(cn[i],n)
    if(sum(a)==1){
      a<-TRUE
    }else{
      a<-FALSE
    }
    
    
    tkp<-append(tkp,a)
    
    
  }
  RNA_Data<-RNA_Data[,tkp]
  cn<-colnames(RNA_Data)
  n<-lapply(cn, function(x){
    n[grepl(x,n)]
  })
  n<-sort(unlist(n))
  
  colnames(RNA_Data)<-n
  
  
  
  return(RNA_Data)
}


TCGA_GENERIC_LoadIlluminamiRseqData <- function(Filename,metadata) {
  
  # read in an illumina methylation file with the following format: 
  # header row with sample labels
  # 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate
  # The first column has the probe names. 
  miRseq_Data=data.table::fread(Filename)
  miRseq_Data=as.matrix(miRseq_Data)
  Probes=miRseq_Data[,1]
  
  
  rownames(miRseq_Data)=Probes
  miRseq_Data=miRseq_Data[,-1]
  class(miRseq_Data)='numeric'
  
  x<-metadata
  n<-(x$cases)
  n<-n[order(n)]
  cn<-colnames(miRseq_Data)
  
  n<-lapply(cn, function(x){
    n[grepl(x,n)]
  })
  trm<-c()
  for(i  in seq_along(n)){
    if(length(n[[i]])==0|length(n[[i]])>1){
      trm<-append(trm,i)
    }
    
  }
  
  if(!is.null(trm)|length(trm)>0){
    n<-n[-trm]
  }
  
  
  n<-(sort(unlist(n)))
  grepl(paste0(cn,collapse = '|'),n)
  
  tkp<-c()
  for (i in seq_along(cn)) {
    
    a<-grepl(cn[i],n)
    if(sum(a)==1){
      a<-TRUE
    }else{
      a<-FALSE
    }
    
    
    tkp<-append(tkp,a)
    
    
  }
  miRseq_Data<-miRseq_Data[,tkp]
  cn<-colnames(miRseq_Data)
  n<-lapply(cn, function(x){
    n[grepl(x,n)]
  })
  n<-sort(unlist(n))
  
  colnames(miRseq_Data)<-n
  
  
  
  return(miRseq_Data)
}


#' The TCGA_GENERIC_GetSampleGroups function
#' 
#' Internal. Looks for the group of the samples (normal/cancer).
#' @param SampleNames vector with sample names.
#' @return a list.
#' @keywords internal
#'
#SampleNames<-x$cases
TCGA_GENERIC_GetSampleGroups <-function(SampleNames) {
  
  # First replace any . with - so the sample groups are uniform. 
  SampleGroups=list()
  
  #1: Primary Tumor
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]01[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Primary=SampleNames[Matches==1]
  
  #2: Recurrent tumor
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]02[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Recurrent=SampleNames[Matches==1]
  
  #3: Primary blood derived cancer - peripheral blood
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]03[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$PeripheralBloodCancer=SampleNames[Matches==1]     
  
  #10: Blood derived normal
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]10[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$BloodNormal=SampleNames[Matches==1]     
  
  #11: Solid tissue derived normal
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]11[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$SolidNormal=SampleNames[Matches==1]     
  
  #20 Cellines
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]20[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$CellLines=SampleNames[Matches==1]    
  
  #06 Cellines
  Matches=regexpr("TCGA[.|-]\\w\\w[.|-]\\w\\w\\w\\w[.|-]06[.|-]*",SampleNames,perl=FALSE,useBytes=FALSE)
  SampleGroups$Metastatic=SampleNames[Matches==1]         
  
  return(SampleGroups)
}

#' The TCGA_GENERIC_CleanUpSampleNames function
#' 
#' Internal. Cleans the samples IDs into the 12 digit format and removes doubles.
#' @param GEN_Data data matrix.
#' @param IDlength length of samples ID.
#' @return data matrix with cleaned sample names.
#' @keywords internal
#'


TCGA_GENERIC_CleanUpSampleNames <-function(GEN_Data, IDlength = 12) {     
  SampleNames=colnames(GEN_Data)
  SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
  if (length(SampleNamesShort)!=length(unique(SampleNamesShort))) {
    # remove the doubles           
    Counts=table(SampleNamesShort)
    Doubles=rownames(Counts)[which(Counts>1)]
    
    cat("Removing doubles for",length(Doubles),"samples.\n")
    for(i in 1:length(Doubles)) {                         
      CurrentDouble=Doubles[i]          
      pos=grep(CurrentDouble,SampleNames)
      #GEN_Data[1:10,pos]
      #cor(GEN_Data[,pos])
      GEN_Data=GEN_Data[,-pos[2:length(pos)]]     
      SampleNames=colnames(GEN_Data) # need to update samplenames because pos is relative to this
    }
    SampleNames=colnames(GEN_Data)
    SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
    
    # now set the samplenames
    colnames(GEN_Data)=SampleNamesShort
  } else {
    colnames(GEN_Data)=SampleNamesShort     
  }     
  return(GEN_Data)
}

Download_TCGA_data <- function(CancerSite, 
                               TargetDirectory, 
                               downloadData = TRUE,
                               data='RNAseq') {    
  
  dir.create(TargetDirectory,showWarnings=FALSE)
  
  # download the 27k data
  dataType='stddata'
  
  if(data=='RNAseq'){
    
    dataFileTag='mRNAseq_Preprocess.Level_3'
    cat('Searching mRNAseq data for:',CancerSite,'\n')
    directories=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    return(directories=list(directories))
  }
  
  if(data=='Methylation'){
    dataFileTag='Merge_methylation__humanmethylation27'
    cat('Searching 27k MET data for:',CancerSite,'\n')
    METdirectory27k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    
    # download the 450k data
    dataFileTag='Merge_methylation__humanmethylation450'
    cat('Searching 450k MET data for:',CancerSite,'\n')
    METdirectory450k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    return(directories=list(METdirectory27k=METdirectory27k,METdirectory450k=METdirectory450k))
    
  }
  
  if(data=='miRseq'){
    dataFileTag='miRseq_Preprocess.Level_3'
    cat('Searching miRNAseq data for:',CancerSite,'\n')
    directories=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    return(directories=list(directories))
  }
  
  
  
  
}



##saveDir="C:/Users/giordano/MethylMix/"
#TCGA_acronym_uppercase = 'READ'
get_firehoseData <- function(downloadData=TRUE,
                             saveDir = "./",
                             TCGA_acronym_uppercase = "BRCA",
                             dataType="stddata",
                             dataFileTag = "mRNAseq_Preprocess.Level_3",
                             FFPE=FALSE,
                             fileType= "tar.gz",
                             gdacURL= "http://gdac.broadinstitute.org/runs/",
                             untarUngzip=TRUE,
                             printDisease_abbr=FALSE){  
  options(timeout=1800)
  # Cases Shipped by BCR  # Cases with Data*  Date Last Updated (mm/dd/yy)
  cancers <- c("Acute Myeloid Leukemia [LAML] \n","Adrenocortical carcinoma [ACC] \n",
               "Bladder Urothelial Carcinoma [BLCA] \n",  "Brain Lower Grade Glioma [LGG] \n",
               "Breast invasive carcinoma [BRCA] \n","Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
               "Cholangiocarcinoma [CHOL] \n",    "Colon adenocarcinoma [COAD] \n",   "Esophageal carcinoma [ESCA] \n",
               "Glioblastoma multiforme [GBM] \n",    "Head and Neck squamous cell carcinoma [HNSC]   \n",
               "Kidney Chromophobe [KICH] \n","Kidney renal clear cell carcinoma [KIRC]   \n",
               "Kidney renal papillary cell carcinoma [KIRP]  \n","Liver hepatocellular carcinoma [LIHC]  \n",
               "Lung adenocarcinoma [LUAD]    \n", "Lung squamous cell carcinoma [LUSC] \n",
               "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]    \n","Mesothelioma [MESO] \n",
               "Ovarian serous cystadenocarcinoma [OV]    \n","Pancreatic adenocarcinoma [PAAD]   \n",
               "Pheochromocytoma and Paraganglioma [PCPG] \n","Prostate adenocarcinoma [PRAD] \n",
               "Rectum adenocarcinoma [READ]  \n","Sarcoma [SARC] \n","Skin Cutaneous Melanoma [SKCM] \n",
               "Stomach adenocarcinoma [STAD] \n","Testicular Germ Cell Tumors [TGCT] \n","Thymoma [THYM] \n",
               "Thyroid carcinoma [THCA]  \n","Uterine Carcinosarcoma [UCS]    \n",
               "Uterine Corpus Endometrial Carcinoma [UCEC]   \n","Uveal Melanoma [UVM] \n", "Colorectal Adenocarcinoma [COADREAD] \n")
  
  if(printDisease_abbr){      
    return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",cancers));      
  }
  gdacURL_orig <- gdacURL
  
  # New code to handle dates - Marcos
  gdacURLnew <- paste0(gdacURL_orig, dataType, "__latest/data/", TCGA_acronym_uppercase, "/")
  urlDataNew <- RCurl::getURL(gdacURLnew)
  urlDataNew <- limma::strsplit2(urlDataNew, "href=\\\"") #regular expressions: need \ to have R recognize any " or \ that's actually in our text
  getlatestdate <- urlDataNew[grep("^201[0-9][01][0-9][0123][0-9]", urlDataNew)]
  getlatestdate <- substring(getlatestdate, 1, 8)
  gdacURLnew <- paste0(gdacURLnew, getlatestdate, "/")
  urlData <- RCurl::getURL(gdacURLnew)
  urlData <- limma::strsplit2(urlData, "href=\\\"") #regular expressions: need \ to have R recognize any " or \ that's actually in our text
  lastDateCompress <- lastDate <- getlatestdate # for compatibility with the rest of old code
  gdacURL <- gdacURLnew # for compatibility with the rest of old code
  # end New code
  
  #remove any FFPE datasets, or only keep those depending on user inputs.
  if (FFPE) { 
    urlData <- urlData[grep("FFPE",urlData)]    
    if(length(urlData)==0){     
      stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")      
    }     
  } else {    
    #we DON'T want FFPE data.
    #but if no FFPE data to begin with: don't subset on this.
    if(length(grep("FFPE",urlData))>0){     
      urlData <- urlData[-grep("FFPE",urlData)]       
    }
    if(length(urlData)==0){     
      stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")       
    }
  }
  #now get full dataset name.
  fileName <- urlData[grep(dataFileTag,urlData)]
  
  if(length(fileName)==0){      
    #warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase," for data type ",dataFileTag ,".No data was downloaded.\n")
    #warning(warnMessage)
    cat("\tThere is no",dataFileTag,"data for",TCGA_acronym_uppercase,"\n")
    return(NA)    
  }
  #some redundancy..but that' OK because we'll add back on the unique tar.gz file tag.
  #first file is one we want - not md5 file.
  fileName <- limma::strsplit2(fileName,"tar.gz")[1,1]
  fileName <- paste(fileName,fileType,sep="")
  
  #final download url
  gdacURL <- paste(gdacURL,fileName,sep="")
  # Directory for downloads
  saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
  
  if (!grepl("Windows", Sys.info()['sysname'])) {
    # Not Windows
    tarfile=paste0(saveDir,fileName)
    finalDir <-  strsplit(tarfile, paste0(".", fileType))[[1]][1]
    if(downloadData){       
      cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")             
      cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
      dir.create(saveDir,showWarnings=FALSE)
      # download file     
      setwd(saveDir)              
      download.file(gdacURL,fileName,quiet=FALSE,mode="wb",method = 'libcurl')
      #this assumes a tar.gz file.
      if(fileType=="tar.gz" && untarUngzip) {                 
        cat("\tUnpacking data.\n")
        tarfile=paste0(saveDir,fileName)
        untar(tarfile)
        #remove tarred file
        fileToRemove <- limma::strsplit2(gdacURL,"/")[ ,ncol(limma::strsplit2(gdacURL,"/"))]
        if(grepl('meth',gdacURL)){
          #remove tarred file
          fileToRemove <- limma::strsplit2(gdacURL,"/")[ ,ncol(limma::strsplit2(gdacURL,"/"))]
          removed <- file.remove(paste0(saveDir,fileToRemove))
        }
        if(grepl('mRNA',gdacURL)){
          f<-list.files(list.dirs(finalDir),full.names = T)
          c<-f[grepl(paste0(CancerSite,'.uncv2.mRNAseq_raw'),f)]
          fr<-f[!grepl(paste0(CancerSite,'.uncv2.mRNAseq_raw'),f)]
          
          removed <- file.remove(fr)
          file.remove(tarfile)
        }
        if(grepl('miRseq',gdacURL)){
          f<-list.files(list.dirs(finalDir),full.names = T)
          c<-f[grepl(paste0('miRseq_raw_counts'),f)]
          fr<-f[!grepl(paste0('miRseq_raw_counts'),f)]
          
          removed <- file.remove(fr)
          file.remove(tarfile)
        }
        
      } else if(untarUngzip) {        
        warning("File expansion/opening only built in for tar.gz files at the moment.\n")       
      }       
      cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
    } else {
      cat("\tdownload data url is :\n ",gdacURL,'\n')
    }
    DownloadedFile=paste0(finalDir,'/')
    return(DownloadedFile)
  } else {
    # new code to handle long names in windows - Marcos
    # WINDOWS: name of file can be too long (and it's repeated in the folder and the file) and windows
    # internally will use another representation for the name of the folder, so then when
    # we want to load the file it says it doesn't exist. So I'm changing the name of the folder
    # to prevent this. We can't change the name of the file as it's used in the Preprocess functions
    idx <- which(sapply(c("methylation27", "methylation450", "mRNAseq", "miRseq"), grepl, dataFileTag))[1]
    if(idx==1){newtag='meth27'}
    if(idx==2){newtag='meth450'}
    if(idx==3){newtag='geneexp'}
    if(idx==4){newtag='miRNA'}
    nameForFolder <- paste(TCGA_acronym_uppercase, dataType, newtag, sep = "_")
    nameForDownloadedFile <- paste0(nameForFolder, ".", fileType)
    nameForDownloadedFileFullPath <- paste0(saveDir, nameForDownloadedFile)
    finalDir <- paste0(saveDir, nameForFolder)
    if (downloadData) {
      cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")             
      cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
      dir.create(saveDir, showWarnings = FALSE)
      
      # Create a virtual drive to overcome long names issue in Windows
      saveDir2 <- gsub("\\\\", "/", saveDir)
      saveDir2 <- substr(saveDir2, 1, nchar(saveDir2) - 1)
      system(paste("subst x:", saveDir2))
      
      download.file(gdacURL, destfile = paste0("x://", nameForDownloadedFile), quiet = FALSE, mode = "wb")
      #this assumes a tar.gz file.
      if(fileType == "tar.gz" && untarUngzip) {                   
        cat("\tUnpacking data.\n")
        untar(nameForDownloadedFileFullPath, exdir = saveDir)
        # untar(paste0("x://", nameForDownloadedFile), exdir = "x://") # doesn't work because it calls system and system doesnt know about the virtual drive
        removed <- file.remove(paste0("x://", nameForDownloadedFile))
        # Anyway I change folder name to make it shorter
        changed <- file.rename(from = paste0("x://", gsub(".tar.gz", "", fileName)), to = paste0("x://", nameForFolder))
        system("subst x: /D") # stop the virtual drive
      } else if(untarUngzip) {        
        warning("File expansion/opening only built in for tar.gz files at the moment.\n")       
      }
      cat("\tFinished downloading", dataFileTag, "data to", finalDir,"\n")
    } else {
      cat("\tdownload data url is :\n ", gdacURL, '\n')
    }
    DownloadedFile = paste0(finalDir, '/')
    return(DownloadedFile)
    # end new code
  }
}


TCGA_Process_EstimateMissingValues <- function(Data, MissingValueThreshold = 0.2,sex=T,SNP=T) {
  
  # FIRST REMOVING BAD PATIENTS
  # removing patients with too many missings values  
  MET_Data=Data
  NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x)))/nrow(MET_Data)
  cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than",MissingValueThreshold*100,"% missing values.\n")
  if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold,drop=FALSE]
  rm(list = c('NrMissingsPerSample'))
  gc()
  
  # removing clones with too many missing values (n>5)
  if(length(colnames(MET_Data))>=5){
    NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x)))/ncol(MET_Data)
    cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than",MissingValueThreshold*100,"% missing values.\n")
    if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,,drop=FALSE]
    rm(list = c('NrMissingsPerGene'))
    gc()
  }else{
    print('Number of samples is too small (N<5) to remove samples with missing probes')
  }
  
  
  
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  data("Locations")
  Locations<-as.data.frame(Locations)
  
  # removing sex chromosomes
  keep <- !(rownames(MET_Data) %in% rownames(Locations)[Locations$chr %in% c("chrX","chrY")])
  rn=rownames(MET_Data)[keep]
  cn=colnames(MET_Data)
  MET_Data<-as.data.frame(MET_Data[keep,])
  rownames(MET_Data)=rn
  colnames(MET_Data)=cn
  # removing SNP probes
  data("SNPprobes")
  SNPprobes<-(SNPprobes)
  keep <- !(rownames(MET_Data) %in% SNPprobes)
  rn=rownames(MET_Data)[keep]
  cn=colnames(MET_Data)
  MET_Data<-as.data.frame(MET_Data[keep,])
  rownames(MET_Data)=rn
  colnames(MET_Data)=cn
  rm(list = c('Locations','keep','SNPprobes'))
  gc()
  
  # knn impute using Tibshirani's method
  if (length(colnames(MET_Data))>1) {
    k=15
    cn=colnames(MET_Data)
    rn=rownames(MET_Data)
    KNNresults=impute::impute.knn(as.matrix(MET_Data),k)
    MET_Data_KNN=KNNresults$data
    colnames(MET_Data_KNN)=cn
    rownames(MET_Data_KNN)=rn
    rm(list = c('KNNresults','MET_Data','k'))
    gc()
    # cleaning up sample names     
    return(MET_Data_KNN)
    
  } else {
    # when only 1 sample,need to make a matrix again
    #MET_Data=as.matrix(MET_Data)
    
    return(MET_Data)    
  }     
}

#tp='sex'
#x=tumor_query.har
#cn=colnames(data)
##x=MET_Data_Cancer
#cn=colnames(MET_Data_Cancer)
snames<-function(x,cn=NULL,tp='sex'){
  if(tp=='batch'|tp=='TSS'|tp=='center'|tp=='sex')  {
    
  }else{
    stop(message('Incorrect batch type'))
  }
  
  if(class(x)=='RangedSummarizedExperiment'){
    assay<-assay(x)
    n<-colnames(assay)[order(colnames(assay))]
    n<-lapply(cn, function(x){
      n[grepl(x,n)]
    })
    n<-sort(unlist(n))
  }else{
    
    x<-x[[1]][[1]]
    n<-(x$cases)
    n<-n[order(n)]
    
    
    #ns<-str_split(n,'-')
    #ns<-lapply(ns, function(x){c(x[1:3],str_sub(x[4],1,2))})
    #ns<-unlist(lapply(ns, function(x){paste0(x,collapse = '-')}))
    

    
    n<-lapply(cn, function(x){
      n[grepl(x,n)]
    })
    n<-(sort(unlist(n)))
  }
  
  #accounting batch for different experiment
  
  
  
  #if(!is.null(cn)){n<-n[n%in%cn]}else{n<-colnames(assay)}
  a=stri_extract_all_regex(str = n, pattern = paste(cn, collapse = "|"))
  a<-unique(a)
  if(!is.null(cn)){
    n<-n
    no<-sort(unlist(a))
  }else{stop()}
  
  
  
  nn<-n
  n<-str_split(n,'-')
  bt<-lapply(n, function(x){
    x[[6]]
  })
  ts<-lapply(n, function(x){
    x[[2]]
  })
  ct<-lapply(n, function(x){
    x[[7]]
  })
  btno<-lapply(unique(bt), function(x){
    x
  })
  tso<-lapply(unique(ts), function(x){
    x
  })
  cto<-lapply(unique(ct), function(x){
    x
  })
  btnf<-rep(0,(length(unique(bt))))
  for (i in seq_along(unique(bt))) {
    btnf[i]<-paste0('B',i)
  }
  dt_name<-as.character(unlist(btno))
  names(dt_name)<-btnf
  rm(list = c('btnf','btno'))
  gc()
  
  tsf<-rep(0,(length(unique(ts))))
  for (i in seq_along(unique(ts))) {
    tsf[i]<-paste0('T',i)
  }
  dt_name2<-as.character(unlist(tso))
  names(dt_name2)<-tsf
  rm(list = c('tsf','tso'))
  gc()
  
  ctf<-rep(0,(length(unique(ct))))
  for (i in seq_along(unique(ct))) {
    ctf[i]<-paste0('C',i)
  }
  dt_name3<-as.character(unlist(cto))
  names(dt_name3)<-ctf
  rm(list = c('ctf','cto'))
  gc()
  
  #accounting for gender induced batch effect
  if(class(x)=='RangedSummarizedExperiment'){
    meta<-as.data.frame(colData(x))
    m<-meta$barcode[meta$gender=='male']
    f<-meta$barcode[meta$gender=='female']
    t<-meta$barcode[meta$sample_type%in%c('Primary Tumor','Primary Blood Derived Cancer - Peripheral Blood')]
    h<-meta$barcode[meta$sample_type%in%c('Solid Tissue Normal')]
  }else{
    barcode=as.character(unlist(nn))
    barcode<-str_split(barcode,'-')
    barcode<-lapply(barcode, function(x){c(x[1:3])})
    barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
    
    meta <- GDCquery_clinic(unique(x$project), type = "clinical")
    meta<-subset(meta,meta$submitter_id%in%barcode)
    if(length(meta$submitter_id[meta$gender=='male'])>0){
      m<-meta$submitter_id[meta$gender=='male']
      m<-unlist(lapply(m, function(x){
        nn[grepl(x,nn)]
      }))
    }else{
      m<-NULL
      print(paste0('No male subjects found for ',x$project[1], ' ', x$platform[1],'!'))
    }
    
    if(length(meta$submitter_id[meta$gender=='female'])>0){
      f<-meta$submitter_id[meta$gender=='female']
      f<-unlist(lapply(f, function(x){
        nn[grepl(x,nn)]
      }))
    }else{
      f<-NULL
      print(paste0('No female subjects found for ',x$project[1], ' ', x$platform[1],'!'))
    }
    
    t<-x[x$sample_type%in%c('Primary Tumor','Primary Blood Derived Cancer - Peripheral Blood'),]
    if(nrow((x[x$sample_type%in%c('Solid Tissue Normal'),]))>0){
      h<-x[x$sample_type%in%c('Solid Tissue Normal'),]
    }else{
      h<-NULL
      print(paste0('No control data found for ',x$project[1], ' ', x$platform[1],'!'))
    }
    
    
    
  }
  
  
  
  n1<-nn
  n2<-rep(0,length(n1))
  c<-rep(0,length(n1))
  for(i in seq_along(n1)){
    
    nts=names(dt_name2)[dt_name2==n[[i]][2]]
    nb<-names(dt_name)[dt_name==n[[i]][6]]
    
    if(!is.null(m)){
      if(n1[i]%in%m){
        ng<-'M'
      }}
    if(!is.null(f)){
          ng<-'F'
        }else{
          ng<-'none'
        }
    
    
    
    if(n1[i]%in%t$cases){
      nt<-'T'
    }else{
      if(!is.null(h)){
        nt<-'H'
      }else{
        nt<-'none'
      }}
    nct<-names(dt_name3)[dt_name3==n[[i]][7]]
    #n2<-append(n2,paste0(nts,nb,ng,nt,nct))
    #n2<-append(n2,paste0(nb))
    c[i]<-no[i]
    if(tp=='TSS'){
      n2[i]<-paste0(nts)
    }
    if(tp=='batch'){
      n2[i]<-paste0(nb)
    }
    if(tp=='sex'){
      n2[i]<-paste0(ng)
    }
    if(tp=='center'){
      n2[i]<-paste0(nct)
    }
    if(tp=='cancer'){
      n2[i]<-paste0(nt)
    }
  }
  nm<-data.frame('Sample'=c,'Name'=n2)
  return(nm)
  rm(list = c('x','meta','nm','m','f','t','h','barcode','n1','n2','nn','c',
              'nts','nb','ng','nt','nct','no','dt_name','dt_name2','dt_name3'))
  gc()
}




#BatchData_mp=BatchData
#BatchData=BatchData[]
##GEN_Data=GEN_Data_Corrected
#BatchData=get('BatchData')

TCGA_GENERIC_CheckBatchEffect <-function(GEN_Data, BatchData) {
  
  # first match the samples to the batch
  Order=match(colnames(GEN_Data),BatchData[,1])
  BatchDataSelected=BatchData[Order,]
  BatchDataSelected$Name <- factor(BatchDataSelected$Name) 
  
  # PCA analysis
  # alternatively use fast.prcomp from package gmodels, but tests do not show this is faster
  PCAanalysis=prcomp(t(GEN_Data))
  PCdata=PCAanalysis$x
  #plot(PCdata[,1]~BatchDataSelected[,3])
  
  if (length(unique(BatchDataSelected$Name[!is.na(BatchDataSelected$Name)]))>1) {
    tmp=aov(PCdata[,1]~BatchDataSelected[,2])          
    return(list(Pvalues=summary(tmp),PCA=PCdata,BatchDataSelected=BatchDataSelected))
  } else {
    return(-1)
  }
}


TCGA_GENERIC_BatchCorrection <-function(GEN_Data,BatchData) {
  
  # select only samples with batch, others get deleted
  WithBatchSamples=is.element(colnames(GEN_Data),BatchData[,1])
  if (length(which(WithBatchSamples==FALSE))>0) GEN_Data=GEN_Data[,-which(WithBatchSamples==FALSE)]
  
  # select only the batch data that is present in the current data set, remove others (remember, the batch data is for all of TCGA)
  PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
  BatchDataSelected=BatchData
  if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
  BatchDataSelected$Name<- factor(BatchDataSelected$Name)
  BatchDataSelected$Sample <- factor(BatchDataSelected$Sample)
  
  # reordening samples (not really necessary as Combat does this too)
  order <- match(colnames(GEN_Data),BatchDataSelected[,1])
  BatchDataSelected=BatchDataSelected[order,]
  BatchDataSelected$Name<- factor(BatchDataSelected$Name)
  
  # running combat
  CombatResults=ComBat_NoFiles(GEN_Data,BatchDataSelected,par.prior = T)
  
  GEN_Data_Corrected=CombatResults[,-1]
  class(GEN_Data_Corrected) <- "numeric"
  return(GEN_Data_Corrected)
}


#GEN_Data=data
TCGA_BatchCorrection_MolecularData <- function (GEN_Data,BatchData,MinInBatch=4,tp='limma') {
  #MinInBatch=MinPerBatch
  # Remove samples with batch number 0
  if (length(-which(BatchData[,2]==0))>0) {
    BatchData=BatchData[-which(BatchData[,2]==0),]
  }
  
  # remove batches that are too small     
  #MinInBatch=5
  PresentSamples=is.element(BatchData[,1],colnames(GEN_Data))
  
  # changes this April 2014, such that the BatchDataSelected only deals with samples in the current GEN_Data
  BatchDataSelected=BatchData[PresentSamples,] 
  if (sum(PresentSamples) != length(colnames(GEN_Data))) BatchDataSelected=BatchData[-which(PresentSamples==FALSE),]
  BatchDataSelected$Name<- factor(BatchDataSelected$Name)
  
  NrPerBatch=table(BatchDataSelected$Name)
  SmallBatches=NrPerBatch<MinInBatch
  BatchesToBeRemoved=names(SmallBatches)[which(SmallBatches==TRUE)]
  SamplesToBeRemoved=as.character(BatchDataSelected[which(BatchDataSelected$Name %in% BatchesToBeRemoved),1])
  
  if (length(colnames(GEN_Data))-length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>5) { # just checking if we have enough samples after removing the too small batches
    if (length(which(colnames(GEN_Data) %in% SamplesToBeRemoved))>0) {
      cat("Removing",length(which(colnames(GEN_Data) %in% SamplesToBeRemoved)),"samples because their batches are too small.\n")
      GEN_Data_rm=GEN_Data[,-which(colnames(GEN_Data) %in% SamplesToBeRemoved)]
    }         
    # batch correction with Combat, incorporate check for only 1 batch
    if(length(SamplesToBeRemoved)==0){
      GEN_Data_rm=GEN_Data
    }
    #BatchCheck=TCGA_GENERIC_CheckBatchEffect(GEN_Data_rm,BatchData)
    
    if (length(unique(BatchData[,2][!BatchData[,1]%in%SamplesToBeRemoved]))<=1) {
      cat("Only one batch after removal, no batch correction possible.\n")
      return(list(GEN_Data_Corrected=GEN_Data_rm,GEN_Data=GEN_Data))
      rm(list = c('GEN_Data_rm','GEN_Data','NrPerBatch','SmallBatches',
                  'BatchesToBeRemoved','SamplesToBeRemoved','BatchData',
                  'MinInBatch','BatchDataSelected','PresentSamples','tp'))
      gc()
    }else{
      
      if (tp == 'ComBat') {
        GEN_Data_Corrected=TCGA_GENERIC_BatchCorrection(
          GEN_Data_rm,BatchData[!BatchData[,1]%in%SamplesToBeRemoved,])
        return(list(GEN_Data_Corrected=GEN_Data_Corrected,GEN_Data=GEN_Data))
        rm(list = c('GEN_Data_Corrected','GEN_Data_rm','GEN_Data','NrPerBatch',
                    'SmallBatches',
                    'BatchesToBeRemoved','SamplesToBeRemoved','BatchData',
                    'MinInBatch','BatchDataSelected','PresentSamples','tp'))
        gc()
      } else{if (tp == 'limma') {
        
        GEN_Data_Corrected=limma::removeBatchEffect(
          GEN_Data_rm,BatchData[,2][!BatchData[,1]%in%SamplesToBeRemoved])
       
        return(list(GEN_Data_Corrected=GEN_Data_Corrected,GEN_Data=GEN_Data))}
        rm(list = c('GEN_Data_Corrected','GEN_Data_rm','GEN_Data','NrPerBatch',
                    'SmallBatches',
                    'BatchesToBeRemoved','SamplesToBeRemoved','BatchData',
                    'MinInBatch','BatchDataSelected','PresentSamples','tp'))
        gc()
    } 
   }
  } else {
    cat("The nr of samples becomes to small, no batch correction possible.\n")
    return(list(GEN_Data_Corrected=GEN_Data,GEN_Data=GEN_Data))
    rm(list = c('GEN_Data_rm','GEN_Data','NrPerBatch',
                'SmallBatches',
                'BatchesToBeRemoved','SamplesToBeRemoved','BatchData',
                'MinInBatch','BatchDataSelected','PresentSamples','tp'))
    gc()
  }
  
}


TCGA_fread_load<-function(data,meth_reduce=T,csv.name,delim="tab",
                          fname, show_progress_bar = T,
                          Probes, sdir, CancerSite, cf, rf,
                          cols, return=T){
  cols<-colnames(data)
  rd<-cols
  if(meth_reduce){
    if(is.null(data)|length(colnames(data))==0){
      stop(message('Error: data argument not found'))
    }
    cols<-colnames(data)
    rd<-cols
    co<-cols[1]
    rd<-rd[-1]
    rd<-rd[seq(1,length(rd),4)]  
    rd<-append(co,rd)
    
    
  }
  
  
  save(cols,file =cf)
  save(rd,file=rf)
  cm<-paste('Rscript /home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Script/aux_scripts/Methylation_fread_script.R', 
            csv.name, delim='tab', fname ,show_progress_bar, Probes, sdir,
            paste0("'",CancerSite,"'"), cf, 
            paste0("'",meth_reduce,"'"), rf)
  system(cm)
  if(return){
    load(paste0(sdir,CancerSite,'_Meth_raw.RData'))
    MET_Data<-(get(paste0(CancerSite,'_Meth_raw')))
    
    MET_Data<-as.data.frame(MET_Data)
    
    
    
    for(i in seq_along(colnames(MET_Data))){
      if(i>1){
        MET_Data[,i]<-as.numeric(MET_Data[,i])
      }  
    }
    MET_Data<-aggregate(MET_Data[, 2:length(colnames(MET_Data))],
                        list(MET_Data$`Hybridization REF`), mean)
    
    a<-as.character(MET_Data[,1])
    MET_Data<-MET_Data[,-1]
    
    rownames(MET_Data)<-a
    if(sum(MET_Data>1,na.rm = T)>0){
      MET_Data<NULL
      stop(message('Error: Methylation beta value cannot be >1!'))
    }else{
      return(MET_Data)
      rm(list = c(paste0(CancerSite,'_raw'),a))
      gc()
    }
    
    
  }
}


### Preprocssing of data, 27 k and 450k separetally


#METdirectory<-METdirectories$METdirectory27k
#Samplegroups<-SampleGroups

Preprocess_CancerSite_Methylation27k <- function(data=NULL,CancerSite, directory, metadata,MissingValueThreshold = 0.2, na=TRUE, batchcor=TRUE, batchtp='limma', BatchData=NULL,sex=F,SNP=T) {
  
  # Settings
  #get("BatchData")
  MinPerBatch=4
  metadata<-metadata[[1]][[1]]
  
  if(is.null(data)) {
    
    if (grepl("Windows", Sys.info()['sysname'])) {
      # If Windows I'll create a virtual drive to handle the long file names issue
      # Create a virtual drive to overcome long names issue in Windows
      METdirectory<-gsub('.{1}$', '', directory)
      virtualDir <- METdirectory
      system(paste("subst x:", virtualDir))
      
      # Load data
      METfiles <- dir("x:")
      
      METfiles = list.files(directory,
                            full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,METfiles)
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('methylation__humanmethylation27',METfiles)  
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('.txt',METfiles)
      METfiles=METfiles[MatchedFilePosition]
      Filename=METfiles    
      s<-file.size(Filename)
      s<-s/1e9
      delim="\t"
      pre_process_size=1000
      df <- read_delim(Filename, delim = delim, n_max = pre_process_size,name_repair = 'minimal')
      nc<-length(colnames(df))
      if(s<1){
        MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
      }else{
        if(nc<=1000){
          cols<-colnames(df)
          csv.name=Filename
          delim="'\t'"
          fname<-paste0(CancerSite,'_Meth','_raw')
          show_progress_bar = T
          Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
          sdir=directory
          CancerSite=CancerSite
          cf<- paste0(directory,'/',CancerSite,'_cols.RData')
          rf<- paste0(directory,'/',CancerSite,'_rdcols.RData')
          MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                    delim = delim, fname = fname, 
                                    show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                    CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
        }
        
        if(nc>1000) {
          if(as.integer(nc/4)+1 <=1000) {
            cols<-colnames(df)
            csv.name=Filename
            delim="'\t'"
            fname<-paste0(CancerSite,'_Meth','_raw')
            show_progress_bar = T
            Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
            sdir=directory
            CancerSite=CancerSite
            cf<- paste0(directory,'/',CancerSite,'_cols.RData')
            rf<- paste0(directory,'/',CancerSite,'_rdcols.RData')
            MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                      delim = delim, fname = fname, 
                                      show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                      CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
          } else {
            stop(message('Number of columns is too large!'))
          }
        }
      }
    } else {
      # Not windows
      # Load data
      METfiles = list.files(directory,
                            full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,METfiles)
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('methylation__humanmethylation27',METfiles)  
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('MANIFEST',METfiles)
      METfiles=METfiles[-MatchedFilePosition]
      MatchedFilePosition=grep('.txt',METfiles)
      METfiles=METfiles[MatchedFilePosition]
      Filename=METfiles    
      df <- read_delim(Filename, delim = '\t', n_max = 100)
      
      #creating new colnames to avoid duplications
      cr<-colnames(df)
      cr<-t(paste0(cr,collapse = '\t'))
      cm<-paste0('sed -i  "1i ',cr,'" ',Filename)
      system(cm)
      df <- read_delim(Filename, delim = '\t', n_max = 100)
      
      MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    }
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    #Samplegroups2=TCGA_GENERIC_GetSampleGroups(x$cases)
    if (CancerSite=='LAML') {
      MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
      rn=rownames(MET_Data_Cancer)
      rn=rn[!duplicated(rn)]
      cn=colnames(MET_Data_Cancer)
      
      MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
      colnames(MET_Data_Cancer)=cn
      rownames(MET_Data_Cancer)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$PeripheralBloodCancer]
      MET_Data_Cancer<-as.data.frame(MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter])
      colnames(MET_Data_Normal)=cn
      rownames(MET_Data_Normal)=rn
      if(!is.data.frame(MET_Data_Cancer)){
        MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
        if(sum(colnames(MET_Data_Cancer)!=Samplegroups$PeripheralBloodCancer[Samplegroups$PeripheralBloodCancer%in%colnames(MET_Data_Cancer)])!=0){
          colnames(MET_Data_Cancer)<-Samplegroups$PeripheralBloodCancer[Samplegroups$PeripheralBloodCancer%in%colnames(MET_Data_Cancer)]
        }
      }
    } else {
      if (CancerSite=='BRCA') {
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary,drop=FALSE]
        rn=rownames(MET_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Cancer)
        
        MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        MET_Data_Cancer<-MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter]
        barcode=colnames(MET_Data_Cancer)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(MET_Data_Cancer))
        MET_Data_Cancer<-MET_Data_Cancer[,id]
        if(!is.data.frame(MET_Data_Cancer)){
          MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
          if(sum(colnames(MET_Data_Cancer)!=Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)])!=0){
            colnames(MET_Data_Cancer)<-Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)]
          }
        }
        
      }else{
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary,drop=FALSE]
        rn=rownames(MET_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Cancer)
        
        MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        MET_Data_Cancer<-as.data.frame(MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter])
        colnames(MET_Data_Normal)=cn
        rownames(MET_Data_Normal)=rn
        if(!is.data.frame(MET_Data_Cancer)){
          MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
          if(sum(colnames(MET_Data_Cancer)!=Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)])!=0){
            colnames(MET_Data_Cancer)<-Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)]
          }
        }
      }
      
    }
    
    
    
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,28)
    
    
    
    if (CancerSite=='LAML') {
      MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal,drop=FALSE]
      rn=rownames(MET_Data_Normal)
      rn=rn[!duplicated(rn)]
      cn=colnames(MET_Data_Normal)
      
      MET_Data_Normal=as.data.frame(MET_Data_Normal[!duplicated(rownames(MET_Data_Normal)),])
      colnames(MET_Data_Normal)=cn
      rownames(MET_Data_Normal)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$BloodNormal]
      MET_Data_Normal<-as.data.frame(MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter])
      colnames(MET_Data_Normal)=Samplegroups$BloodNormal
      if(!is.data.frame(MET_Data_Normal)){
        MET_Data_Normal<-as.data.frame(MET_Data_Normal)
        if(sum(colnames(MET_Data_Normal)!=Samplegroups$BloodNormal[Samplegroups$BloodNormal%in%colnames(MET_Data_Normal)])!=0){
          colnames(MET_Data_Normal)<-Samplegroups$BloodNormal[Samplegroups$BloodNormal%in%colnames(MET_Data_Normal)]
        }
      }
    } else {
      if (CancerSite=='BRCA') {
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
        rn=rownames(MET_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Normal)
        
        MET_Data_Normal=as.data.frame(MET_Data_Normal[!duplicated(rownames(MET_Data_Normal)),])
        colnames(MET_Data_Normal)=cn
        rownames(MET_Data_Normal)=rn
        
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        MET_Data_Normal<-as.data.frame(MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter])
        colnames(MET_Data_Normal)=cn
        rownames(MET_Data_Normal)=rn
        barcode=colnames(MET_Data_Normal)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(MET_Data_Normal))
        MET_Data_Normal<-as.data.frame(MET_Data_Normal[,id])
        colnames(MET_Data_Normal)=id
        if(!is.data.frame(MET_Data_Normal)){
          MET_Data_Normal<-as.data.frame(MET_Data_Normal)
          if(sum(colnames(MET_Data_Normal)!=Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)])!=0){
            colnames(MET_Data_Normal)<-Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)]
          }
        }
        
      }else{
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
        rn=rownames(MET_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Normal)
        
        MET_Data_Normal=as.data.frame(MET_Data_Normal[!duplicated(rownames(MET_Data_Normal)),])
        colnames(MET_Data_Normal)=cn
        rownames(MET_Data_Normal)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        MET_Data_Normal<-as.data.frame(MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter])
        colnames(MET_Data_Normal)=Samplegroups$SolidNormal
        rownames(MET_Data_Normal)=rn
        if(!is.data.frame(MET_Data_Normal)){
          MET_Data_Normal<-as.data.frame(MET_Data_Normal)
          if(sum(colnames(MET_Data_Normal)!=Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)])!=0){
            colnames(MET_Data_Normal)<-Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)]
          }
        }
      }
      
    }
    
    
    
    if (length(MET_Data_Normal)>0) {
      MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,28)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    rm(list = c('MET_Data','meta_filter','METfiles','Filename','metadata','Samplegroups'))
    gc()
  }
  
  # Missing value estimation
  if(na){
    if(is.null(data)) {
      cat("\tMissing value estimation for the cancer samples.\n")
      MET_Data_Cancer=TCGA_Process_EstimateMissingValues(Data=MET_Data_Cancer,MissingValueThreshold,sex,SNP)
      if (length(MET_Data_Normal)>0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal=TCGA_Process_EstimateMissingValues(Data=MET_Data_Normal,MissingValueThreshold,sex,SNP)
      }
    }
    
    if(!is.null(data)){
      cat("\tMissing value estimation for the cancer samples.\n")
      MET_Data_Cancer=TCGA_Process_EstimateMissingValues(data,MissingValueThreshold,sex,SNP)
      
    }
  }
  
  # Batch correction for cancer and normal. 
  if(batchcor) {
    if(is.null(data)) {
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          ##BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch,tp=batchtp)
          
          if (length(MET_Data_Normal)>0) {
            cat("\tComBat Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            
          } else {
            MET_Data_Normal=c()
          }
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(GEN_Data=MET_Data_Cancer,BatchData = BatchData,MinInBatch = MinPerBatch,tp=batchtp)
          MET_Data_Cancer=MET_Data_Cancer[[1]]
          if (length(MET_Data_Normal)>0) {
            cat("\tLimma Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            MET_Data_Normal=MET_Data_Normal[[1]]
          } else {
            MET_Data_Normal=c()
          }
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          MET_Data_Cancer<-range01(MET_Data_Cancer)
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal<-range01(MET_Data_Normal)
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
          
          
        }
      }else{
        stop(message('BatchData argument cannot be null'))
      }
      
      
      
    }
    if(!is.null(data)){
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          #BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(data,BatchData)
          data=TCGA_BatchCorrection_MolecularData(
            data,BatchData,MinPerBatch,tp=batchtp)
          
          
          
          data[data<0]=0
          data[data>1]=1
          
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          data_ls=TCGA_BatchCorrection_MolecularData(
            GEN_Data=data,BatchData = BatchData,MinInBatch = MinPerBatch,
            tp=batchtp)
          data01=data_ls$GEN_Data_Corrected
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          if(is.data.frame(data01)|is.matrix(data01)){
            data01<-range01(data01)
            
            data01[data01<0]=0
            data01[data01>1]=1
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }else{
            data01<-data_ls$GEN_Data
            data01[data01<0]=0
            data01[data01>1]=1
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }
          
        }
      }
    }
  }
  
  
  
  # Reducing to 12 ids. 
  #MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,12)
  #MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,12)     
  
  
  if(is.null(data) & batchcor==F){
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
    rm(list = c('MET_Data_Cancer','MET_Data_Normal'))
    gc()
  }
  
  if(is.null(data) & batchcor==T){
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal,data_ls))
    rm(list = c('MET_Data_Cancer','=MET_Data_Normal','data_ls','data01'))
  }
  
  if(!is.null(data) & batchcor==T){
    return(data_ls)
    rm(list = c('data_ls','data01'))
    gc()
  }
  
}

##METdirectory<-METdirectories$METdirectory450k

Preprocess_CancerSite_Methylation450k <- function(data=NULL,CancerSite, directory, metadata = NULL, MissingValueThreshold = 0.2,na=TRUE, batchcor=TRUE, batchtp='limma', BatchData=NULL,sex=F,SNP=T,meth_reduce=TRUE) {
  
  ##get("BatchData")
  MinPerBatch=4
  metadata<-metadata[[1]][[1]]
  #METdirectory<-gsub('.{1}$', '', METdirectory)
  if(is.null(data)) {
    if (grepl("Windows", Sys.info()['sysname'])) {
      # If Windows I'll create a virtual drive to handle the long file names issue
      # Create a virtual drive to overcome long names issue in Windows
      METdirectory<-gsub('.{1}$', '', directory)
      virtualDir <- METdirectory
      system(paste("subst x:", virtualDir))
      
      # Load data
      METfiles <- dir("x:")
      
      METfiles = list.files(directory,
                            full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,METfiles)
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)  
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('.txt',METfiles)
      METfiles=METfiles[MatchedFilePosition]
      Filename=METfiles    
      s<-file.size(Filename)
      s<-s/1e9
      delim="\t"
      pre_process_size=1000
      df <- read_delim(Filename, delim = delim, n_max = pre_process_size,name_repair = 'minimal')
      nc<-length(colnames(df))
      if(s<1){
        MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
      }else{
        if(nc<=1000){
          cols<-colnames(df)
          csv.name=Filename
          delim="'\t'"
          fname<-paste0(CancerSite,'_Meth','_raw')
          show_progress_bar = T
          Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
          sdir=directory
          CancerSite=CancerSite
          cf<- paste0(directory,'/',CancerSite,'_cols.RData')
          rf<- paste0(directory,'/',CancerSite,'_rdcols.RData')
          MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                    delim = delim, fname = fname, 
                                    show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                    CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
        }
        
        if(nc>1000) {
          if(as.integer(nc/4)+1 <=1000) {
            cols<-colnames(df)
            csv.name=Filename
            delim="'\t'"
            fname<-paste0(CancerSite,'_Meth','_raw')
            show_progress_bar = T
            Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
            sdir=directory
            CancerSite=CancerSite
            cf<- paste0(directory,'/',CancerSite,'_cols.RData')
            rf<- paste0(directory,'/',CancerSite,'_rdcols.RData')
            MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                      delim = delim, fname = fname, 
                                      show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                      CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
          } else {
            stop(message('Number of columns is too large!'))
          }
        }
      }
    } else {
      # Not windows
      # Load data
      METfiles = list.files(directory,
                            full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,METfiles)
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)  
      METfiles=METfiles[MatchedFilePosition]
      MatchedFilePosition=grep('MANIFEST',METfiles)
      METfiles=METfiles[-MatchedFilePosition]
      MatchedFilePosition=grep('.txt',METfiles)
      METfiles=METfiles[MatchedFilePosition]
      Filename=METfiles    
      s<-file.size(Filename)
      s<-s/1e9
      delim="\t"
      pre_process_size=1000
      df <- read_delim(Filename, delim = delim, n_max = pre_process_size)
      
      #creating new colnames to avoid duplications
      cr<-colnames(df)
      cr<-t(paste0(cr,collapse = '\t'))
      cm<-paste0('sed -i  "1i ',cr,'" ',Filename)
      system(cm)
      df <- read_delim(Filename, delim = delim, n_max = pre_process_size)
      nc<-length(colnames(df))
      if(s<1){
        MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
      }else{
        if(nc<=1000){
          cols<-colnames(df)
          csv.name=Filename
          delim="'\t'"
          fname<-paste0(CancerSite,'_Meth','_raw')
          show_progress_bar = T
          Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
          sdir=directory
          CancerSite=CancerSite
          cf<- paste0(directory,CancerSite,'_cols.RData')
          rf<- paste0(directory,CancerSite,'_rdcols.RData')
          MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                    delim = delim, fname = fname, 
                                    show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                    CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
        }
        
        if(nc>1000) {
          if(as.integer(nc/4)+1 <=1000) {
            cols<-colnames(df)
            csv.name=Filename
            delim="'\t'"
            fname<-paste0(CancerSite,'_Meth','_raw')
            show_progress_bar = T
            Probes='/home/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap/PreProcess/Resources/p27k.RData'
            sdir=directory
            CancerSite=CancerSite
            cf<- paste0(directory,CancerSite,'_cols.RData')
            rf<- paste0(directory,CancerSite,'_rdcols.RData')
            MET_Data<-TCGA_fread_load(data = df,meth_reduce = T,csv.name = Filename,
                                      delim = delim, fname = fname, 
                                      show_progress_bar = T, Probes = Probes, sdir = sdir, 
                                      CancerSite = CancerSite, cf = cf, rf=rf,cols = cols,return = T)
          } else {
            stop(message('Number of columns is too large!'))
          }
        }
      }
    }
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    #Samplegroups2=TCGA_GENERIC_GetSampleGroups(x$cases)
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    #Samplegroups2=TCGA_GENERIC_GetSampleGroups(x$cases)
    if (CancerSite=='LAML') {
      MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
      
      rn=rownames(MET_Data_Cancer)
      rn=rn[!duplicated(rn)]
      cn=colnames(MET_Data_Cancer)
      
      MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
      colnames(MET_Data_Cancer)=cn
      rownames(MET_Data_Cancer)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$PeripheralBloodCancer]
      MET_Data_Cancer<-as.data.frame(MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter])
      colnames(MET_Data_Cancer)=cn
      rownames(MET_Data_Cancer)=rn
      if(!is.data.frame(MET_Data_Cancer)){
        MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
        if(sum(colnames(MET_Data_Cancer)!=Samplegroups$PeripheralBloodCancer[Samplegroups$PeripheralBloodCancer%in%colnames(MET_Data_Cancer)])!=0){
          colnames(MET_Data_Cancer)<-Samplegroups$PeripheralBloodCancer[Samplegroups$PeripheralBloodCancer%in%colnames(MET_Data_Cancer)]
        }
      }
    } else {
      if (CancerSite=='BRCA') {
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(MET_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Cancer)
        
        MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        MET_Data_Cancer<-as.data.frame(MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        barcode=colnames(MET_Data_Cancer)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(MET_Data_Cancer))
        MET_Data_Cancer<-MET_Data_Cancer[,id]
        if(!is.data.frame(MET_Data_Cancer)){
          MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
          if(sum(colnames(MET_Data_Cancer)!=Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)])!=0){
            colnames(MET_Data_Cancer)<-Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)]
          }
        }
        
      }else{
        MET_Data_Cancer=MET_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(MET_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(MET_Data_Cancer)
        
        MET_Data_Cancer=as.data.frame(MET_Data_Cancer[!duplicated(rownames(MET_Data_Cancer)),])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        MET_Data_Cancer<-as.data.frame(MET_Data_Cancer[,colnames(MET_Data_Cancer)%in%meta_filter])
        colnames(MET_Data_Cancer)=cn
        rownames(MET_Data_Cancer)=rn
        if(!is.data.frame(MET_Data_Cancer)){
          MET_Data_Cancer<-as.data.frame(MET_Data_Cancer)
          if(sum(colnames(MET_Data_Cancer)!=Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)])!=0){
            colnames(MET_Data_Cancer)<-Samplegroups$Primary[Samplegroups$Primary%in%colnames(MET_Data_Cancer)]
          }
        }
      }
      
    }
    
    
    
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,28)
    
    
    
    
    if (CancerSite=='LAML') {
      MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal,drop=FALSE]
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$BloodNormal]
      MET_Data_Normal<-as.data.frame(MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter])
      colnames(MET_Data_Normal)=Samplegroups$BloodNormal
      if(!is.data.frame(MET_Data_Normal)){
        MET_Data_Normal<-as.data.frame(MET_Data_Normal)
        if(sum(colnames(MET_Data_Normal)!=Samplegroups$BloodNormal[Samplegroups$BloodNormal%in%colnames(MET_Data_Normal)])!=0){
          colnames(MET_Data_Normal)<-Samplegroups$BloodNormal[Samplegroups$BloodNormal%in%colnames(MET_Data_Normal)]
        }
      }
    } else {
      if (CancerSite=='BRCA') {
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        MET_Data_Normal<-MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter]
        barcode=colnames(MET_Data_Normal)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(MET_Data_Normal))
        MET_Data_Normal<-as.data.frame(MET_Data_Normal[,id])
        colnames(MET_Data_Normal)=id
        if(!is.data.frame(MET_Data_Normal)){
          MET_Data_Normal<-as.data.frame(MET_Data_Normal)
          if(sum(colnames(MET_Data_Normal)!=Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)])!=0){
            colnames(MET_Data_Normal)<-Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)]
          }
        }
        
      }else{
        MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        MET_Data_Normal<-as.data.frame(MET_Data_Normal[,colnames(MET_Data_Normal)%in%meta_filter])
        colnames(MET_Data_Normal)=Samplegroups$SolidNormal
        if(!is.data.frame(MET_Data_Normal)){
          MET_Data_Normal<-as.data.frame(MET_Data_Normal)
          if(sum(colnames(MET_Data_Normal)!=Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)])!=0){
            colnames(MET_Data_Normal)<-Samplegroups$SolidNormal[Samplegroups$SolidNormal%in%colnames(MET_Data_Normal)]
          }
        }
      }
      
    }
    
  
    if (length(MET_Data_Normal)>0) {
      MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,28)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    
    
    # Clear space
    rm(list = c('MET_Data','meta_filter','METfiles','Filename','metadata','Samplegroups'))
    gc()
  }
  
  
  
  
  
  if(na){
    if(is.null(data)) {
      cat("\tMissing value estimation for the cancer samples.\n")
      MET_Data_Cancer=TCGA_Process_EstimateMissingValues(Data=MET_Data_Cancer,MissingValueThreshold,sex,SNP)
      if (length(MET_Data_Normal)>0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal=TCGA_Process_EstimateMissingValues(Data=MET_Data_Normal,MissingValueThreshold,sex,SNP)
      }
    }
    
    if(!is.null(data)){
      cat("\tMissing value estimation for the cancer samples.\n")
      data=TCGA_Process_EstimateMissingValues(data,MissingValueThreshold,sex,SNP)
      
    }
  }
  
  
  # Batch correction for cancer and normal. 
  if(batchcor) {
    if(is.null(data)) {
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch,tp=batchtp)
          
          if (length(MET_Data_Normal)>0) {
            cat("\tComBat Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            
          } else {
            MET_Data_Normal=c()
          }
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(GEN_Data=MET_Data_Cancer,BatchData = BatchData,MinInBatch = MinPerBatch,tp=batchtp)
          MET_Data_Cancer=MET_Data_Cancer[[1]]
          if (length(MET_Data_Normal)>0) {
            cat("\tLimma Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            MET_Data_Normal=MET_Data_Normal[[1]]
          } else {
            MET_Data_Normal=c()
          }
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          MET_Data_Cancer<-range01(MET_Data_Cancer)
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal<-range01(MET_Data_Normal)
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
          
          
        }
      }else{
        stop(message('BatchData argument cannot be null'))
      }
      
      
      
    }
    if(!is.null(data)){
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(data,BatchData)
          data=TCGA_BatchCorrection_MolecularData(data,BatchData,MinPerBatch,tp=batchtp)
          
          
          
          data[data<0]=0
          data[data>1]=1
          
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          data_ls=TCGA_BatchCorrection_MolecularData(GEN_Data=data,BatchData = BatchData,MinInBatch = MinPerBatch,tp=batchtp)
          data01=data_ls$GEN_Data_Corrected
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          if(is.data.frame(data01)|is.matrix(data01)){
            data01<-range01(data01)
            
            data01[data01<0]=0
            data01[data01>1]=1
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }else{
            data01<-data_ls$GEN_Data
            data01[data01<0]=0
            data01[data01>1]=1
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }
        }
      }
    }
  }
  
  
  if(is.null(data) & batchcor==F){
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
    rm(list = c('MET_Data_Cancer','=MET_Data_Normal'))
    gc()
  }
  
  if(is.null(data) & batchcor==T){
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal,data_ls))
    rm(list = c('MET_Data_Cancer','=MET_Data_Normal','data_ls','data01'))
  }
  
  if(!is.null(data) & batchcor==T){
    return(data_ls)
    rm(list = c('data_ls','data01'))
    gc()
  }
}

Preprocess_CancerSite_RNAseq <- function(data=NULL,CancerSite, directory, metadata, batchcor=TRUE, batchtp='limma', BatchData=NULL) {
  
  # Settings
  #get("BatchData")
  MinPerBatch=4
  metadata<-metadata[[1]][[1]]
  
  if(is.null(data)) {
    
    if (grepl("Windows", Sys.info()['sysname'])) {
      # If Windows I'll create a virtual drive to handle the long file names issue
      # Create a virtual drive to overcome long names issue in Windows
      directory<-directories[[1]]
      virtualDir <- directory
      virtualDir <- gsub("\\\\", "/", virtualDir)
      virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
      system(paste("subst x:", virtualDir))
      
      # Load data
      files <- dir("x:")
      files = list.files(directory,
                         full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,files)
      files=files[MatchedFilePosition]
      MatchedFilePosition=grep(paste0(CancerSite,'.uncv2.mRNAseq_raw'),files) 
      
      Filename <- paste0(files[MatchedFilePosition])          
      RNA_Data <- TCGA_GENERIC_LoadIlluminaRNAseqData(Filename,metadata)
      
      system("subst x: /D") #stop virtual drive
    } else {
      # Not windows
      # Load data
      directory<-directories[[1]]
      files = list.files(directory,
                         full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,files)
      files=files[MatchedFilePosition]
      MatchedFilePosition=grep(paste0(CancerSite,'.uncv2.mRNAseq_raw'),files)  
      files=files[MatchedFilePosition]
      
      MatchedFilePosition=grep('.txt',files)
      files=files[MatchedFilePosition]
      Filename=files         
      RNA_Data=TCGA_GENERIC_LoadIlluminaRNAseqData(Filename,metadata)
    }
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(RNA_Data))
    #Samplegroups2=TCGA_GENERIC_GetSampleGroups(x$cases)
    if (CancerSite=='LAML') {
      RNA_Data_Cancer=RNA_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
      
      rn=rownames(RNA_Data_Cancer)
      rn=rn[!duplicated(rn)]
      cn=colnames(RNA_Data_Cancer)
      
      RNA_Data_Cancer=as.data.frame(RNA_Data_Cancer[!duplicated(rownames(RNA_Data_Cancer)),])
      colnames(RNA_Data_Cancer)=cn
      rownames(RNA_Data_Cancer)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$PeripheralBloodCancer]
      RNA_Data_Cancer<-as.data.frame(RNA_Data_Cancer[,colnames(RNA_Data_Cancer)%in%meta_filter])
    } else {
      if (CancerSite=='BRCA') {
        RNA_Data_Cancer=RNA_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(RNA_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(RNA_Data_Cancer)
        
        RNA_Data_Cancer=as.data.frame(RNA_Data_Cancer[!duplicated(rownames(RNA_Data_Cancer)),])
        colnames(RNA_Data_Cancer)=cn
        rownames(RNA_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        RNA_Data_Cancer<-RNA_Data_Cancer[,colnames(RNA_Data_Cancer)%in%meta_filter]
        barcode=colnames(RNA_Data_Cancer)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(RNA_Data_Cancer))
        RNA_Data_Cancer<-RNA_Data_Cancer[,id]
        
      }else{
        RNA_Data_Cancer=RNA_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(RNA_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(RNA_Data_Cancer)
        
        RNA_Data_Cancer=as.data.frame(RNA_Data_Cancer[!duplicated(rownames(RNA_Data_Cancer)),])
        colnames(RNA_Data_Cancer)=cn
        rownames(RNA_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        RNA_Data_Cancer<-RNA_Data_Cancer[,colnames(RNA_Data_Cancer)%in%meta_filter]
      }
      
    }
    RNA_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(RNA_Data_Cancer,28)
    if (CancerSite=='LAML') {
      RNA_Data_Normal=RNA_Data[,Samplegroups$BloodNormal,drop=FALSE]
      
      rn=rownames(RNA_Data_Normal)
      rn=rn[!duplicated(rn)]
      cn=colnames(RNA_Data_Normal)
      
      RNA_Data_Normal=as.data.frame(RNA_Data_Normal[!duplicated(rownames(RNA_Data_Normal)),])
      colnames(RNA_Data_Normal)=cn
      rownames(RNA_Data_Normal)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$BloodNormal]
      RNA_Data_Normal<-RNA_Data_Normal[,colnames(RNA_Data_Normal)%in%meta_filter]
    } else {
      if (CancerSite=='BRCA') {
        RNA_Data_Normal=RNA_Data[,Samplegroups$SolidNormal,drop=FALSE]
        
        rn=rownames(RNA_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(RNA_Data_Normal)
        
        RNA_Data_Normal=as.data.frame(RNA_Data_Normal[!duplicated(rownames(RNA_Data_Normal)),])
        colnames(RNA_Data_Normal)=cn
        rownames(RNA_Data_Normal)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        RNA_Data_Normal<-RNA_Data_Normal[,colnames(RNA_Data_Normal)%in%meta_filter]
        barcode=colnames(RNA_Data_Normal)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(RNA_Data_Normal))
        RNA_Data_Normal<-RNA_Data_Normal[,id]
        
      }else{
        RNA_Data_Normal=RNA_Data[,Samplegroups$SolidNormal,drop=FALSE]
        
        rn=rownames(RNA_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(RNA_Data_Normal)
        
        RNA_Data_Normal=as.data.frame(RNA_Data_Normal[!duplicated(rownames(RNA_Data_Normal)),])
        colnames(RNA_Data_Normal)=cn
        rownames(RNA_Data_Normal)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        RNA_Data_Normal<-as.data.frame(RNA_Data_Normal[,colnames(RNA_Data_Normal)%in%meta_filter])
        colnames(RNA_Data_Normal)=cn
        rownames(RNA_Data_Normal)=rn
      }
      
    }
    if (length(RNA_Data_Normal)>0) {
      RNA_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(RNA_Data_Normal,28)
    }
    
    if (length(RNA_Data_Normal)>0) {
      RNA_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(RNA_Data_Normal,28)
    }
    cat("There are",length(colnames(RNA_Data_Cancer)),"cancer samples and",length(colnames(RNA_Data_Normal)),"normal samples in",CancerSite,"\n")
    rm(list = c('RNA_Data','meta_filter','files','Filename','metadata','Samplegroups'))
    gc()
  }
  
  # Missing value estimation
  
  
  # Batch correction for cancer and normal. 
  if(batchcor) {
    if(is.null(data)) {
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          ##BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch,tp=batchtp)
          
          if (length(MET_Data_Normal)>0) {
            cat("\tComBat Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            
          } else {
            MET_Data_Normal=c()
          }
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(GEN_Data=MET_Data_Cancer,BatchData = BatchData,MinInBatch = MinPerBatch,tp=batchtp)
          MET_Data_Cancer=MET_Data_Cancer[[1]]
          if (length(MET_Data_Normal)>0) {
            cat("\tLimma Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            MET_Data_Normal=MET_Data_Normal[[1]]
          } else {
            MET_Data_Normal=c()
          }
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          MET_Data_Cancer<-range01(MET_Data_Cancer)
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal<-range01(MET_Data_Normal)
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
          
          
        }
      }else{
        stop(message('BatchData argument cannot be null'))
      }
      
      
      
    }
    if(!is.null(data)){
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          #BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(data,BatchData)
          data=TCGA_BatchCorrection_MolecularData(
            data,BatchData,MinPerBatch,tp=batchtp)
          
          
          
          data[data<0]=0
          data[data>1]=1
          
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          data_ls=TCGA_BatchCorrection_MolecularData(
            GEN_Data=data,BatchData = BatchData,MinInBatch = MinPerBatch,
            tp=batchtp)
          data01=data_ls$GEN_Data_Corrected
          
          
          if(is.data.frame(data01)|is.matrix(data01)){
            
            
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }else{
            data01<-data_ls$GEN_Data
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }
          
        }
      }
    }
  }
  
  
  
  # Reducing to 12 ids. 
  #MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,12)
  #MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,12)     
  
  
  if(is.null(data) & batchcor==F){
    return(list(RNA_Data_Cancer=RNA_Data_Cancer,RNA_Data_Normal=RNA_Data_Normal))
    rm(list = c('RNA_Data_Cancer','RNA_Data_Normal'))
    gc()
  }
  
  if(is.null(data) & batchcor==T){
    return(list(RNA_Data_Cancer=RNA_Data_Cancer,RNA_Data_Normal=RNA_Data_Normal,data_ls))
    rm(list = c('RNA_Data_Cancer','RNA_Data_Normal','data_ls','data01'))
  }
  
  if(!is.null(data) & batchcor==T){
    return(data_ls)
    rm(list = c('data_ls','data01'))
    gc()
  }
  
}


Preprocess_CancerSite_miRNAseq <- function(data=NULL,CancerSite, directory, metadata, batchcor=TRUE, batchtp='limma', BatchData=NULL) {
  
  # Settings
  #get("BatchData")
  MinPerBatch=4
  metadata<-metadata[[1]][[1]]
  
  if(is.null(data)) {
    
    if (grepl("Windows", Sys.info()['sysname'])) {
      # If Windows I'll create a virtual drive to handle the long file names issue
      # Create a virtual drive to overcome long names issue in Windows
      directory<-directories[[1]]
      virtualDir <- directory
      virtualDir <- gsub("\\\\", "/", virtualDir)
      virtualDir <- substr(virtualDir, 1, nchar(virtualDir) - 1)
      system(paste("subst x:", virtualDir))
      
      # Load data
      files <- dir("x:")
      files = list.files(directory,
                         full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,files)
      files=files[MatchedFilePosition]
      MatchedFilePosition=grep(paste0(CancerSite,'.miRseq_raw_counts'),files) 
      
      Filename <- paste0(files[MatchedFilePosition]) 
      miRseq_Data <- TCGA_GENERIC_LoadIlluminamiRseqData(Filename,metadata)
      
      system("subst x: /D") #stop virtual drive
    } else {
      # Not windows
      # Load data
      directory<-directories[[1]]
      files = list.files(directory,
                         full.names = T,all.files = T,recursive = T)
      MatchedFilePosition=grep(CancerSite,files)
      files=files[MatchedFilePosition]
      MatchedFilePosition=grep(paste0(CancerSite,'.miRseq_raw_counts'),files)  
      files=files[MatchedFilePosition]
      
      MatchedFilePosition=grep('.txt',files)
      files=files[MatchedFilePosition]
      Filename=files         
      miRseq_Data <- TCGA_GENERIC_LoadIlluminamiRseqData(Filename,metadata)
    }
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(miRseq_Data))
    #Samplegroups2=TCGA_GENERIC_GetSampleGroups(x$cases)
    if (CancerSite=='LAML') {
      miRseq_Data_Cancer=miRseq_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
      
      rn=rownames(miRseq_Data_Cancer)
      rn=rn[!duplicated(rn)]
      cn=colnames(miRseq_Data_Cancer)
      
      miRseq_Data_Cancer=as.data.frame(miRseq_Data_Cancer[!duplicated(rownames(miRseq_Data_Cancer)),])
      colnames(miRseq_Data_Cancer)=cn
      rownames(miRseq_Data_Cancer)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$PeripheralBloodCancer]
      miRseq_Data_Cancer<-as.data.frame(miRseq_Data_Cancer[,colnames(miRseq_Data_Cancer)%in%meta_filter])
    } else {
      if (CancerSite=='BRCA') {
        miRseq_Data_Cancer=miRseq_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(miRseq_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(miRseq_Data_Cancer)
        
        miRseq_Data_Cancer=as.data.frame(miRseq_Data_Cancer[!duplicated(rownames(miRseq_Data_Cancer)),])
        colnames(miRseq_Data_Cancer)=cn
        rownames(miRseq_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        miRseq_Data_Cancer<-miRseq_Data_Cancer[,colnames(miRseq_Data_Cancer)%in%meta_filter]
        barcode=colnames(miRseq_Data_Cancer)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(miRseq_Data_Cancer))
        miRseq_Data_Cancer<-miRseq_Data_Cancer[,id]
        
      }else{
        miRseq_Data_Cancer=miRseq_Data[,Samplegroups$Primary,drop=FALSE]
        
        rn=rownames(miRseq_Data_Cancer)
        rn=rn[!duplicated(rn)]
        cn=colnames(miRseq_Data_Cancer)
        
        miRseq_Data_Cancer=as.data.frame(miRseq_Data_Cancer[!duplicated(rownames(miRseq_Data_Cancer)),])
        colnames(miRseq_Data_Cancer)=cn
        rownames(miRseq_Data_Cancer)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$Primary]
        miRseq_Data_Cancer<-as.data.frame(miRseq_Data_Cancer[,colnames(miRseq_Data_Cancer)%in%meta_filter])
        colnames(miRseq_Data_Cancer)=cn
        rownames(miRseq_Data_Cancer)=rn
      }
      
    }
    miRseq_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(miRseq_Data_Cancer,28)
    
    if (CancerSite=='LAML') {
      miRseq_Data_Normal=miRseq_Data[,Samplegroups$BloodNormal,drop=FALSE]
      
      rn=rownames(miRseq_Data_Normal)
      rn=rn[!duplicated(rn)]
      cn=colnames(miRseq_Data_Normal)
      
      miRseq_Data_Normal=as.data.frame(miRseq_Data_Normal[!duplicated(rownames(miRseq_Data_Normal)),])
      colnames(miRseq_Data_Normal)=cn
      rownames(miRseq_Data_Normal)=rn
      
      meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$BloodNormal]
      miRseq_Data_Normal<-miRseq_Data_Normal[,colnames(miRseq_Data_Normal)%in%meta_filter]
    } else {
      if (CancerSite=='BRCA') {
        miRseq_Data_Normal=miRseq_Data[,Samplegroups$SolidNormal,drop=FALSE]
        
        rn=rownames(miRseq_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(miRseq_Data_Normal)
        
        miRseq_Data_Normal=as.data.frame(miRseq_Data_Normal[!duplicated(rownames(miRseq_Data_Normal)),])
        colnames(miRseq_Data_Normal)=cn
        rownames(miRseq_Data_Normal)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        miRseq_Data_Normal<-miRseq_Data_Normal[,colnames(miRseq_Data_Normal)%in%meta_filter]
        barcode=colnames(miRseq_Data_Normal)
        barcode<-str_split(barcode,'-')
        barcode<-lapply(barcode, function(x){c(x[1:3])})
        barcode<-lapply(barcode, function(x){paste0(x,collapse = '-')})
        
        meta <- GDCquery_clinic(unique(metadata$project), type = "clinical")
        meta<-subset(meta,meta$submitter_id%in%barcode)
        meta_filter<-meta$submitter_id[meta$gender=='female']
        
        id<-grepl(paste0(meta_filter,collapse = '|'),colnames(miRseq_Data_Normal))
        miRseq_Data_Normal<-miRseq_Data_Normal[,id]
        
      }else{
        miRseq_Data_Normal=miRseq_Data[,Samplegroups$SolidNormal,drop=FALSE]
        
        rn=rownames(miRseq_Data_Normal)
        rn=rn[!duplicated(rn)]
        cn=colnames(miRseq_Data_Normal)
        
        miRseq_Data_Normal=as.data.frame(miRseq_Data_Normal[!duplicated(rownames(miRseq_Data_Normal)),])
        colnames(miRseq_Data_Normal)=cn
        rownames(miRseq_Data_Normal)=rn
        
        meta_filter<-metadata$cases[metadata$cases%in%Samplegroups$SolidNormal]
        miRseq_Data_Normal<-miRseq_Data_Normal[,colnames(miRseq_Data_Normal)%in%meta_filter]
      }
      
    }
    if (length(miRseq_Data_Normal)>0) {
      miRseq_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(miRseq_Data_Normal,28)
    }
    
   
    cat("There are",length(colnames(miRseq_Data_Cancer)),"cancer samples and",length(colnames(miRseq_Data_Normal)),"normal samples in",CancerSite,"\n")
    rm(list = c('miRseq_Data','meta_filter','files','Filename','metadata','Samplegroups'))
    gc()
  }
  
  # Missing value estimation
  
  
  # Batch correction for cancer and normal. 
  if(batchcor) {
    if(is.null(data)) {
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          ##BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch,tp=batchtp)
          
          if (length(MET_Data_Normal)>0) {
            cat("\tComBat Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            
          } else {
            MET_Data_Normal=c()
          }
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(GEN_Data=MET_Data_Cancer,BatchData = BatchData,MinInBatch = MinPerBatch,tp=batchtp)
          MET_Data_Cancer=MET_Data_Cancer[[1]]
          if (length(MET_Data_Normal)>0) {
            cat("\tLimma Batch correction for the normal samples.\n")
            BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
            MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch,tp=batchtp)
            MET_Data_Normal=MET_Data_Normal[[1]]
          } else {
            MET_Data_Normal=c()
          }
          
          range01<-function(x){(x-min(x))/(max(x)-min(x))}
          MET_Data_Cancer<-range01(MET_Data_Cancer)
          
          MET_Data_Cancer[MET_Data_Cancer<0]=0
          MET_Data_Cancer[MET_Data_Cancer>1]=1
          if (length(MET_Data_Normal)>0) {
            MET_Data_Normal<-range01(MET_Data_Normal)
            MET_Data_Normal[MET_Data_Normal<0]=0
            MET_Data_Normal[MET_Data_Normal>1]=1
          }
          
          
        }
      }else{
        stop(message('BatchData argument cannot be null'))
      }
      
      
      
    }
    if(!is.null(data)){
      if(!is.null(BatchData)) {
        if(batchtp=='ComBat'){
          cat("\tComBat Batch correction for the cancer samples.\n")
          
          #BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(data,BatchData)
          data=TCGA_BatchCorrection_MolecularData(
            data,BatchData,MinPerBatch,tp=batchtp)
          
          
          
          data[data<0]=0
          data[data>1]=1
          
        }
        
        if(batchtp=='limma'){
          
          cat("\tLimma Batch correction for the cancer samples.\n")
          data_ls=TCGA_BatchCorrection_MolecularData(
            GEN_Data=data,BatchData = BatchData,MinInBatch = MinPerBatch,
            tp=batchtp)
          data01=data_ls$GEN_Data_Corrected
          
          
          if(is.data.frame(data01)|is.matrix(data01)){
            
            
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }else{
            data01<-data_ls$GEN_Data
            data_ls$GEN_Data_Corrected<-data01
            data_ls$GEN_Data_Corrected<-data_ls$GEN_Data_Corrected[,order(
              colnames(data_ls$GEN_Data_Corrected))]
            data_ls$GEN_Data<-data[,order(colnames(data))]
          }
          
        }
      }
    }
  }
  
  
  
  # Reducing to 12 ids. 
  #MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,12)
  #MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,12)     
  
  
  if(is.null(data) & batchcor==F){
    return(list(miRseq_Data_Cancer=miRseq_Data_Cancer,miRseq_Data_Normal=miRseq_Data_Normal))
    rm(list = c('mirRseq_Data_Cancer','mirRseq_Data_Normal'))
    gc()
  }
  
  if(is.null(data) & batchcor==T){
    return(list(miRseq_Data_Cancer=miRseq_Data_Cancer,miRseq_Data_Normal=miRseq_Data_Normal,data_ls))
    rm(list = c('RNA_Data_Cancer','RNA_Data_Normal','data_ls','data01'))
  }
  
  if(!is.null(data) & batchcor==T){
    return(data_ls)
    rm(list = c('data_ls','data01'))
    gc()
  }
  
}
#data_matrix=tmp
#batch_label=batch$Name
DSC <- function(data_matrix, batch_label) {
  # Dispersion within batches
  Sw <- scale(t((data_matrix)), center = TRUE, scale = FALSE)**2 / dim(data_matrix)[2] # vectorized
  Dw <- sqrt(sum(as.numeric(Sw)))
  #Dw <- sqrt(sum(sapply(split(data.frame(t(data_matrix)), batch_label), function(x) sum(diag(cov(x))))))
  # Dispersion between batches
  M = apply(data_matrix, 1, mean)
  mean_var <- function(x) {o <- apply(x, 2, mean) - M; return(sum(o**2))} # vectorized
  #mean_var <- function(x) {o <- apply(x, 2, mean) - M; return(sum(diag(o%*%t(o))))}
  Sb <- sapply(split(data.frame(t(data_matrix)), batch_label), mean_var)
  Db = sqrt(sum(Sb))
  dsc<-(Db/Dw)
  rm(list = c('Sw','Dw','M','mean_var','Sb','Db'))
  gc()
  return(dsc)
}

#data=MET_Data_Cancer
#metainfo=metainfo
batch_detect<-function(tp=c('batch','center','sex','TSS'), data=NULL , metainfo=NULL){
  if(is.null(data)){
    stop(message('Data argument must be a data.frame or matrix'))
  }
  if(is.null(metainfo)){
    stop(message('metainfo argument cannot be null'))
  }
  dsc<-c()
  
  if(length(colnames(data))>=2){
    for (i in seq_along(tp)) {
      print(tp[i])
      tmp<-data
      tmp<-na.omit(tmp) 
      batch<-snames(metainfo,colnames(data),tp=tp[i])
      colnames(tmp)<-as.character(batch$Name)
      if(length(colnames(tmp))>=2){
        d<-DSC(tmp,colnames(tmp))
        dsc<-append(dsc,d)
      }
    }
    dsc<-data.frame('Covariate'=tp,'DSC'=dsc)
  }else{
    print('Data has only 1 sample, skipping DSC calculation.')
    dsc<-data.frame('Covariate'=tp,'DSC'=rep(0,length(tp)))
  }
  
  
  return(dsc)
}


#data=x
#data=MET_Data_Normal
#data2<-data[,!colnames(data)%in%'T11']
#data2<-data[,-267]
PCA<-function(data,batch='batch',main='Title',meta=metainfo,outrm=T){
  
  if(length(colnames(data))<10){
    return(list('Dataset is too small',Outlier_rm=0))
  }else {
    
    ori<-colnames(data)
    colnames(data)<-snames(meta,colnames(data),batch)$Name
    
    pcomp<-function(x){
      x<-x[ , which(apply(x, 2, var) != 0)]
      x<-x[ which(apply(x, 1, var) != 0),]
      
      pca <- prcomp(t(x), scale=T, center=TRUE) 
      
      ## calculate the percentage of variation that each PC accounts for...
      pca.var <- pca$sdev^2
      pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
      pca.var.per
      
      ## now make a fancy looking plot that shows the PCs and the variation:
      pca.data <- data.frame(Sample=rownames(pca$x),
                             X=pca$x[,1],
                             Y=pca$x[,2])
      
      return(list(pca.data=pca.data,pca.var.per=pca.var.per))
    rm(list = c('pca','pca.var','pca.var.per','pca.data'))
    gc()
    }
    
    
    
    out_rm<-function(x){
      s<-abs(x[,2])+abs(x[,3])
      
      k = round(length(x[,1])/100)
      if(k>0){
        test <- rosnerTest(s,
                           k = k
        )
        
        bad<-test$all.stats$Obs.Num[test$all.stats$Outlier==T]
        return(bad)
      }else{
        test <- rosnerTest(s,
                           k = 1
        )
        
        bad<-test$all.stats$Obs.Num[test$all.stats$Outlier==T]
        return(bad)
        rm(list = c('s','k','test','bad'))
        gc()
      }
      
      
      
    }
    
    if(outrm) {
      pca<-pcomp(data)
      pca.data<-pca$pca.data
      pca.var.per<-pca$pca.var.per
      bad<-out_rm(pca.data)
      orm<-0
      while(length(bad)>0){
        data<-data[,-bad]
        ori<-ori[-bad]
        orm<-orm+length(bad)
        cat(paste0(length(bad),' outliers were removed!'))
        pca<-pcomp(data)
        pca.data<-pca$pca.data
        pca.var.per<-pca$pca.var.per
        bad<-out_rm(pca.data)
      }
      
      colnames(data)<-ori  
      
      return(list(Data=data,Outlier_rm=orm))
      rm(list = c('pca','pca.data','pca.var.ver','bad','orm','ori','data'))
      gc()
      
    }else{
      
      pca<-pcomp(data)
      pca.data<-pca$pca.data
      pca.var.per<-pca$pca.var.per
      a<-ggplot(data=pca.data, aes(x=X, y=Y, color=Sample)) +
        geom_point() +
        xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
        ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
        theme_bw() + 
        ggtitle(main)
      return(a)
      rm(list = c('pca','pca.data','pca.var.ver','a'))
      gc()
      
    }
  }
  
  
}

  
  
  
  
  
  
  
  


#data=GEN_Data
den<-function(data,main='Title'){
  data=as.matrix(data)
  d <- density(data,na.rm=TRUE)
  plot(d, main=main)
  polygon(d, col="red", border="blue")
}

tcga_projects<-function(){
  url <- "https://api.gdc.cancer.gov/projects?size=1000&format=json"
  json <- fromJSON(content(GET(url), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
  projects <- json$data$hits
  projects$tumor <- unlist(lapply(projects$project_id, function(x) {
    unlist(str_split(x, "-"))[2]
  }))
  
  
  projects <- projects$project_id
  projects <- projects[grepl('TCGA',projects)]
  projects <- sort(projects)
  return(projects)
}

