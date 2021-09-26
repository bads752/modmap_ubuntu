library(readr)
library(data.table)


cm<-commandArgs(trailingOnly=TRUE)

csv.name=cm[1]
print('cm[1]')
print(cm[1])

delim=cm[2]
print('cm[2]')
print(cm[2])
if(delim=='tab'){
  delim='\t'
}

fname<-cm[3]
print('cm[3]')
print(cm[3])

show_progress_bar = eval(parse(text = cm[4]))
print('cm[4]')
print(cm[4])

Probes=cm[5]
print('cm[5]')
print(cm[5])

sdir=cm[6]
print('cm[6]')
print(cm[6])

CancerSite=cm[7]
print('cm[7]')
print(cm[7])

cf=cm[8]
print('cm[8]')
print(cm[8])

meth_reduce = eval(parse(text = cm[9]))
print('cm[9]')
print(cm[9])

rf=cm[10]
print('cm[10]')
print(cm[10])

load(cf)
load(rf)
load(Probes)


if(meth_reduce){
  ref<-rd
}else{
  ref<-cols
}

dt_scm = fread(csv.name, select = rd[1], 
               showProgress = show_progress_bar,sep = delim)

ids<-match(p27k,dt_scm$`Hybridization REF`)


n<-as.integer(length(ref)/100)










for(i in seq_len(n+1)) {
  if(i==1){
    tmp <- fread(csv.name, select = rd[1:100], 
                 showProgress = show_progress_bar,sep = delim)
    
    tmp<-subset(tmp,tmp$`Hybridization REF` %in%p27k)
    
    dt<-tmp
  }
  
  if(i>1&i<(n+1)){
    s<-((i-1)*100)+1
    e=i*100
    
    tmp <- fread(csv.name, select = rd[s:e], 
                 showProgress = show_progress_bar,sep = delim)

    tmp$`Hybridization REF`<-dt_scm$`Hybridization REF`
    tmp<-subset(tmp,tmp$`Hybridization REF` %in%p27k)
    tmp2<-subset(tmp,select=-c(`Hybridization REF`))
    tmp<-tmp2
  }
  
  if(i==(n+1)){
    s<-((i-1)*100)+1
    e=length(ref)
    tmp <- fread(csv.name, select = rd[s:e], 
                 showProgress = show_progress_bar,sep = delim)
    
    tmp$`Hybridization REF`<-dt_scm$`Hybridization REF`
    
    tmp<-subset(tmp,tmp$`Hybridization REF` %in%p27k)
    tmp<-subset(tmp,select=-c(`Hybridization REF`))
  }
  
  
  if(i>1){
    dt<-cbind(dt,tmp)
  }
  print(paste0('Finished parsing ',i*100,' columns!'))
  
}






n<-fname
assign(c(n),dt)
save(list = c(n),file=paste0(sdir,'/',n,'.RData')) 
