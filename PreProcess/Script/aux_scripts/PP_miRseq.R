source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_function.R')
source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')
load('C:/Users/giordano/PP_params.RData')
ProcessedDatamiRseq<-miRseq_data_make(CancerSite,directories,lg,
                                       dtype='miRseq',
                                       mvt=mvt,bct,normal)
save(ProcessedDatamiRseq,file=paste0(svtmp,'ProcessedDatamiRseq.Rdata'))
