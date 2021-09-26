source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')
source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_function.R')
load('C:/Users/giordano/PP_params.RData')
ProcessedDataRNA<-RNAseq_data_make(CancerSite,directories,lg,
                                       dtype='RNAseq',
                                       mvt=mvt,bct,normal)
save(ProcessedDataRNA,file=paste0(svtmp,'ProcessedDataRNA.Rdata'))
