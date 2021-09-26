source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')
source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_function.R')
load('C:/Users/giordano/PP_params.RData')

ProcessedData27k<-Methylation_data_make(CancerSite,directories,lg,
                                          dtype='Meth_27k',
                                          mvt=mvt,bct,normal)


save(ProcessedData27k,file=paste0(svtmp,'ProcessedData27k.Rdata'))


