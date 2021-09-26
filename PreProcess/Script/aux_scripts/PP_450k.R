source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_lib.R')
source('C:/Users/giordano/Trabalho/2019-2022/Doutorado/Development/ws-modmap/modmap_02/PreProcess/Script/aux_scripts/source_function.R')
load('C:/Users/giordano/PP_params.RData')
ProcessedData450k<-Methylation_data_make(CancerSite = CancerSite, 
                                       directories = directories,
                                       lg = lg,
                                       dtype='Meth_450k',
                                       mvt=mvt,
                                       bct = bct,
                                       normal = normal)
save(ProcessedData450k,file=paste0(svtmp,'ProcessedData450k.Rdata'))