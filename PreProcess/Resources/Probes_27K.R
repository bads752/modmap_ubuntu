
load('/home/giordano/Methylation/Pre_Process/ACC_Data.RData')
load('/home/giordano/Methylation/Pre_Process/READ_Data.RData')



p1<-rownames(ACC_Data$ProcessedData450k$n1$ACC_Cancer_data_meth_450k)
p2<-rownames(READ_Data$ProcessedData450k$n1$READ_Cancer_data_meth_450k)
p3<-rownames(READ_Data$ProcessedData450k$n2$READ_Normal_data_meth_450k)
p4<-rownames(READ_Data$ProcessedData27k$n1$READ_Cancer_data_meth_27k)
p5<-rownames(READ_Data$ProcessedData27k$n2$READ_Normal_data_meth_27k)

p27k<-Reduce(intersect,list(p1,p2,p3,p4,p5))

save(p27k,file = '/home/giordano/Methylation/p27k.RData')






