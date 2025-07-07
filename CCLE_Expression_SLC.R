require(tidyr)
require(dplyr)
require(ggplot2)
require(stringr)


setwd("E:/XSu/CCLE")

DataInput.Expression <- read.csv(file="Expression_RPKM.csv",check.names = FALSE)
DataInput.Expression.SLC <- DataInput.Expression %>% filter(str_detect(Description, "SLC."))
DataInput.Expression.SLC.Long <- DataInput.Expression.SLC %>%
  pivot_longer(cols=3:ncol(DataInput.Expression),names_to="Cell.Line",values_to="RPKM")

#Target.Gene <- "IDO2"

#Target.Gene.Expression <- DataInput.Expression %>% 
#  filter(Description == Target.Gene) %>% 
#  pivot_longer(cols=3:ncol(DataInput.Expression),names_to="CCLE_ID",values_to="RPKM") %>%
#  filter(RPKM >0) %>%
#  mutate(Scaled_Exp=scale(log(RPKM,2)))

DataInput.Metabolites <- read.csv(file="CCLE_metabolomics_20190502.csv",check.names = FALSE)
DataInput.Metabolites.Long <- DataInput.Metabolites %>% pivot_longer(cols=3:ncol(DataInput.Metabolites),names_to="Metabolites",values_to="Abundance")

###############


Gene.List <- DataInput.Expression.SLC$Description
Association.Table <- matrix(0,nrow=length(Gene.List),ncol=length(unique(DataInput.Metabolites.Long$Metabolites)))

system.time(
for(k in 1:length(unique(DataInput.Metabolites.Long$Metabolites))) {
  Target.Metabolite <- colnames(DataInput.Metabolites)[k+2]
  
  Target.Metabolite.Abundance <- DataInput.Metabolites.Long %>% 
    filter(Metabolites == Target.Metabolite) %>% 
    select(CCLE_ID,Metabolites,Abundance)

    for(i in 1:length(Gene.List)) {
      Target.Gene <- Gene.List[i]
      
      Target.Gene.Expression <- DataInput.Expression.SLC %>% 
        filter(Description == Target.Gene) %>% 
        pivot_longer(cols=3:ncol(DataInput.Expression.SLC),names_to="CCLE_ID",values_to="RPKM") %>%
        filter(RPKM >0) %>%
        mutate(Scaled_Exp=scale(log(RPKM,2)))
      if(nrow(Target.Gene.Expression)>10) {
        Target.Association <- left_join(Target.Metabolite.Abundance,Target.Gene.Expression,by="CCLE_ID") %>%
          filter(!is.na(RPKM))
        
        Association <- summary(lm(Abundance ~ Scaled_Exp, data=Target.Association))
        Association.Table[i,k] <- Association$coefficients[2,3]
      }
    }
})

Metabolite.Associations <- data.frame(Gene.List,Association.Table)
colnames(Metabolite.Associations) <- c("Gene",colnames(DataInput.Metabolites)[3:227])
write.csv(Metabolite.Associations,file=paste("Metabolite_Association_","SLC",".csv",sep=""))

Metabolite.Associations.Long <- Metabolite.Associations %>% pivot_longer(cols=2:ncol(Metabolite.Associations),names_to="Metabolites",values_to="Association")
write.csv(Metabolite.Associations.Long,file=paste("Metabolite_Association_","SLC.Long",".csv",sep=""))
