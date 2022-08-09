library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
require(DOSE)
library(enrichplot)
library(pathview)
library(stringr)
library(openxlsx)
library(msigdbr)
library(stringr)

hs1<- createStyle(halign = "CENTER", textDecoration = "bold",border = "TopBottomLeftRight",fontSize=12)

organism="org.Hs.eg.db"
library(organism, character.only = TRUE)

######################################################################################################

Datos=c()
Comparativas=c() 
TestComparativas=c()

######################################### T1 vs T2

Comparativa = "T1_vs_T2"
Comparando = "T2_over_T1"
TestComparativas[Comparativa] = Comparando
Datos[[Comparativa]]=c()
d = read.delim("../DGE/Time_1_VS_Time_2_All_Corrected.txt", "\t",header=TRUE)
Datos[[Comparativa]][["Upregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05
Comparativas[length(Comparativas)+1]=Comparativa

d$EntrezID <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
d$Uniprot <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "UNIPROT", keytype = "ENSEMBL", multiVals = "first")
d<-d[!duplicated(d[,"EntrezID"]),]
d<-d[!is.na(d[,"EntrezID"]),]
row.names(d)<-d[,"EntrezID"]
d<-d[order(d$Log2.FC.,decreasing=T),]

Datos[[Comparativa]][["All_Entrez"]]<-d
Datos[[Comparativa]][["Significative_Entrez"]]<-d[(d[,"pValue"]<=0.05),]
Datos[[Comparativa]][["Upregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05

######################################### T3 vs T2

Comparativa = "T3_vs_T2"
Comparando = "T2_over_T3"
TestComparativas[Comparativa] = Comparando
Datos[[Comparativa]]=c()
d = read.delim("../DGE/Time_3_VS_Time_2_All_Corrected.txt", "\t",header=TRUE)
Datos[[Comparativa]][["Upregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05
Comparativas[length(Comparativas)+1]=Comparativa

d$EntrezID <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
d$Uniprot <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "UNIPROT", keytype = "ENSEMBL", multiVals = "first")
d<-d[!duplicated(d[,"EntrezID"]),]
d<-d[!is.na(d[,"EntrezID"]),]
row.names(d)<-d[,"EntrezID"]
d<-d[order(d$Log2.FC.,decreasing=T),]

Datos[[Comparativa]][["All_Entrez"]]<-d
Datos[[Comparativa]][["Significative_Entrez"]]<-d[(d[,"pValue"]<=0.05),]
Datos[[Comparativa]][["Upregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05

######################################### T1 vs T3

Comparativa = "T1_vs_T3"
Comparando = "T3_over_T1"
TestComparativas[Comparativa] = Comparando
Datos[[Comparativa]]=c()
d = read.delim("../DGE/Time_1_VS_Time_3_All_Corrected.txt", "\t",header=TRUE)
Datos[[Comparativa]][["Upregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05
Comparativas[length(Comparativas)+1]=Comparativa

d$EntrezID <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
d$Uniprot <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "UNIPROT", keytype = "ENSEMBL", multiVals = "first")
d<-d[!duplicated(d[,"EntrezID"]),]
d<-d[!is.na(d[,"EntrezID"]),]
row.names(d)<-d[,"EntrezID"]
d<-d[order(d$Log2.FC.,decreasing=T),]

Datos[[Comparativa]][["All_Entrez"]]<-d
Datos[[Comparativa]][["Significative_Entrez"]]<-d[(d[,"pValue"]<=0.05),]
Datos[[Comparativa]][["Upregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]>=0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Downregulated_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (d[,"Log2.FC."]<=-0.58),] ##pValue o FDR...0.05
Datos[[Comparativa]][["Altered_Entrez"]]<-d[(d[,"pValue"]<=0.05) & (abs(d[,"Log2.FC."])>=0.58),] ##pValue o FDR...0.05



###################################################################################################################

######################################### Qué hacer

keyType_Ont="ENSEMBL"

toTest=c("Upregulated", "Downregulated")

GO_ORA=T
onthologies=c("BP", "MF", "CC")
nameOnthologies=c()
nameOnthologies["BP"]="Biological Process"
nameOnthologies["CC"]="Celular Component"
nameOnthologies["MF"]="Molecular Function"

KEGG_GSE=F
MKEGG_GSE=T
KEGG_ORA=F
MKEGG_ORA=T
KEGG_Organism="hsa"
pathways_KEGG=T
Kegg_Pathways=c("hsa00900")

Reactome_GSE=F
Reactome_ORA=T
organism_Reactome="human" #"human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"


WP_GSE=F
WP_ORA=F
WP_Organism="Homo sapiens"

MSIGDB_ORA=F
MSIGDB_GSE=F
#msigdbr_show_species()
MSIG_Organism="Homo sapiens"
#Aux<-msigdbr_collections()
#View(Aux)
MSIG_Cats<-c("C1","C2","C3","C3", "C8", "C7")
MSIG_SubCats<-c("","CP:BIOCARTA", "TFT:GTRD", "TFT:TFT_Legacy","","IMMUNESIGDB")



######## Ejecucion
for(k in 1:length(Comparativas)){
  Comparativa = Comparativas[k]
  Comparando = TestComparativas[[Comparativa]]
  if (! file.exists(Comparativa)){
    dir.create(Comparativa)
  }
  
  fileName <- paste(Comparativa,".xlsx",sep="")
  excel <- createWorkbook(fileName,creator="Enrique Vazquez")
  
  if(GO_ORA){
    for(j in 1:length(toTest)){
      testing = toTest[j]
      
      d2 = Datos[[Comparativa]][[testing]]
      geneList = d2[,5]
      names(geneList) = as.character(d2[,1])
      geneList = sort(geneList, decreasing = TRUE)
      
      for(i in 1:length(onthologies)){
        onthology = onthologies[i]
        
        print(paste(onthology,testing,Comparativa))
        
        Aux_Enrichment <- enrichGO(names(geneList), OrgDb = organism, ont=onthology, readable=TRUE, keyType = keyType_Ont)
        Aux_Enrichment_Simplified <- simplify(Aux_Enrichment)
        
        nTerms=nrow(as.data.frame(Aux_Enrichment))
        nTerms_Simplified=nrow(as.data.frame(Aux_Enrichment_Simplified))
        
        sheet=str_sub(paste(testing,onthology),1,31)
        
        if(nTerms > 0){
          if(nTerms_Simplified > 0){
            addWorksheet(excel, sheet)
            writeData(excel, Aux_Enrichment_Simplified, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
            
            paraGuardar<-dotplot(Aux_Enrichment_Simplified,showCategory=10, title=paste(testing," ",Comparando,"\n",nameOnthologies[onthology],sep=""), orderBy="GeneRatio") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + theme(plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1), plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1), axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.text.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.title.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))
          }
          else {
            sheet=paste(testing,onthology,"_NS")
            addWorksheet(excel, sheet)
            writeData(excel, as.data.frame(Aux_Enrichment), sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
            
            paraGuardar<-dotplot(Aux_Enrichment,showCategory=10, title=paste(testing," ",Comparando,"\n",nameOnthologies[onthology],sep=""), orderBy="GeneRatio") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + theme(plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1), plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1), axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.text.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.title.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))
            
          }
          png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 3000, height = 2000, units = "px",res=300)
          print(paraGuardar)
          dev.off()
          
          m=pairwise_termsim(Aux_Enrichment,showCategory=5000)
          if(nrow(as.data.frame(m)) > 1){
            paraGuardar<-emapplot(m,showCategory=50)
            png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,"_emaplot.HiRes.png",sep=""), width = 3000, height = 2500, units = "px",res=300)
            print(paraGuardar)
            dev.off()
            
            paraGuardar <- tryCatch(
              {
                paraGuardar <- emapplot_cluster(m,showCategory=100)
              },
              error = function(cond){
                paraGuardar <- "NO"
              },
              finally = {
                #pass
              })
            
            
            if(! (paraGuardar == "NO")){
              png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,"_emaplot_cluster.HiRes.png",sep=""), width = 3000, height = 2500, units = "px",res=300)
              print(paraGuardar)
              dev.off()
            }
          }
        } else {
          print("No hay terminos que imprimir")
        }
        
      }
    }
  }
  
  if(pathways_KEGG){
    geneList<-Datos[[Comparativa]][["All_Entrez"]][,"Log2.FC."]
    names(geneList)<-row.names(Datos[[Comparativa]][["All_Entrez"]])
    geneList2<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
    names(geneList2)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
    dir.create(paste(Comparativa, "/Kegg",sep=""))
    for(i in 1:length(Kegg_Pathways)){
      wd <- getwd()
      od <- paste(wd,"/",Comparativa,"/Kegg/",sep="")
      setwd(od)
      
      pathway<-Kegg_Pathways[i]
      
      aux <- pathview(gene.data  = geneList, pathway.id = pathway, species = KEGG_Organism, out.suffix="All", limit = list(gene=2, cpd=1), bins = list(gene=30, cpd=10))
      aux <- pathview(gene.data  = geneList2, pathway.id = pathway, species = KEGG_Organism, out.suffix="Significative", limit = list(gene=2, cpd=1), bins = list(gene=30, cpd=10))
      
      setwd(wd)
    }
    
  }
  
  if(KEGG_GSE){
    testing="GSEA"
    onthology="KEGG"
    
    print(paste(testing,onthology, Comparativa))
    geneList<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
    names(geneList)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
    geneList = sort(geneList, decreasing = TRUE)
    
    Aux_Enrichment <- gseKEGG(geneList = geneList, organism = KEGG_Organism, minGSSize = 50, pvalueCutoff = 0.05, verbose = FALSE)
    nTerms=nrow(as.data.frame(Aux_Enrichment))
    
    sheet=str_sub(paste(testing,onthology),1,31)
    
    if(nTerms > 0){
      
      Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
      
      addWorksheet(excel, sheet)
      writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
      paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
      png(paste(Comparativa,"/",onthology,"_",testing,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
      print(paraGuardar)
      dev.off()
    }
  }
  
  if(MKEGG_GSE){
    testing="GSEA"
    onthology="MKEGG"
    
    print(paste(testing,onthology, Comparativa))
    geneList<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
    names(geneList)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
    geneList = sort(geneList, decreasing = TRUE)
    
    Aux_Enrichment <- gseMKEGG(geneList = geneList, organism = KEGG_Organism, pvalueCutoff = 0.05)
    nTerms=nrow(as.data.frame(Aux_Enrichment))
    
    sheet=str_sub(paste(testing,onthology),1,31)
    
    if(nTerms > 0){
      
      Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
      
      addWorksheet(excel, sheet)
      writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
      
      paraGuardar<-dotplot(Aux_Enrichment,showCategory=10, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + theme(plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1), plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1), axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.text.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.title.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))
      png(paste(Comparativa,"/",onthology,"_",testing,".HiRes.png",sep=""), width = 3000, height = 2000, units = "px",res=300)
      print(paraGuardar)
      dev.off()
    }
  }
  
  if(Reactome_GSE){
    testing="GSEA"
    onthology="Reactome"
    
    print(paste(testing,onthology, Comparativa))
    geneList<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
    names(geneList)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
    geneList = sort(geneList, decreasing = TRUE)
    
    Aux_Enrichment <- gsePathway(geneList, pvalueCutoff = 0.2, pAdjustMethod = "BH", organism=organism_Reactome, verbose = FALSE)
    nTerms=nrow(as.data.frame(Aux_Enrichment))
    
    sheet=str_sub(paste(testing,onthology),1,31)
    
    if(nTerms > 0){
      
      Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
      
      addWorksheet(excel, sheet)
      writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
      
      paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
      png(paste(Comparativa,"/",onthology,"_",testing,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
      print(paraGuardar)
      dev.off()
    }
  }
  
  if(WP_GSE){
    testing="GSEA"
    onthology="Wikipathways"
    
    print(paste(testing,onthology, Comparativa))
    geneList<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
    names(geneList)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
    geneList = sort(geneList, decreasing = TRUE)
    
    Aux_Enrichment <- gseWP(geneList, organism = WP_Organism)
    nTerms=nrow(as.data.frame(Aux_Enrichment))
    
    sheet=str_sub(paste(testing,onthology),1,31)
    
    if(nTerms > 0){
      
      Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
      
      addWorksheet(excel, sheet)
      writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
      
      paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
      png(paste(Comparativa,"/",onthology,"_",testing,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
      print(paraGuardar)
      dev.off()
    }
  }
  
  if(KEGG_ORA){
    testing="ORA"
    onthology="KEGG"
    for(j in 1:length(toTest)){
      testing = toTest[j]
      print(paste(testing,onthology, Comparativa))
      
      d2 = Datos[[Comparativa]][[paste(testing,"_Entrez",sep="")]]
      geneList = d2[,5]
      names(geneList) = row.names(d2)
      geneList = sort(geneList, decreasing = TRUE)
    
      Aux_Enrichment <- enrichKEGG(gene = names(geneList), organism = KEGG_Organism, pvalueCutoff = 0.05)
      nTerms=nrow(as.data.frame(Aux_Enrichment))
      
      sheet=str_sub(paste(testing,onthology),1,31)
      
      if(nTerms > 0){
        
        Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
        
        addWorksheet(excel, sheet)
        writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
        paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
        png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
        print(paraGuardar)
        dev.off()
      }
    }
  }
    
  if(MKEGG_ORA){
    testing="ORA"
    onthology="MKEGG"
    for(j in 1:length(toTest)){
      testing = toTest[j]
      print(paste(testing,onthology, Comparativa))
      
      d2 = Datos[[Comparativa]][[paste(testing,"_Entrez",sep="")]]
      geneList = d2[,5]
      names(geneList) = row.names(d2)
      geneList = sort(geneList, decreasing = TRUE)
      
      Aux_Enrichment <- enrichMKEGG(gene = names(geneList), organism = KEGG_Organism, pvalueCutoff = 0.05, qvalueCutoff = 1)
      nTerms=nrow(as.data.frame(Aux_Enrichment))
      
      sheet=str_sub(paste(testing,onthology),1,31)
      
      if(nTerms > 0){
        
        Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
        
        addWorksheet(excel, sheet)
        writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
        paraGuardar<-dotplot(Aux_Enrichment,showCategory=10, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + theme(plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1), plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1), axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.text.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.title.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))
        png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 3000, height = 2000, units = "px",res=300)
        print(paraGuardar)
        dev.off()
      }
    }
  }  
    
  if(Reactome_ORA){
    testing="ORA"
    onthology="Reactome"
    for(j in 1:length(toTest)){
      testing = toTest[j]
      print(paste(testing,onthology, Comparativa))
      
      d2 = Datos[[Comparativa]][[paste(testing,"_Entrez",sep="")]]
      geneList = d2[,5]
      names(geneList) = row.names(d2)
      geneList = sort(geneList, decreasing = TRUE)
      
      Aux_Enrichment <- enrichPathway(gene = names(geneList), pvalueCutoff = 0.05, organism=organism_Reactome)
      nTerms=nrow(as.data.frame(Aux_Enrichment))
      
      sheet=str_sub(paste(testing,onthology),1,31)
      
      if(nTerms > 0){
        
        Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
        
        addWorksheet(excel, sheet)
        writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
        paraGuardar<-dotplot(Aux_Enrichment,showCategory=10, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio") + scale_y_discrete(labels=function(x) str_wrap(x, width=40)) + theme(plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1), plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1), axis.text.y = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.text.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), axis.title.x = element_text(angle = 0, size = 16, face = 'bold', vjust = 0.5), legend.title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 14))
        png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 3000, height = 2000, units = "px",res=300)
        print(paraGuardar)
        dev.off()
        
        svg(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.svg",sep=""), width = 11.7, height = 8.3)
        print(paraGuardar)
        dev.off()
        
        pdf(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.pdf",sep=""), width = 11.7, height = 8.3)
        print(paraGuardar)
        dev.off()
        
      }
    }
  }  
   
  if(WP_ORA){
    testing="ORA"
    onthology="WikiPathways"
    for(j in 1:length(toTest)){
      testing = toTest[j]
      print(paste(testing,onthology, Comparativa))
      
      d2 = Datos[[Comparativa]][[paste(testing,"_Entrez",sep="")]]
      geneList = d2[,5]
      names(geneList) = row.names(d2)
      geneList = sort(geneList, decreasing = TRUE)
      
      Aux_Enrichment <- enrichWP(gene = names(geneList), organism = WP_Organism, pvalueCutoff = 0.05)
      nTerms=nrow(as.data.frame(Aux_Enrichment))
      
      sheet=str_sub(paste(testing,onthology),1,31)
      
      if(nTerms > 0){
        
        Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
        
        addWorksheet(excel, sheet)
        writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
        paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
        png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
        print(paraGuardar)
        dev.off()
      }
    }
  } 
  
  if(MSIGDB_ORA){
    testing="ORA"
    for(h in 1:length(MSIG_Cats)){
      cat<-MSIG_Cats[h]
      subcat<-MSIG_SubCats[h]
      onthology=str_replace(paste("MSIG",cat,subcat,sep="_"), ":", "_")
      
      Aux_MSig<-msigdbr(species = MSIG_Organism, category = cat, subcategory = subcat) %>% dplyr::select(gs_name, entrez_gene)
      
      for(j in 1:length(toTest)){
        testing = toTest[j]
        print(paste(testing,onthology, Comparativa))
        d2 = Datos[[Comparativa]][[paste(testing,"_Entrez",sep="")]]
        geneList = d2[,5]
        names(geneList) = row.names(d2)
        geneList = sort(geneList, decreasing = TRUE)
        
        Aux_Enrichment <- enricher(gene = names(geneList), TERM2GENE=Aux_MSig, pvalueCutoff = 0.05)
        nTerms=nrow(as.data.frame(Aux_Enrichment))
        
        sheet=str_sub(paste(testing,onthology),1,31)
        
        if(nTerms > 0){
          
          Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
          
          addWorksheet(excel, sheet)
          writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
          
          paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
          png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
          print(paraGuardar)
          dev.off()
        }
      }
    }
  }
  
  if(MSIGDB_GSE){
    testing="GSEA"
    for(h in 1:length(MSIG_Cats)){
      cat<-MSIG_Cats[h]
      subcat<-MSIG_SubCats[h]
      onthology=str_replace(paste("MSIG",cat,subcat,sep="_"), ":", "_")
      
      Aux_MSig<-msigdbr(species = MSIG_Organism, category = cat, subcategory = subcat) %>% dplyr::select(gs_name, entrez_gene)
      
      testing="GSE"
      print(paste(testing,onthology, Comparativa))
      
      geneList<-Datos[[Comparativa]][["Significative_Entrez"]][,"Log2.FC."]
      names(geneList)<-row.names(Datos[[Comparativa]][["Significative_Entrez"]])
      geneList = sort(geneList, decreasing = TRUE)
      
      Aux_Enrichment <- GSEA(gene = geneList, TERM2GENE=Aux_MSig, pvalueCutoff = 0.05)
      nTerms=nrow(as.data.frame(Aux_Enrichment))
      
      sheet=str_sub(paste(testing,onthology),1,31)
      
      if(nTerms > 0){
        
        Aux_Enrichment <- setReadable(Aux_Enrichment, OrgDb = organism,keyType = "ENTREZID")
        
        addWorksheet(excel, sheet)
        writeData(excel, Aux_Enrichment, sheet=sheet, rowNames = FALSE, borders = "columns", headerStyle = hs1)
        
        paraGuardar<-dotplot(Aux_Enrichment,showCategory=50, title=paste(testing," ",Comparando,"\n",onthology,sep=""), orderBy="GeneRatio")
        png(paste(Comparativa,"/",onthology,"_",testing,"_",Comparando,".HiRes.png",sep=""), width = 5000, height = 3500, units = "px",res=300)
        print(paraGuardar)
        dev.off()
      }
        
    }
  }
  
  
  saveWorkbook(excel, file = fileName, overwrite = TRUE)
  
}

columns(org.Hs.eg.db)
d$Accnum <- mapIds(x = org.Hs.eg.db, keys = d[,1], column = "UNIPROT", keytype = "ENSEMBL", multiVals = "first")
for(cat in columns(org.Hs.eg.db)){
  d[[cat]]<- mapIds(x = org.Hs.eg.db, keys = d[,1], column = cat, keytype = "ENSEMBL", multiVals = "first")
}
View(d)
cat="ALIAS"
d[[cat]]<- mapIds(x = org.Hs.eg.db, keys = d[,1], column = cat, keytype = "ENSEMBL", multiVals = "list")
cat="UNIPROT"
d[[cat]]<- mapIds(x = org.Hs.eg.db, keys = d[,1], column = cat, keytype = "ENSEMBL", multiVals = "list")
cat="SYMBOL"
d[[cat]]<- mapIds(x = org.Hs.eg.db, keys = d[,1], column = cat, keytype = "ENSEMBL", multiVals = "list")
cat="GENENAME"
d[[cat]]<- mapIds(x = org.Hs.eg.db, keys = d[,1], column = cat, keytype = "ENSEMBL", multiVals = "list")
