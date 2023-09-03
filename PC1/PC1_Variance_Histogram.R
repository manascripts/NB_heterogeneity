library(factoextra)
library(ggplot2)
library(data.table)
library(readxl)
library(ggpubr)
library(grid)

    type = "Neuro"
    #Directories for plots to be generated, datasets, and gene lists
    plots = paste0("./", type,"/plots")
    Datasets = paste0("./", type, "/Datasets")
    gene_lists = paste0("./", type,"/gene_lists")
    
    #Input list of datasets from folder
    setwd(Datasets)
    ds = list.files(pattern =".txt")
    num = "NOR"
    #Input list of housekeeping genes
    setwd(gene_lists)
    hk_gene = toupper(read.delim("Housekeeping genes.txt")[,1])
    hk_gene = gsub(" ", "", hk_gene)
    
    eig_list = list()
    eig_list2 = list()
    NOR_gene = as.data.frame(read_excel("./genes.xlsx",sheet = 2))
    NOR_gene = NOR_gene[NOR_gene$Phenotype == "NOR",]
    NOR_gene = gsub(" ", "", NOR_gene[,1])
    MES_gene = as.data.frame(read_excel("./genes.xlsx",sheet = 2))
    MES_gene = MES_gene[MES_gene$Phenotype == "MES",]
    MES_gene = gsub(" ", "", MES_gene[,1])
    
    
    g = num 
    #Perform clustering 
    for (g in num){
      y = which(num == g)    #for each genelist
      eig_list[[y]] = matrix(nrow = 8000, ncol = (length(ds) + 1))
      eig_list2[[y]] = matrix(nrow = 8000, ncol = (length(ds) + 1))
      
      for (i in ds){
        
        k = which(ds==i)
        j = gsub("_gene-exp.txt", "", i)
        
        #set wd to folder with datasets
        setwd(Datasets)
        
        data<-read.table(paste0(i), header = T, sep = "\t")
        if(colnames(data)[2]=="Probe_ID")
        {
          print("Yes")
          data<-data[,-c(2)]
        }
        data = data[!duplicated(data[,1]),]
        row.names(data) = data[,1]
        data = data[,-1]
        
        setwd(gene_lists)
        gene_list = as.data.frame(read_excel(paste0(g, ".xlsx"),sheet=1))[,1]
        
        data_wo_gene = data[-which(row.names(data) %in% gene_list),]
        data_gene = data[which(row.names(data) %in% gene_list),]
        data_hk_gene = data[which(row.names(data) %in% hk_gene),]
        
        
        for(z in 1:1000) #loop for 1000 unique swaps 
        {
          
          RandomNum = sample(1:nrow(data_hk_gene), 26, replace = FALSE) #select a random index of a gene from HK matrix
          
          data_hk_gene1 = data_hk_gene[RandomNum,]
          data1 = na.omit(data_hk_gene1); rm(data_hk_gene1)
          data1 <- as.data.frame(t(data1))
          data1 = data1[, colSums(data1 != 0) > 0]
          data1 <- scale(data1)
          data.pca <- prcomp(data1)
          
          # add variance for each swap
          eig_list[[y]][z,k] =  (data.pca$sdev[1]/sum(data.pca$sdev))*100 
          # number of genes swapped 
          eig_list[[y]][z,(length(ds) + 1)] = z 
        }

        dfr = na.omit(as.data.frame(eig_list[[y]])[,c(k,(length(ds) + 1))]) 
        names(dfr) = c("PC1V", "genes")

        setwd(gene_lists)
        nor_gene = as.data.frame(read_excel("./NOR.xlsx",sheet = 1))
        nor_gene = nor_gene[nor_gene$Phenotype == "NOR",]
        nor_gene = gsub(" ", "", nor_gene[,1])
        mes_gene = as.data.frame(read_excel("./NOR.xlsx",sheet = 1))
        mes_gene = mes_gene[mes_gene$Phenotype == "MES",]
        mes_gene = gsub(" ", "", mes_gene[,1])

        data_nor = data[which(row.names(data) %in% nor_gene),]
        data_mes = data[which(row.names(data) %in% mes_gene),]
        data_NOR = data[which(row.names(data) %in% NOR_gene),]
        data_MES = data[which(row.names(data) %in% MES_gene),]
        
        
        for(z in 1:1000) #loop for 1000 unique swaps 
        {
          
          
          RandomNOR = sample(1:nrow(data_NOR), nrow(data_nor), replace = FALSE)  
          RandomMES = sample(1:nrow(data_MES), nrow(data_mes), replace = FALSE) 
          data_NOR1 = data_NOR[RandomNOR,]   
          data_MES1 = data_MES[RandomMES,] 
          
          data_NOR1 = na.omit(data_NOR1); data_MES1 = na.omit(data_MES1)
          data2 = rbind(data_NOR1, data_MES1)
          data2 <- as.data.frame(t(data2))
          data2 = data2[, colSums(data2 != 0) > 0]
          data2 <- scale(data2) 
          data.pca <- prcomp(data2)
          
          eig_list2[[y]][z,k] =  (data.pca$sdev[1]/sum(data.pca$sdev))*100
          eig_list2[[y]][z,(length(ds) + 1)] = z 
        }
        
        
        dfr2 = na.omit(as.data.frame(eig_list2[[y]])[,c(k,(length(ds) + 1))]) 
        names(dfr2) = c("PC1V", "genes")
        
        data_gene1 = na.omit(data_gene)
        data_gene1 <- as.data.frame(t(data_gene1))
        data_gene1 = data_gene1[, colSums(data_gene1 != 0) > 0]
        data_gene1 <- scale(data_gene1)
        data.pca1 <- prcomp(data_gene1)
        PC1V1 <- (data.pca1$sdev[1]/sum(data.pca1$sdev))*100 
        
        percentile <- sum(dfr$PC1V <= PC1V1)/nrow(dfr)
        dfr$genes <- as.factor(dfr$genes)
        dfr2$genes <- as.factor(dfr2$genes)
        pval <- t.test(dfr2$PC1V, dfr$PC1V)$p.value
        unlab = gsub("_gene-exp.txt", "", i)
        significant_text <- ifelse(p_value < 2.2e-16, "p < 2.2e-16", "ns")
        
        ggplot(NULL) +
          geom_histogram(data = dfr, aes(x = PC1V, y = ..density..), colour = "blue", fill = "#0276ab", alpha = 0.55, bins = 100) +
          geom_histogram(data = dfr2, aes(x = PC1V, y = ..density..), colour = "red", fill = "#cb3b38", alpha = 0.55, bins = 100) +
          geom_density(alpha = 0.86, fill ="#89CFF0", outline.type = "full", colour = NA) +
          ggtitle(j) +
          xlab("PC1 Variance") +    
          ylab("Density") +
          geom_vline(xintercept = PC1V1, colour = "red", size = 1) +
          theme_classic() +
          theme(plot.title = element_text(size = 26, hjust = 0.5, face = "bold", colour = "#222222")) +
          theme(axis.text = element_text(size = 18, face = "bold"),
                panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                axis.title = element_text(size = 24,face = "bold")) +
          coord_cartesian(ylim =c(0, 0.45), xlim = c(9, 26), expand = F) +
          annotate("text", x = (PC1V1-0.1*PC1V1), y = 0.39, label = significant_text, color = "black", size = 7.2, fontface = "bold")
        
        
        setwd(plots)
        ggsave(paste0(j,"_", " PC1_1000_swaps.png"), height = 7, width = 13)
      }
      
    }
