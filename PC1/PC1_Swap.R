library(factoextra)
library(ggplot2)
library(data.table)
library(readxl)
library(ggpubr)
library(grid)

    type = "Neuro"
    
    
    
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
    
    
    #Perform clustering 
    g = num
    
    for (g in num){
      y = which(num == g)    #for each genelist
      eig_list[[y]] = matrix(nrow = 8000, ncol = (length(ds) + 1))
      
      for (i in ds){
        
        k = which(ds==i)
        j = gsub("_gene-exp.txt", "", i)
        
        #set wd to folder with datasets
        setwd(Datasets)
        data<-read.table(paste0(i), header = T, sep = "\t")
        
        if(colnames(data)[2] == "Probe_ID")
        {
          print("Yes")
          data <- data[,-c(2)]
        }
        
        data = data[!duplicated(data[,1]),]
        row.names(data) = data[,1]
        data = data[,-1]
        
        setwd(gene_lists)
        gene_list = as.data.frame(read_excel(paste0(g, ".xlsx"),sheet=1))[,1]
        
        data_wo_gene = data[-which(row.names(data) %in% gene_list),]
        data_gene = data[which(row.names(data) %in% gene_list),]
        data_hk_gene = data[which(row.names(data) %in% hk_gene),]
        
        
        for(p in (1:nrow(data_gene))){          #loop for all genes in the geneset  
          for(z in ((((p-1)*100)+1) :(((p-1)*100)+100))){     #loop for 1 to 100 swaps 
            
            RandomNum = sample(1:nrow(data_hk_gene), p, replace = FALSE)    #select a random index of a gene from HK matrix
            RN = sample(1:nrow(data_gene), (nrow(data_gene)-(p)), replace = FALSE) #select a random number of indices from the actual matrix leaving one
            
            data_hk_gene1 = data_hk_gene[RandomNum,]   
            data_gene1 = data_gene[RN,]
            data1 = rbind(data_gene1, data_hk_gene1)   
            
            data1 = na.omit(data1)
            data1 <- as.data.frame(t(data1))
            data1 = data1[, colSums(data1 != 0) > 0]
            data1 <- scale(data1) 
            data.pca <- prcomp(data1)
            
            eig_list[[y]][z,k] =  (data.pca$sdev[1]/sum(data.pca$sdev))*100  # add variance for each swap
            eig_list[[y]][z,(length(ds) + 1)] = p   # number of genes swapped  # for 3rd dataset z,4
          }
        }
        
        dfr = na.omit(as.data.frame(eig_list[[y]])[,c(k,(length(ds) + 1))])   # for 3rd dataset z,4
        names(dfr) = c("PC1V", "genes")
        dfr$PC1V <- dfr$PC1V/100
        dfr$genes <- round(dfr$genes/26,6)
        means <- aggregate(dfr, by = list(dfr[,2]), FUN = mean)
        names(means)[1] <- "Swaps"
        corr <- cor.test(means$Swaps, means$PC1V, method = "pearson")
        means <- means[,-3]
        
        dfr$genes <- as.factor(dfr$genes)
        unlab = gsub("_gene-exp.txt", "", i)

        corr$p.value
        corr$estimate
        
        ggplot(dfr, aes(x = genes, y = PC1V, fill = genes)) + 
          geom_boxplot(show.legend = FALSE)+
          ggtitle(j)+
          xlab("Swaps")+
          ylab("Variance explained (PC1)")+
          theme_classic()+
          theme(plot.title = element_text(size = 26, hjust = 0.5, face = "bold", colour = "#222222"))+
          theme(axis.text = element_text(size = 14, face = "bold"),
                panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
                axis.title = element_text(size = 24,face = "bold"))
        
        setwd(plots)
        ggsave(paste0(unlab,"_", " PC1_Boxplot.png"), height = 7, width = 12)
      }
      
    }
    
    
    
