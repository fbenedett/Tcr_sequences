
# Function to process bulk and single cell TCR sequencing
Process_TCR <- function(TCR_files){
  
  Bulk_or_SC <- function(file){
    test1<-read.delim(file)
    test2<-read.csv(file)
    if(ncol(test1)>ncol(test2)){return("bulk")}
    if(ncol(test2)>ncol(test1)){return("sc")}
  }
  format_bulk <- function(bulk){
    bulk[]<-lapply(bulk,as.character) # Transform factors into character
    ##### the following line should be always replaced with the full path
    Conversion<-read.delim("/media/fbenedet/Elements/scratch/my_cyto_r_utils/tcr_dbase/Conversion_Bulk_10X_Database.txt")
    
    # Remove Outofframe and Stop
    if(length(which(bulk$Frame %in% c("OUT","STOP")))>0){bulk<-bulk[-which(bulk$Frame %in% c("OUT","STOP")),]}
     
    #remove undefined V and J chains
    if(length(grep("_|undefined",bulk$TRBJ))>0){bulk<-bulk[-grep("_|undefined",bulk$TRBJ),]}
    if(length(grep("_|undefined",bulk$TRBV))>0){bulk<-bulk[-grep("_|undefined",bulk$TRBV),]}
    
    #Format the chain as in the single-cell data
    bulk$TRBV<-Conversion$Sc_sequencing[match(bulk$TRBV,Conversion$Bulk_sequencing)]
    bulk$TRBJ<-Conversion$Sc_sequencing[match(bulk$TRBJ,Conversion$Bulk_sequencing)]
    
    #Remove last character of the CD3 sequence to use the same nomenclature than 10X
    bulk$CDR3_aaseq<-substr(bulk$CDR3_aaseq,1,nchar(bulk$CDR3_aaseq)-1)
    
    # Create unique ID (V + J + CD3)
    bulk$unique_id<-paste(bulk$TRBV,bulk$TRBJ,bulk$CDR3_aaseq,sep="_")
    
    # Remove Chains with counts equaling 1 (like Raphael Genolet does)
    bulk$Count<-as.numeric(bulk$Count)
    if(length(which(bulk$Count==1))>0){bulk<-bulk[-which(bulk$Count==1),]}
    
    return(bulk)
  }
  format_10X <- function(VDJ){
    
    VDJ[]<-lapply(VDJ,as.character) # Transform factors into character
    
    # Remove the few multi-mapping chains
    if(length(which(VDJ$chain=="Multi"))>0){VDJ<-VDJ[-which(VDJ$chain=="Multi"),]} 
    
    # Remove non-productive CD3 (out of frame)
    if(length(which(VDJ$productive=="False"))>0){VDJ<-VDJ[-which(VDJ$productive=="False"),]} 
    
    # Keep those that are attributed to a cell
    if(length(which(VDJ$is_cell=="False"))>0){VDJ<-VDJ[-which(VDJ$is_cell=="False"),]} 
    
    # Remove the cells that have 4 chains
    duplets<-names(which(table(VDJ$barcode)==4))
    if(length(duplets)>0){VDJ<-VDJ[-which(VDJ$barcode %in% duplets),]}
    
    # Keep the total number of cells for percentage
    VDJ$total_number_cells<-length(unique(VDJ$barcode))
  
    # Separate alpha and beta chains in two files
    VDJ_a<-VDJ[which(VDJ$chain=="TRA"),]
    VDJ_b<-VDJ[which(VDJ$chain=="TRB"),]
    
    # Create unique ID
    VDJ_a$unique_id<-paste(VDJ_a$v_gene,VDJ_a$j_gene,VDJ_a$cdr3,sep="_")
    VDJ_b$unique_id<-paste(VDJ_b$v_gene,VDJ_b$j_gene,VDJ_b$cdr3,sep="_")
    return(list("alpha"=VDJ_a,"beta"=VDJ_b))
  }
  
  n<-length(TCR_files)
  technology<-c()
  for(i in 1:n){
    techno<-Bulk_or_SC(TCR_files[[i]])
    technology<-c(technology,techno)
    names(technology)[i]<-names(TCR_files)[i]
    if(techno=="bulk"){
      tmp<-read.delim(TCR_files[[i]])
      tmp<-format_bulk(tmp)
      assign(paste0(names(TCR_files)[i]),tmp)
    }
    if(techno=="sc"){
      tmp<-read.csv(TCR_files[[i]])
      tmp<-format_10X(tmp)
      assign(paste0(names(TCR_files)[i]),tmp$beta)
      assign(paste0(names(TCR_files)[i],"_alpha"),tmp$alpha)
    }
  }
  
  ####################################
  ### Merge Beta-chains data #########
  ####################################
  
  unique_chain<-c()
  for(i in 1:n){unique_chain<-c(unique_chain,get(names(TCR_files)[i])$unique_id)}
  unique_chain<-unique(unique_chain)
  
  #Create Data Frame for comparisons
  TCR_b_comp<-data.frame("unique_id"=unique_chain)
  TCR_b_comp[]<-lapply(TCR_b_comp,as.character)
  
  # Initializing columns
  for(i in 1:n){TCR_b_comp$tmp<-"no";colnames(TCR_b_comp)[which(colnames(TCR_b_comp)=="tmp")]<-paste0("is_in_",names(TCR_files)[i])}
  for(i in 1:n){TCR_b_comp$tmp<-0;colnames(TCR_b_comp)[which(colnames(TCR_b_comp)=="tmp")]<-paste0("counts_in_",names(TCR_files)[i])}

  for(i in 1:nrow(TCR_b_comp)){
    for(j in 1:n){
      tmp<-get(names(TCR_files)[j])
      tmp<-tmp[which(tmp$unique_id==TCR_b_comp$unique_id[i]),]
      if(nrow(tmp)>0){
        TCR_b_comp[i,which(colnames(TCR_b_comp)==paste0("is_in_",names(TCR_files)[j]))]<-"yes"
        if(technology[j]=="bulk"){
          TCR_b_comp[i,which(colnames(TCR_b_comp)==paste0("counts_in_",names(TCR_files)[j]))]<-sum(tmp$Count,na.rm = T)
        }
        if(technology[j]=="sc"){
          TCR_b_comp[i,which(colnames(TCR_b_comp)==paste0("counts_in_",names(TCR_files)[j]))]<-length(unique(tmp$barcode))
        }
      }
    }
  }
  
  #Compute Percentages
  for(i in 1:n){
    if(technology[i]=="bulk"){
      TCR_b_comp$tmp<-(TCR_b_comp[,which(colnames(TCR_b_comp)==paste0("counts_in_",names(TCR_files)[i]))]/sum(TCR_b_comp[,which(colnames(TCR_b_comp)==paste0("counts_in_",names(TCR_files)[i]))]))*100
      colnames(TCR_b_comp)[which(colnames(TCR_b_comp)=="tmp")]<-paste0("Perc_in_",names(TCR_files)[i])
    }
    if(technology[i]=="sc"){
      tmp<-get(names(TCR_files)[i])
      TCR_b_comp$tmp<-(TCR_b_comp[,which(colnames(TCR_b_comp)==paste0("counts_in_",names(TCR_files)[i]))]/unique(tmp$total_number_cells))*100
      colnames(TCR_b_comp)[which(colnames(TCR_b_comp)=="tmp")]<-paste0("Perc_in_",names(TCR_files)[i])
    }
  }

  #Create Object to be returned
  Processed_TCR<-list()
  for(i in 1:n){
    if(technology[i]=="bulk"){
      Processed_TCR$tmp<-get(names(TCR_files)[i])
      names(Processed_TCR)[which(names(Processed_TCR)=="tmp")]<-names(TCR_files)[i]
    }
    if(technology[i]=="sc"){
      Processed_TCR$tmp<-get(paste0(names(TCR_files)[i],"_alpha"))
      names(Processed_TCR)[which(names(Processed_TCR)=="tmp")]<-paste0(names(TCR_files)[i],"_alpha")
      Processed_TCR$tmp<-get(paste0(names(TCR_files)[i]))
      names(Processed_TCR)[which(names(Processed_TCR)=="tmp")]<-paste0(names(TCR_files)[i],"_beta")
    }
  }
  Processed_TCR$TCR_b_comp<-TCR_b_comp
  
  return(Processed_TCR)
}
