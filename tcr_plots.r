library(ggplot2)
library(dplyr)
library("stringr")
library(data.table)
library(scales)
library(Rcpp)
library(devtools)
pkgbuild::has_rtools(TRUE)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -std=c++14")
Sys.setenv("PKG_LIBS" = "-fopenmp")


## CPP function to calculate the Gini Index (much faster than lawstat)
cppFunction("
double my_gini(const NumericVector& gph){

size_t size=gph.size();

double avr=0, gini=0;
for(size_t i=0; i<size; ++i)
  avr+=gph[i];
avr/=double(size);

for(size_t i=0; i<size; ++i){
  for(size_t j=0; j<size; ++j){
    gini+= fabs(gph[i]-gph[j]);
  }
}

return(gini/(2.0*size*size*avr));

}")


## read a TCR file taking only certain lines, in particular the line with "IN"
## and where the count are bigger that 1
cleanse=function(aname){
  df=read.csv(aname, sep="\t")
  df=df[df[,2]>1 & df[,5]=="IN" & df[,4]!="undefined,_,_",]
  df$CDR3_aaseq=as.character(df[,6])
  df=df[df$CDR3_aaseq!="", ]
  df=df[!grepl("X", df$CDR3_aaseq), ]
  
  df$frequency=df$Count/sum(df$Count)
  df$trxvseqtrxj=paste(df[,3], df$CDR3_aaseq, df[,4], sep="_")
  
  # order them from the most frequent to the less frequent
  df=df[order(df$frequency, decreasing=TRUE), ]
  
  return(df)
}

## trace a list of sequences among different dataframe associating,
## to each sequence, the corresponding frequence.
## If the sequence is not found the frequence is ZERO
## Return a dataframe containing the  1..n sequence, 1..n patient_id, 1..n frequency value
## the columns of the dataframe are: "patient, frequency, and trxcseqtrxj
trace_seq=function(ldataf, lseq, sampletype){
  tddf=data.frame(patient=NULL, frequency=NULL, trxvseqtrxj=NULL)
  for(i in 1:length(ldataf)){
    for(j in lseq){
      
      if(length(unique(ldataf[[i]]$trxvseqtrxj==j))==2 ){
        tddf=rbind(tddf, ldataf[[i]][ldataf[[i]]$trxvseqtrxj==j,])
      }
      else{
        # message(j)
        t=data.frame(patient=ldataf[[i]]$patient[1], frequency=0, trxvseqtrxj=as.character(j))
        tddf=rbind(tddf, t)
      }
    }
  }
  tddf$patient=factor(tddf$patient, levels=sampletype, ordered=TRUE)
  return(tddf)
}


## trace the first "n" sequences of a given file trough different files
## if a sequence is not found in one of the file we try the next one to
## keep "n" constant
## ldataf is the list of dataframes where we should look
## startdf is the CURRENT dataframe (the one where we try to find the first "nseq")
## nseq is the number of sequences we would like to find
find_n_common_seq=function(ldataf, nseq, startdf){
  curr=1
  found_seq=0
  list_found_seq=NULL
  while( found_seq<nseq & curr<nrow(startdf) ){  # while the found_seq are less than the nseq
    #taking the current sequence
    aseq=startdf$trxvseqtrxj[curr]
    # we need to check if the current sequence is in all the files
    all_found=TRUE
    for(j in 1:length(ldataf)){
      if(aseq %in% ldataf[[j]]$trxvseqtrxj==FALSE ){
        all_found=FALSE
      }
    }
    if(all_found){
      list_found_seq=c(list_found_seq, aseq)
      found_seq=found_seq+1
    }
    curr=curr+1
  }
  
  return(list_found_seq)
}




## makes violin plots tracing the sequences from tddf
violin_tcr=function(databind, tddf, title){
  counts=databind %>% group_by(patient) %>% tally
  ggplot(databind, aes(x=patient, y=frequency))+
    xlab("Sample") +
    ylab("Frequency (log10)") +
    geom_violin(scale="width")+
    geom_text(data=counts, aes(label=n, y=max(databind$frequency)*2),
              position=position_dodge(0.9) )+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_line(data=tddf, aes(x=patient, y=frequency, group=trxvseqtrxj, colour=trxvseqtrxj))+
    scale_color_discrete(name = "TCR sequence") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)), 
      limits=c(0.000001,0.5)
    ) +
    annotation_logticks(sides="l") +
    geom_point(size=0.5) +
    stat_summary(fun.y=median, geom="point", size=4, color="red")
}

## makes violin plots tracing the sequences from tddf, but in this case
## the lines are smaller and transparent.  Usually used when there are a lot of sequences
violin_tcr_lot=function(databind, tddf, title){
  counts=databind %>% group_by(patient) %>% tally
  ggplot(databind, aes(x=patient, y=frequency))+
    xlab("Sample") +
    ylab("Frequency (log10)") +
    geom_violin(scale="width")+
    geom_text(data=counts, aes(label=n, y=max(databind$frequency)*2),
              position=position_dodge(0.9) )+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_line(data=tddf, aes(x=patient, y=frequency, group=trxvseqtrxj), size=0.3, alpha=0.3, color="red")+
    scale_color_discrete(name = "TCR sequence") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)), 
      limits=c(0.000001,0.5)
    ) +
    annotation_logticks(sides="l") +
    geom_point(size=0.5) +
    stat_summary(fun.y=median, geom="point", size=4, color="red")
}


# violin plot of the data plotting also the statistics
base_violin=function(databind, title){
  counts=databind %>% group_by(patient) %>% tally
  p=ggplot(databind, aes(x=patient, y=frequency))+
    xlab("Sample") +
    ylab("Frequency (log10)") +
    geom_violin(scale="width")+
    geom_text(data=counts, aes(label=n, y=max(databind$frequency)*2),
              position=position_dodge(0.9) )+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides="l") +
    geom_point(size=0.5) +
    stat_summary(fun.y=median, geom="point", size=4, color="red")+
    theme(text = element_text(size=20))
    return(p)
}

## generate a violin plot with no points in the middle
base_violin_nopoints=function(databind, title){
  counts=databind %>% group_by(patient) %>% tally
  p=ggplot(databind, aes(x=patient, y=frequency))+
    xlab("Sample") +
    ylab("Frequency (log10)") +
    geom_violin(scale = "width")+
    geom_text(data=counts, aes(label=n, y=max(databind$frequency)*2),
              position=position_dodge(0.9) )+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)), 
      limits=c(0.000001,0.5)
    ) +
    annotation_logticks(sides="l") +
    stat_summary(fun.y=median, geom="point", size=4, color="red")+
    theme(text = element_text(size=20))
  return(p)
}

# violin plot of the data, no statistics
base_violin_nostat=function(databind, title){
  p=ggplot(databind, aes(x=patient, y=frequency))+
    xlab("Sample") +
    ylab("Frequency (log10)") +
    geom_violin(scale="width")+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)), 
      limits=c(0.000001,0.5)
    ) +
    annotation_logticks(sides="l") +
    geom_point(size=0.5)+
    stat_summary(fun.y=median, geom="point", size=4, color="red")
  return(p)
}


## from a list of files and a list of identifier (the identifier explain the origin of the file and what they are)
## generate a dataframe that will be used for the violin plots, it will also generate
## the columns for the frequency and the concatenation between trxv-cdr3-trxj
## it remove the sequences with only one entry, it keeps only the "IN" frame

create_df=function(lfiles, sampletype){
  databind=NULL
  for(i in 1:length(lfiles)){
    t=cleanse(lfiles[i])
    t=t[, c("trxvseqtrxj", "frequency")]
    t$patient=rep(sampletype[i], nrow(t))
    databind=rbind(databind,t)
  }
  databind$patient=factor(databind$patient, levels=sampletype, ordered=TRUE)
  return(databind)
}




## calculate the Shannon entropy of a TCR dataframe
## for dataframe with ONE SINGLE SAMPLE
entropy=function(adf){
  res=-sum(adf$frequency*log2(adf$frequency))
  return(res)
}

## calculate how many TCR sequences have a frequency "thresholds times" above the median
## for dataframe with ONE SINGLE SAMPLE
n_over_median=function(adf, athreshold){
  ## these dataframe should be already ordered but just to be sure...
  order_res=sort(adf$frequency, decreasing = TRUE)
  ntotal=nrow(adf)
  amedian=median(order_res)
  return(sum(order_res>athreshold*amedian)/ntotal)
}

## calculate how many TCR sequences (fraction of TCR) are needed to reach a given threshold
## starting from the most frequent to the lowest   (the threshold should not be bigger than 1)
## for dataframe with ONE SINGLE SAMPLE
how_many_for=function(adf, athreshold){
  order_res=sort(adf$frequency, decreasing = TRUE)
  ntotal=nrow(adf)
  res=0
  for(i in 1:length(order_res)){
    res=res+order_res[i]
    if(res>athreshold){
      return(i/ntotal)
    }
  }
  return(-1)
}

## compute the clonality of a dataframe containing ONE SINGLE SAMPLE
## this is also called 1-Pielou index
## definition from:
## "TCR Repertoire Analysis Reveals Mobilization of Novel CD8+ T Cell
##  Clones Into the Cancer-Immunity Cycle Following Anti-CD4 Antibody Administration"
clonality=function(adf){
  res=sum(adf$frequency*log(adf$frequency))
  res=  1+ (res/log(nrow(adf)))
  return(res)
}


convert_seq=function(x){
  astr=strsplit(x, "_")[[1]]
  return(paste(astr[1], astr[3], astr[2], sep="_")) 
}

## convert David Barras dataframe to ggplot format for the violin plots
convert_barras=function(barr_df){
  
  ## count how many samples there are in the dataframe
  sampletypes=(grepl("is_in_", colnames(barr_df)))
  num_samples=sum(sampletypes)
  sampletypes=colnames(barr_df)[sampletypes]
  sampletypes=str_replace(sampletypes, "is_in_", "")
  
  
  ## counting columns start at num_samples+2, frequencies start at num_samples*2 + 2
  fdf=NULL
  for(i in 1:num_samples ){
    type_col=i+1
    count_col=num_samples+type_col
    
    tdf=barr_df[ barr_df[,count_col]>0,c(1, count_col) ]
    tdf[,2]=tdf[,2]/sum(tdf[,2])
    tdf$patient=rep(sampletypes[i],nrow(tdf) )
    tdf$unique_id=sapply(tdf$unique_id, convert_seq)
    colnames(tdf)=c("trxvseqtrxj","frequency","patient")
    fdf=rbind(fdf, tdf)
  }
  return(fdf)
}



### estimate the number of common sequences between N samples
### show their relative weight in each sample
weight_of_shared=function(databind, sampletype){
  ## finding the common sequences among all samples
  seq_vect=list()
  count=1
  for(k in sampletype){
    tdf=databind[databind$patient==k,"trxvseqtrxj"]
    seq_vect[[k]]=tdf
  }
  
  seq_inters=Reduce(intersect, seq_vect)
  
  ## computing their relative weight on the sample
  rel_weight=data.frame(samplename=character(), rel_freq=numeric())
  for(k in sampletype){
    tdf=databind[databind$patient==k,]
    ares=data.frame(samplename=k,rel_fre=sum(tdf[tdf$trxvseqtrxj %in% seq_inters, "frequency"]))
    rel_weight=rbind(rel_weight, ares)
  }
  return(rel_weight)
}



### plot Entropy, Clonality and Gini index for all the different tissue/patient contained
### in a dataframe
dev_tcr_estimators=function(databind,sampletype){
  astatdf=setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("entropy", "clonality", "gini", "sampletype"))

  for(k in sampletype){
    adf=databind[databind$patient==k,]
    en=entropy(adf)
    cl=clonality(adf)
    gi=my_gini(adf$frequency)
    adf=data.frame(en, cl, gi, k)
    astatdf=rbind(astatdf, adf)
  }
  
  rel_weight=weight_of_shared(databind, sampletype)
  
  astatdf=cbind(relative_weight=rel_weight$rel_fre, astatdf)
  
  colnames(astatdf)=c("Shared freq. (weighted)", "Entropy", "Clonality", "Gini", "sampletype")
  astatdf2=melt(astatdf, id.vars="sampletype")
  p=ggplot(astatdf2, aes(x=sampletype, y=value, group=variable, color=variable))+
    geom_line()+
    geom_point()+
    geom_text(aes(label=round(value, 2)),hjust=0.5, vjust=-1, show.legend = FALSE)+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits=c(0.01,50) )+
    theme(text = element_text(size=20))+
    labs(color="Estimators")
  return(p)
}


