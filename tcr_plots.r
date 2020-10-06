library(ggplot2)
library(dplyr)
library("stringr")
library(data.table)
library(scales)
library(Rcpp)
library(devtools)
library(ggrepel)
pkgbuild::has_rtools(TRUE)
Sys.setenv("PKG_CXXFLAGS" = "-fopenmp -std=c++14")
Sys.setenv("PKG_LIBS" = "-fopenmp")


## CPP function to calculate the Gini Index (much faster than lawstat)
cppFunction("
double my_gini(const NumericVector& gph){

size_t size=gph.size();

double avr=0, gini=0;
#pragma omp parallel for reduction(+:avr)
for(size_t i=0; i<size; ++i)
  avr+=gph[i];
avr/=double(size);

#pragma omp parallel for reduction(+:gini)
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
    res=lapply(lseq, function(x, t_ldataf){
      if(length(unique(t_ldataf$trxvseqtrxj==x))==2 ){
        return(t_ldataf[t_ldataf$trxvseqtrxj==x,])
      }
      else{
        return(data.frame(patient=t_ldataf$patient[1], frequency=0, trxvseqtrxj=as.character(x)))
      }
    }, t_ldataf=ldataf[[i]] )
    
    tddf=rbind(tddf, do.call(rbind, res))
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
  aseq=startdf$trxvseqtrxj
  
  for(i in 1:length(ldataf)){
    aseq=ldataf[[i]]$trxvseqtrxj[ldataf[[i]]$trxvseqtrxj %in% aseq]
  }
  if(length(aseq)>nseq){
    return(aseq[1:nseq])
  }
  else{
    return(aseq)
  }
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
      limits=c(0.00000001,0.5)
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
      limits=c(0.00000001,0.5)
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
      limits=c(0.00000001,0.5)
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
      limits=c(0.00000001,0.5)
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

## calculate the Gini-Simpson index
gini_simpson=function(adf){
  res=1.0 - sum(adf$frequency^2)
  return(res)
}


## we can use the morisita horn formula even using the probability
##  the formula is    2*sum_i_to_s Count_A_i * Count_B_i / ((sum_i_to_s f_Ai^2 * sum_i_to_s f_B_i^2)*total_A*total_B)


morisita_horn=function(databind, type1, type2){
  adf1=databind[databind$patient==type1,]
  adf2=databind[databind$patient==type2,]
  common_seq=intersect(adf1$trxvseqtrxj, adf2$trxvseqtrxj)
  adf1=adf1[ adf1$trxvseqtrxj %in% common_seq, ]
  adf2=adf2[ adf2$trxvseqtrxj %in% common_seq, ]
  total_n1=nrow(adf1)
  total_n2=nrow(adf2)
  # sorting the dataframe in the same order
  rownames(adf1)=common_seq
  rownames(adf2)=common_seq
  adf1=adf1[common_seq,]
  adf2=adf2[common_seq,]
  return(  2*sum(adf1$frequency*adf2$frequency)/( sum(adf1$frequency^2)+sum(adf2$frequency^2))  )
}


bhattacharyya=function(databind, type1, type2){
  adf1=databind[databind$patient==type1,]
  adf2=databind[databind$patient==type2,]
  common_seq=intersect(adf1$trxvseqtrxj, adf2$trxvseqtrxj)
  adf1=adf1[ adf1$trxvseqtrxj %in% common_seq, ]
  adf2=adf2[ adf2$trxvseqtrxj %in% common_seq, ]
  total_n1=nrow(adf1)
  total_n2=nrow(adf2)
  # sorting the dataframe in the same order
  rownames(adf1)=common_seq
  rownames(adf2)=common_seq
  adf1=adf1[common_seq,]
  adf2=adf2[common_seq,]
  return( sum(sqrt(adf1$frequency*adf2$frequency )) )
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

## convert David Barras dataframe to a list of dataframe that can be used in the violin plots
## each element of the list is a dataframe similar to the ones obtained in "cleanse"
convert_barras=function(abarr_list){
  lastb=abarr_list[[length(abarr_list)]]
  ## first we need to find the names of the various objects
  name_obj=colnames(lastb)[grepl("is_", colnames(lastb))]
  name_obj=sapply(name_obj, function(x){
    return(strsplit(x, "is_", fixed=TRUE)[[1]][2])})
  name_obj2=sapply(name_obj, function(x){
    return(strsplit(x, "in_", fixed=TRUE)[[1]][2])})
  
  l_df=list()
  cl=1
  ## create a number of dataframes 
  for(n in 1:length(name_obj)){
    cn1=paste0("is_",name_obj[n])
    tf_vect=lastb[,cn1]=="yes"
    cn2=paste0("Perc_", name_obj[n])
    
    pname=rep(name_obj2[n], sum(tf_vect))  
    tdf=data.frame(trxvseqtrxj=lastb[tf_vect,"unique_id"] , frequency=lastb[tf_vect,cn2]  , patient=pname  )
    l_df[[cl]]=tdf
    cl=cl+1
  }
  return(l_df)
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
  astatdf=setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("entropy", "clonality", "gini", "gini-simpson", "richness", "sampletype"))
  
  for(k in sampletype){
    adf=databind[databind$patient==k,]
    en=entropy(adf)
    cl=clonality(adf)
    gi=my_gini(adf$frequency)
    gs=gini_simpson(adf)
    rc=nrow(adf)
    adf=data.frame(en, cl, gi, gs, rc, k)
    astatdf=rbind(astatdf, adf)
  }
  
  rel_weight=weight_of_shared(databind, sampletype)
  
  astatdf=cbind(relative_weight=rel_weight$rel_fre, astatdf)
  
  colnames(astatdf)=c("Shared freq. (weighted)", "Entropy", "Clonality", "Gini", "Gini-Simpson", "Richness", "sampletype")
  astatdf2=melt(astatdf, id.vars="sampletype")
  p=ggplot(astatdf2, aes(x=sampletype, y=value, group=variable, color=variable))+
    geom_line()+
    geom_point()+
    geom_text_repel(aes(x=sampletype, y=value, label=signif(value,2)), show.legend = FALSE) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)) )+
    theme(text = element_text(size=20))+
    labs(color="Estimators")
  return(p)
}


