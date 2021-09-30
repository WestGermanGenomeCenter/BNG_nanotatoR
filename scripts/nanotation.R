

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("At least three arguments must be supplied (data folder path) (input file) (output file) (own_databale_file)", call.=FALSE)
}



# test()
#args[1]="/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/data"
#args[2]="/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/input/test_convert_1_out.smap"
#args[3]="/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/annotated/test1_test_converted_1_annotated.tsv"



packages_needed=c("BiocManager","nanotatoR","dplyr","stringr")
package.check<-lapply(packages_needed,
                      FUN=function(x){
                        if(!require(x,character.only=T)){
                          install.packages(x,dependencies=T)
                          library(x,character.only = T)
                        }
                      })

requireNamespace("nanotatoR","dplyr")


'%!in%' <- function(x,y)!('%in%'(x,y))

# params for each run of the fun
bedfile=paste0(args[1],"/bioMart_hg38_annotation.bed")
smap=args[2]
#bedfile="/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/data/bioMart_hg38_annotation.bed"
#smap="/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/input/test_convert_1_out.smap"
               
newName="Genes"
n_to_rename=1




# first annotation full to get all data
datcomp_annot=compSmapbed(smap=args[2],
                                     bed=paste0(args[1],"/bioMart_hg38_annotation.bed"),
                                     inputfmtBed="BED",outpath = ".",n = 10, returnMethod_bedcomp = c("dataFrame" ))



bedOverlap_RenameCols_returnCols <-function(bedfile,n_to_take,newName,n_hits,colsToTake){
  fullOutDF=compSmapbed(smap=smap,bed=bedfile,inputfmtBed="BED",outpath = ".",n = n_hits, returnMethod_bedcomp = c("dataFrame" ))
  # rename
  coln_old=colnames(fullOutDF)
  colnames_to_edit=tail(coln_old,n=3) # always rename the last 3 cols
  colnames_new=paste0(newName,colnames_to_edit)
  cols_complete=c(coln_old[1:(length(coln_old)-3)],colnames_new)
  colnames(fullOutDF)=cols_complete
  df_to_return=fullOutDF[,names(fullOutDF) %in% colsToTake]
  return(df_to_return)
  
}

bnbedOverlap_RenameCols_returnCols <-function(bedfile,n_to_take,newName,n_hits,colsToTake){
  bedfile=ncbi_refseq
  fullOutDF=compSmapbed(smap=smap,bed=bedfile,inputfmtBed="BNBED",outpath = ".",n = n_hits, returnMethod_bedcomp = c("dataFrame" ))
  # rename
  coln_old=colnames(fullOutDF)
  colnames_to_edit=tail(coln_old,n=3) # always rename the last 3 cols
  colnames_new=paste0(newName,colnames_to_edit)
  cols_complete=c(coln_old[1:(length(coln_old)-3)],colnames_new)
  colnames(fullOutDF)=cols_complete
  to_return_cols=fullOutDF[,((length(colnames(fullOutDF)))-2):(length(colnames(fullOutDF)))]
  return(to_return_cols)
  
}



############ adding correct coordinates, UCSC clickable link ################
full_sv_coords=paste0("chr",datcomp_annot$RefcontigID1,":",as.integer(datcomp_annot$RefStartPos),"-",as.integer(datcomp_annot$RefEndPos))
print("gene annotation finished...")
# now add clickable link
# the session cookie can be changed at any time, as well as the session settings
site_start="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position="
site_end="&hgsid=856936115_4aA1WK7ohwZycxkESPOCzhHUjzEm"
full_link=paste0(site_start,full_sv_coords,site_end)

sv_size=as.integer(datcomp_annot$Size)
svsize_clean=sv_size

# genomic features annotation
print("genomic features overlap...")
genomi_feats=paste0(args[1],"/genomic_feats_hg38_full.bed")
genomic_feats_annot=bedOverlap_RenameCols_returnCols(genomi_feats,1,"GenFeats",10,"GenFeatsOverlapGenes_strand_perc")


ncbi_refseq=paste0(args[1],"/nochr_hg38_annot_file.bed")
ncbi_refseq_annot=bnbedOverlap_RenameCols_returnCols(ncbi_refseq,1,"NCBIGenFeats",10,"NCBIGenFeatsOverlapGenes_strand_perc")



# databases

print("SV databases overlap...")

# 1000 genomes project variants
onekg_svs=paste0(args[1],"/1kg_sv_translated.bed")
onekg_sv_ovlp=bedOverlap_RenameCols_returnCols(onekg_svs,1,"thousandGenomes",10,"thousandGenomesOverlapGenes_strand_perc")

# dbVar
dbvar=paste0(args[1],"/dbvar_vars_hg38.bed")
dbvar_ovrlp=bedOverlap_RenameCols_returnCols(dbvar,1,"dbvar",10,"dbvarOverlapGenes_strand_perc")

# hgsvc, freeze  3 data
hgsvc_data=paste0(args[1],"/even_more_cleanhgsv2_svs_hg38.bed")
hgsvc_ovlp=bedOverlap_RenameCols_returnCols(hgsvc_data,1,"hgsvc_SVs",10,"hgsvc_SVsOverlapGenes_strand_perc")
#hgsvc_ovlp

# dgv, another database of genomic varians for humans
dgv_data=paste0(args[1],"/hg38_dgv_clean_try4.bed")
dgv_ovrlp=bedOverlap_RenameCols_returnCols(dgv_data,1,"dgv_SVs_",10,"dgv_SVs_OverlapGenes_strand_perc")


# own database overlap
own_bed=args[4]
own_overlap=bedOverlap_RenameCols_returnCols(own_bed,1,"own_cohort_",10,"own_cohort_OverlapGenes_strand_perc")



print("cleanup...")

########### cleanup functions to be implemented later
cleanup_overlap_column<-function(col_with_all_hits,min_ovrl=1,max_orvl=101){
  results=rep("",length(col_with_all_hits))
  to_fix=col_with_all_hits
  for (k in 1:length(to_fix)){
    line_complete=to_fix[k]
    single_entry=""
    if(grepl(";",line_complete,fixed=TRUE)){ # only useful if we have more than 1 hit
      all_singles=unlist(strsplit(line_complete,"[;]"))
      for (en in 1:length(all_singles)){
        single_entry=all_singles[en]
        entry_parts=unlist(strsplit(single_entry,"[:)]"))
        dist_kb=entry_parts[2]
        if((typeof(dist_kb)=="double")||(typeof(dist_kb)=="integer")||(typeof(dist_kb)=="character")){
          if(!(is.na(dist_kb))){
            number_io=suppressWarnings(round(as.numeric(as.character(dist_kb)),digits = 0))
            if(!(is.na(number_io))){
              if(number_io>min_ovrl){
                if(number_io<max_orvl){
                  results[k]<-paste0(results[k],single_entry,"_;")
                }
              }
            }
          }
        }
      }
    }
    else{
      single_entry=line_complete
      entry_parts=unlist(strsplit(single_entry,"[:)]"))
      dist_kb=entry_parts[2]
      if((typeof(dist_kb)=="double")||(typeof(dist_kb)=="integer")||(typeof(dist_kb)=="character")){
        if(!(is.na(dist_kb))){
          number_io=suppressWarnings(round(as.numeric(as.character(dist_kb)),digits = 0))
          if(!(is.na(number_io))){
            # print(number_io)
            if(number_io>min_ovrl){
              if(number_io<max_orvl){
                results[k]<-paste0(results[k],single_entry,"_;")
              }
            }
          }
        }
      }
    }
  }
  return (results)
}

####
dgv_cleanup<-function(dgv_ovrlp){
  svsize_clean=svsize_clean
  to_cleanup=dgv_ovrlp
  svtype_searching=datcomp_annot$Type
  library(stringr)
  col_with_good_hits=c(rep("",length(to_cleanup)))# to be filled with good overlaps
  for(i in seq(1,length(to_cleanup),1)){# for each row
    full_line=to_cleanup[i]
    needed_type=svtype_searching[i]
    ned_split=unlist(strsplit(needed_type, "[_]"))
    part_underst=ned_split[1]
    needed_svsize=as.integer(svsize_clean[i]/2)
    all_parts_of_entry=unlist(strsplit(full_line, "[;]"))
    for (o in seq(1,length(all_parts_of_entry),1)){
      part_curr=all_parts_of_entry[o]
      potential_overlap_perc=5000000# placeholder if no number gets found
      parts_of_entry=unlist(strsplit(part_curr,"[:)]"))
      potential_overlap_perc=parts_of_entry[2]
      new_splitentyr=unlist(strsplit(part_curr,"[_]"))
      hit_svsize=suppressWarnings(as.integer(new_splitentyr[2]))
      if(!is.na(suppressWarnings(as.integer(potential_overlap_perc)))){
        if((as.integer(potential_overlap_perc)) <= 100){ # at some point this need s to be <=100
          if((as.integer(potential_overlap_perc)) >= 50){
            clean_entry=str_replace(part_curr,"essv.*\\(","") # cleanup mor, than use the clean entry instead, also clean by svtype overlap
            if( grepl( part_underst, part_curr, fixed = TRUE)){
              if(hit_svsize >=needed_svsize){
                more_splitted=unlist(strsplit(parts_of_entry[1],"[-]"))
                part_w_n=unlist(strsplit(more_splitted[2],"[_]"))
                freq_str="freq=NA"
                if(!is.na(suppressWarnings(as.integer(part_w_n[2])))){
                  if(!is.na(suppressWarnings(as.integer(part_w_n[3])))){
                    if(!is.na(suppressWarnings(as.integer(part_w_n[4])))){
                      n_samples=as.numeric(as.character(part_w_n[2]))
                      n_max=max(as.numeric(as.character(c(part_w_n[3],part_w_n[4]))))
                      freq_str=paste0("freq=",round((n_max/n_samples),digits=5))
                    }
                  }
                }
                col_with_good_hits[i]<-paste0(col_with_good_hits[i],part_curr,freq_str,"_;")
              }
            }
          }
        }
      }
    }
    
  }
  dgv_hits=col_with_good_hits
  return(dgv_hits)
}
### 


cleanup_overlap_column_sadsplit_with_length<-function(col_with_all_hits,min_ovrl=1,max_orvl=101){
  results=rep("",length(col_with_all_hits))
  to_fix=col_with_all_hits
  for (k in 1:length(to_fix)){
    line_complete=to_fix[k]
    single_entry=""
    svsize_needed=svsize_clean[k]/2
    if(grepl(";",line_complete,fixed=TRUE)){ # only useful if we have more than 1 hit
      all_singles=unlist(strsplit(line_complete,");"))
      for (en in 1:length(all_singles)){
        single_entry=all_singles[en]
        entry_parts=unlist(strsplit(single_entry,"[:)]"))
        dist_kb=entry_parts[2]
        if((typeof(dist_kb)=="double")||(typeof(dist_kb)=="integer")||(typeof(dist_kb)=="character")){
          if(!(is.na(dist_kb))){
            number_io=suppressWarnings(round(as.numeric(as.character(dist_kb)),digits = 0))
            if(!(is.na(number_io))){
              if(number_io>min_ovrl){
                if(number_io<max_orvl){
                  parts_full=unlist(strsplit(single_entry,"[_]"))
                  size=suppressWarnings(round(as.numeric(as.character(parts_full[6])),digits = 0))
                  if(size>=svsize_needed){
                    results[k]<-paste0(results[k],single_entry,"_;")
                  }
                }
              }
            }
          }
        }
        
      }
    }
    else{
      single_entry=line_complete
      entry_parts=unlist(strsplit(single_entry,"[:)]"))
      dist_kb=entry_parts[2]
      if((typeof(dist_kb)=="double")||(typeof(dist_kb)=="integer")||(typeof(dist_kb)=="character")){
        if(!(is.na(dist_kb))){
          number_io=suppressWarnings(round(as.numeric(as.character(dist_kb)),digits = 0))
          if(!(is.na(number_io))){
            if(number_io>min_ovrl){
              if(number_io<max_orvl){
                parts_full=unlist(strsplit(single_entry,"[_]"))
                size=suppressWarnings(round(as.numeric(as.character(parts_full[6])),digits = 0))
                if(size>=svsize_needed){
                  results[k]<-paste0(results[k],single_entry,"_;")
                }
              }
            }
          }
        }
      }
    }
  }
  return (results)
}

##
database_hit_cleanup<-function(col_with_hgsvc_hits){
  hgsv_overlap_cleanup=col_with_hgsvc_hits
  col_with_good_hits=c(rep("",length(hgsv_overlap_cleanup)))# to be filled with good overlaps
  for(i in seq(1,length(hgsv_overlap_cleanup),1)){# for each row
    full_line=hgsv_overlap_cleanup[i]
    needed_type=datcomp_annot$Type[i]
    ned_split=unlist(strsplit(needed_type, "[_]"))
    part_underst=ned_split[1]
    needed_svsize=as.integer(svsize_clean[i]/2)
    all_parts_of_entry=unlist(strsplit(full_line, "[;]"))
    for (o in seq(1,length(all_parts_of_entry),1)){
      part_curr=all_parts_of_entry[o]
      potential_overlap_perc=5000000# placeholder if no number gets found
      parts_of_entry=unlist(strsplit(part_curr,"[:)]"))
      potential_overlap_perc=parts_of_entry[2]
      new_splitentyr=unlist(strsplit(part_curr,"[_=]"))
      hit_svsize=suppressWarnings(as.integer(new_splitentyr[4]))
      if(!is.na(suppressWarnings(as.integer(potential_overlap_perc)))){
        if((as.integer(potential_overlap_perc)) <= 100){ # at some point this need s to be <=100
          if((as.integer(potential_overlap_perc)) >= 50){
            if( grepl( needed_type, part_curr, fixed = TRUE)){
              if(hit_svsize >=needed_svsize){
                col_with_good_hits[i]<-paste0(col_with_good_hits[i],part_curr,"_;")
              }
            }
          }
        }
      }
    }
  }
  SV_hgsv_overlp=col_with_good_hits
  return(SV_hgsv_overlp)
}


# tzryout
#head(hgsvc_ovlp)
hgsvc_overlap=database_hit_cleanup(hgsvc_ovlp)
#head(hgsvc_overlap)
dbvar_overlap=cleanup_overlap_column_sadsplit_with_length(dbvar_ovrlp)

# dgv needs also a cleanup
# head(dgv_ovrlp)
DGV_overlap=dgv_cleanup(dgv_ovrlp)




gene_annotation_cleanup<-function(col_with_genes){
  # cleanup with the genes of 0 overlap
  cleanup_singletons_zero=gsub("[:alnum:]{2,17}[:punct:]{4}0[:punct:]{1}","",col_with_genes)
  return(cleanup_singletons_zero)
}


internal_overlap_cleanup<-function(vec_w_own_overlaps){
  # cleanup stuff like: +(\"insertion\":100)
  cleaner_col=gsub("+(deletion:NA)","",vec_w_own_overlaps,fixed=T)
  cleaner_col=gsub("+(inversion:NA)","",cleaner_col,fixed=T)
  cleaner_col=gsub("+(duplication:NA)","",cleaner_col,fixed=T)
  cleaner_col=gsub("+(insertion:NA)","",cleaner_col,fixed=T)
  cleaner_col=gsub("+(translocation:NA)","",cleaner_col,fixed=T)
  return(cleaner_col)
}
ovrlp_internal_cohort=internal_overlap_cleanup(own_overlap)
overlap_internal_clean=cleanup_overlap_column(own_overlap,max_orvl = 99)
# cleanup of other columns

#genomic_feats_clean=genomic_feats_cleanup(genomic_feats_annot)
annotations_clean=gene_annotation_cleanup(ncbi_refseq_annot$NCBIGenFeatsOverlapGenes_strand_perc)



print("creating output file...")


# cleanup every column by one of thiese functions
full_df=cbind(datcomp_annot,sv_size,genomic_feats_annot,full_sv_coords,ncbi_refseq_annot,hgsvc_overlap,onekg_sv_ovlp,dbvar_overlap,DGV_overlap,overlap_internal_clean,full_link)
print("writing full output file... ")
write.table(full_df,args[3],row.names = F, sep="\t")
# preparation for sub-set dataframe
# get all colnames from datcomp_annot
original_smap_cols=colnames(datcomp_annot)
# output is : 
#[1] "SmapEntryID"                                              "QryContigID"                                              "RefcontigID1"                                            
#[4] "RefcontigID2"                                             "QryStartPos"                                              "QryEndPos"                                               
#[7] "RefStartPos"                                              "RefEndPos"                                                "Confidence"                                              
#[10] "Type"                                                     "Type1"                                                    "Type2"                                                   
#[13] "XmapID1"                                                  "XmapID2"                                                  "LinkID"                                                  
#[16] "QryStartIdx"                                              "QryEndIdx"                                                "RefStartIdx"                                             
#[19] "RefEndIdx"                                                "RawConfidence"                                            "RawConfidenceLeft"                                       
#[22] "RawConfidenceRight"                                       "RawConfidenceCenter"                                      "Sample"                                                  
#[25] "Algorithm"                                                "Size"                                                     "Present_in_._of_BNG_control_samples"                     
#[28] "Present_in_._of_BNG_control_samples_with_the_same_enzyme" "Fail_BSPQI_assembly_chimeric_score"                       "Fail_BSSSI_assembly_chimeric_score"                      
#[31] "OverlapGenes"                                             "NearestNonOverlapGene"                                    "NearestNonOverlapGeneDistance"                           
#[34] "PutativeGeneFusion"                                       "Found_in_self_molecules"                                  "Self_molecule_count"                                     
#[37] "OverlapGenes_strand_perc"                                 "Upstream_nonOverlapGenes_dist_kb"                         "Downstream_nonOverlapGenes_dist_kb"                      


print("done.")