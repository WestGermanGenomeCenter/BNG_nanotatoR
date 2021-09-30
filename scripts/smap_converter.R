# parse a bionano .smap file and press it into the needed output format for nanotatoR
# to get the header:
# snakemake database_add_test_5.log --cores 2 --latency-wait 600
#   head    raw_keimbahn_smap_infile.smap | grep "^#h" | sed 's/#h//'|>header_line.smap
requireNamespace("dplyr")
library("dplyr")
'%!in%' <- function(x,y)!('%in%'(x,y))


args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least two arguments must be supplied (input file) (input_header file)  (output file)", call.=FALSE)
}

print("reading in files...")
#raw_input_file=read.csv("/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/input/raw_keimbahn_smap_infile.smap",header=F,comment.char = "#",sep="\t",row.names = NULL)
raw_input_file=read.csv(args[1],header=F,comment.char = "#",sep="\t",row.names = NULL)

#head(raw_input_file)
#raw_input_file_header=read.csv("/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/input/header_line.smap",header=T,sep="\t")
raw_input_file_header=read.csv(args[2],header=T,sep="\t",row.names=NULL)
#raw_input_file_header=read.csv(args[1],header=T,sep="\t",row.names = NULL)

#head(raw_input_file_header,n=35)
colnames(raw_input_file)=colnames(raw_input_file_header)

super_df=raw_input_file
colnames(super_df)=colnames(raw_input_file_header)

#head(super_df)
print(colnames(raw_input_file_header))

#needed_header_smap=read.csv("/home/daric/work_WGGC/SV_annotation/sv_nanotation_snake/data/needed_header.smap",header=T,sep="\t",row.names = NULL)


#colnames(needed_header_smap)

# R does convert "%" in a colname into a "." 

header_with_family_info=c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2","QryStartPos","QryEndPos","RefStartPos","RefEndPos","Confidence","Type","Type1","Type2","XmapID1","XmapID2","LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx","Zygosity","Genotype","GenotypeGroup","RawConfidence","RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter","Sample","Algorithm","Size","Present_in_._of_BNG_control_samples","Present_in_._of_BNG_control_samples_with_the_same_enzyme","Fail_BSPQI_assembly_chimeric_score","Fail_BSSSI_assembly_chimeric_score","OverlapGenes","NearestNonOverlapGene","NearestNonOverlapGeneDistance","Found_in_parents_assemblies","Found_in_parents_molecules","Found_in_self_BSPQI_molecules","Mother_molecule_count","Father_molecule_count","Found_in_self_molecules")
header_without_family_info=c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2","QryStartPos","QryEndPos","RefStartPos","RefEndPos","Confidence","Type","Type1","Type2","XmapID1","XmapID2","LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx","Zygosity","Genotype","GenotypeGroup","RawConfidence","RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter","Sample","Algorithm","Size","Present_in_._of_BNG_control_samples","Present_in_._of_BNG_control_samples_with_the_same_enzyme","Fail_BSPQI_assembly_chimeric_score","Fail_BSSSI_assembly_chimeric_score","OverlapGenes","NearestNonOverlapGene","NearestNonOverlapGeneDistance","Found_in_self_molecules")
header_rvp_file=c("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2","QryStartPos","QryEndPos","RefStartPos","RefEndPos","Confidence","Type","Type1","Type2","XmapID1","XmapID2","LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx","RawConfidence","RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter","Sample","Algorithm","Size","Present_in_._of_BNG_control_samples","Present_in_._of_BNG_control_samples_with_the_same_enzyme","Fail_BSPQI_assembly_chimeric_score","Fail_BSSSI_assembly_chimeric_score","OverlapGenes","NearestNonOverlapGene","NearestNonOverlapGeneDistance","PutativeGeneFusion","Found_in_self_molecules","Self_molecule_count")


# find out what we want- we can check for cols in the infile 

if("Found_in_parents_assemblies" %in% colnames(super_df)){
  needed_header_smap = header_with_family_info 
  print("Detected germ line input .smap with parents info, converting...")
}
if (("Genotype" %in% colnames(super_df))&&("Found_in_parents_assemblies" %!in% colnames(super_df))){
  needed_header_smap = header_without_family_info
  print("Detected de novo input .smap without parents info, converting...")
  
}
if (("PutativeGeneFusion" %in% colnames(super_df))&&("Found_in_parents_assemblies" %!in% colnames(super_df))){
  needed_header_smap = header_rvp_file
  print("Detected RVP pipeline input .smap, converting...")
  
}

#needed_header_smap=read.csv(args[3],header=T,sep="\t",row.names = NULL)

# check how much we need to synthesize
deletable_cols=colnames(raw_input_file)[colnames(raw_input_file) %!in% needed_header_smap]
usable_cols=intersect(needed_header_smap,colnames(raw_input_file))
needed_cols=needed_header_smap[needed_header_smap %!in% usable_cols]
#needed_cols
# rounding refstartpos , refend pos and others


# this time its: Type1 Type2 "Fail_BSPQI_assembly_chimeric_score" "Fail_BSSSI_assembly_chimeric_score" "Found_in_self_BSPQI_molecules"
cols_reconstructable=c("Type1","Type2","Fail_BSPQI_assembly_chimeric_score","Fail_BSSSI_assembly_chimeric_score","Found_in_self_BSPQI_molecules")
can_be_reconstructed=intersect(needed_cols,cols_reconstructable)
not_reconstructable=needed_cols[needed_cols %!in% can_be_reconstructed]
#print("difference: can be recon:")
#print(can_be_reconstructed)
#print("reconstructable")
#print(cols_reconstructable)
# needs to be built out: some information might change and we might not need to rebuild everyting here
print("trying to reconstruct...")
if(length(can_be_reconstructed)>= 0){ # if any hits, recreate whatever can be made
  # rebuild all needed ones
  #Type1=select(raw_input_file,c("Type"))
  #Type2=select(raw_input_file,c("Type"))
  print("reconstructing...")
  Type1=super_df$Type
  Type2=Type1
  

  Fail_BSSSI_assembly_chimeric_score= rep("not_applicable",nrow(raw_input_file))
  Fail_BSPQI_assembly_chimeric_score= rep("not_applicable",nrow(raw_input_file))
  Found_in_self_BSPQI_molecules=rep("yes",nrow(raw_input_file))
  
  #Size=select(raw_input_file,c("SVsize"))
  Size=super_df$SVsize
  if(length(Size)==0){
    Size=super_df$Size
    # in assemblys with parental info, the size column is different
  }
  
  print("check dims of columns:")
  print(paste("raw_input:",nrow(raw_input_file),"Type2",length(Type2),"Fail_BSSSI_assembly_chimeric_score",length(Fail_BSSSI_assembly_chimeric_score),"size:",length(Size)))
  
  
  unordered_new_df=cbind(raw_input_file,Type1,Type2,Fail_BSSSI_assembly_chimeric_score,Found_in_self_BSPQI_molecules,Fail_BSPQI_assembly_chimeric_score,Size)
  }


df_new=select(unordered_new_df,all_of(needed_header_smap))
#print("new frame colanmes:")
#print(colnames(df_new))
#print("dims of df new:")
#print(dim(df_new))
#print("number of colnames taken:")
#print(length(needed_header_smap))


# re-name the header accordingly
colnames(df_new)[1]="# SmapEntryID"

# create following output:
print("Columns in infile that will be deleted:")
print(deletable_cols)
print("Columns that will be added:")
print(needed_cols)
print("from these the following Columns can be reconstructed:")
print(can_be_reconstructed)
print("Dummy columns that will be added:")
print(not_reconstructable)

if(length(not_reconstructable)>0){
  # aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
  for (column_re in not_reconstructable){
    column_re =rep("-",nrows(raw_input_file))
    df_new=cbind(df_new,column_re)
    # rename the column accordingly
    colnames(df_new)[length(colnames(df_new))]<-column_re
    print(paste("reconstructed column",column_re))
  }
  # use the name as the name, add "-" or something to adjust
}



write.table(df_new,file=args[3],row.names = F,sep="\t")

