#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# convert smap files to add to a steadily growing .bed file with custom cohort creation

my$linfile= "../annotated/smap_annotated.csv";
my$database="../data/current_cohort.bed";


GetOptions('database|d=s' => \$database,'infile|i=s' => \$linfile) or warn "Using default parameters: $0 --from NAME\n";

chomp $linfile;
chomp $database;
my%svs_coords_hash=();
my@results;

# to include the samplename in the SV for each one, we need to extract that from the file, then simply add this to the type string
my$sample_name_from_file="unkn"; # default if no samplename is found
$linfile =~/[0-9]{2,13}/;
$sample_name_from_file = $&;
print("attaching SVs from sample $sample_name_from_file into $database\n");

# read infile into hash
open(DB,$database)|| die "$!";
my@database_lines= <DB>;
print ("checking the current database entries..\n");
# database is expected to be proper .bed format
for(my$o=0;$o<scalar(@database_lines);$o++){
  my$line_db=$database_lines[$o];
  chomp $line_db;
  if(!($line_db=~/#/ || $line_db=~/SmapEntryID/)){
    # no header line
      my @all_lineparts=split("\t",$line_db);
      # we need to build the hash key using 3 parts- chrom start end.
      my$chr=$all_lineparts[0];
      my$start=$all_lineparts[1];
      my$end=$all_lineparts[2];
      my$strand=$all_lineparts[4];
      # otherwise nanotatoR will spit out an error
      if($end < $start){
        next;
      }
      my$seen_coords="$chr"."_"."$start"."_"."$end"."_"."$strand";# to be later translated into .bed
      # enter all into the hash

      #print("old file: $seen_coords key: value : $all_lineparts[3]\n");
      $svs_coords_hash{$seen_coords}=$all_lineparts[3]; # saving the sample name as value here because why not
  }
}
close(DB);
print("adding the entries from $linfile\n");
########## now read the new smap into the same hash if we have new stuff
open(IN,$linfile)|| die "$!";
my@allelines= <IN>;
my$sv_stelle;
# default  on rvp files
foreach my $line (@allelines){
  chomp $line;
  if($line=~/[a-z]/){# check if there is not only whitespace in a row
  # find the index of Size or SVsize in the header

    if($line=~/SmapEntryID/ig){
    # split intop parts
      my@header_line_parts=split("\t",$line);
      $sv_stelle=0;
      # new way to find it
      my$nummer_grad=0;
      foreach my $header_part (@header_line_parts){

        if($header_part=~/Size/){
          $sv_stelle=$nummer_grad;
        }
        $nummer_grad++;
      }

    #  print("found out size header index: $sv_stelle\n");
      # check later, maybe we need +1 here
    }


  if(!($line=~/#/ || $line=~/SmapEntryID/ || $line=~/RefStart/)){ # check for header here- we do not need this
    # now convert smap annotated coords to .bed with additional info in 4th column
    my@line_smap_parts=split("\t",$line);
    my$refchr=int($line_smap_parts[2]);
    my$start=int($line_smap_parts[6]);
    my$end=int($line_smap_parts[7]);
    my$stra="+";# dummy since proper .smap apparently does not include strand info
    # check if end coordinates are bigger then start
  #  print("index with SVsize in smap: $sv_stelle\n");
  #  print("size id here: $line_smap_parts[$sv_stelle]\n");
    my$sv_size=int($line_smap_parts[$sv_stelle]);
    if($end < $start){
      #my$start_new= $start;
      my$end_new=$end;
      $end=$start;
      $start =$end_new;
    }
    if($start eq "-1"){
      $start = $end -1;
    }
    my$coords_key="$refchr"."_"."$start"."_"."$end"."_"."$stra";
    if(!(exists $svs_coords_hash{$coords_key})){
      # if not already in DB, add it
      $sv_size=~s/\-1//;
    #  print(" new file: now getting the hash value of $coords_key to be $line_smap_parts[9]"."_"."$sv_size\n");

      $svs_coords_hash{$coords_key}="$line_smap_parts[9]"."_"."$sv_size"."bp_"."sample="."$sample_name_from_file";# the type of the SV is here the 4th field and _ and size of the SV
    }
  }
}
}
# assemble all information


my@coordinates_complete=keys(%svs_coords_hash);
open(DB,">$database")|| die "$!";

foreach my$line_comp (@coordinates_complete){
  # get all values
  # split by _
  # make tab-separated entries with strand
  my@parts_of_coords=split("_",$line_comp);
  # change the structure: chr start end info_with_sample_ strand
  #$parts_of_coords=~s/"//g;
  #my$line_bed_nostrand=join("\t",@parts_of_coords);
  #print("now adding to db file: $parts_of_coords[0]\t$parts_of_coords[1]\t$parts_of_coords[2]\t$svs_coords_hash{$line_comp}\t$parts_of_coords[3]\n");
  my$final_bed_line="$parts_of_coords[0]\t$parts_of_coords[1]\t$parts_of_coords[2]\t$svs_coords_hash{$line_comp}\t$parts_of_coords[3]";
  #my$final_bed_line="$line_bed_nostrand"."\t"."$svs_coords_hash{$line_comp}";
  $final_bed_line=~s/"//g;


  # we expect: chr start end info_type_length strand (+)
  # sanity checks of final line
  if(($final_bed_line=~s/\+\t+[A-z]//)||($final_bed_line=~s/\-1\t+//)||($final_bed_line=~s/inversion_partial//)){
  #  print "skipping $final_bed_line because of line error!\n";
    next;
  }


  # only clean lines enter the DB file
  if($final_bed_line=~/[a-z]/){# do not produce emtpy lines
    print DB "$final_bed_line\n";
  }
}
print("done.\n");
1;
