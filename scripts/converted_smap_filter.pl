#!/usr/bin/perl -w
use strict;

# filter smap files for 1% in BNG ctrls
# and split translocations into two entries
use Getopt::Long qw(GetOptions);
my$linfile= "converted_smap.smap";
my$to_replace="1";
my$min_confidence="0.8";

GetOptions('to_replace|r=s' => \$to_replace,'infile|i=s' => \$linfile,'confidence|c=s'=>\$min_confidence) or warn "Using default parameters: $0 --from NAME\n";

chomp $linfile;
chomp $to_replace;
chomp $min_confidence;
# positions(starting from 0.th position) of refcontigid1 and 2, needs to be known for translocation re-editing
my$chr_1_index=2;
my$chr_2_index=3;
my$start_base_index=4;
my$end_base_index=5;
my$ref_start_pos=6;
my$ref_end_pos=7;


my$assumed_index=29; # default in germline smap files
my$assumed_same_index=30; # for Present_in_%_of_BNG_control_samples_with_the_same_enzyme
my$asumed_type_index=9; # if we find a translocation, we need to print the line twice: one line for each breakpoint
my$assumed_confidence_filter=8;
#print "reading input file $linfile ...\n";
open(IN,$linfile)|| die "$!";
my@allelines= <IN>;
if(scalar(@allelines)<1){
  die "seems like the input file is empty, cancelling...\n";
}
foreach my $line (@allelines){
  # first find out what item in the list the %bionano is
  chomp $line;
  $line=~s/"//g;
  my@line_p=split(/\t/,$line);


  if($line =~/\#/){

    my @i = grep { $line_p[$_] eq 'Present_in_._of_BNG_control_samples' } 0 .. $#line_p;
    $assumed_index=$i[0];
    my @k = grep { $line_p[$_] eq 'Present_in_._of_BNG_control_samples_with_the_same_enzyme' } 0 .. $#line_p;
    $assumed_same_index=$k[0];
    my @o = grep { $line_p[$_] eq 'Type' } 0 .. $#line_p;
    $asumed_type_index=$o[0];

    my@c = grep { $line_p[$_] eq 'Confidence' } 0 .. $#line_p;
    $assumed_confidence_filter=$c[0];
    # cleanup the line
    $line=~s/Present_in_._of_BNG_control_samples_with_the_same_enzyme/Present_in_%_of_BNG_control_samples_with_the_same_enzyme/;
    $line=~s/Present_in_._of_BNG_control_samples/Present_in_%_of_BNG_control_samples/;

    print "$line\n";




  }
  else{
    if($line=~/[a-z]/){# check for empty line
    # isolate the % and type values from each row
    my$col_to_filter=$line_p[$assumed_index];
    my$type_ident=$line_p[$asumed_type_index];
    my$col2_to_filter=$line_p[$assumed_same_index];

    if($line_p[$assumed_confidence_filter]< $min_confidence){
      next; #  SV is not confident enough
    }

    # need to cleanup the inversion_partial things:
      #- check chrom, querypos : end needs to ne bigger than start + chrom cannot be negative!

    if(($line_p[$chr_2_index] < $line_p[$chr_1_index]) ||($line_p[$chr_2_index] eq "-1") ){
      $line_p[$chr_2_index] = $line_p[$chr_1_index] ;
    }
    if(($line_p[$end_base_index] < $line_p[$start_base_index])||($line_p[$end_base_index] eq "-1")){
      $line_p[$end_base_index]=$line_p[$start_base_index] ;
    }

    if(($line_p[$ref_end_pos] < $line_p[$ref_start_pos])||($line_p[$ref_end_pos] eq "-1")){
      $line_p[$ref_end_pos] = $line_p[$ref_start_pos] ;
    }

    if($line=~/translocation_interchr/){
      # even for translocations the freq needs to be filtered
      if(($col_to_filter <= $to_replace)&&($col2_to_filter <= $to_replace)){
        # now re_edit the breakpoints
        my$new_line=$line;

        # first, make refcontigid2= refcontigid1+1
        my@all_lineparts=split(/\t/,$new_line);
        my@second_line=@all_lineparts;
        $all_lineparts[$chr_2_index] = $all_lineparts[$chr_1_index];
        $second_line[$chr_1_index]= $second_line[$chr_2_index];
        # now we have edited the chromosomes- now the start position
        $all_lineparts[$end_base_index] = $all_lineparts[$start_base_index] +1;
        $second_line[$start_base_index]= $second_line[$end_base_index] -1;
        # adding all parts together again
        my$line_2=join("\t",@second_line);
        my$line_1=join("\t",@all_lineparts);

        print "$line_1\n";
        print "$line_2\n";

      }
    }
    elsif(($col_to_filter <= $to_replace)&&($col2_to_filter <= $to_replace)){
      # here we have the SVs that have a BNG freq below 1 % - with same enzyme and overall aswell.
      my$new_line=join("\t",@line_p);
      print "$new_line\n";

    }
  }
  }

}
1;
