#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

my$infile="in.smap";
my$outfile="$infile.converted.smap";
my$mode=0;# by default, convert a keimbahn sample- including mother/father info
# setting it to 1 will make a tumor sample conversion
#my$strict=0; # set to 1 to only accept 40x40 alignment quality

GetOptions('infile|i=s' => \$infile,'outfile|o=s'=> \$outfile,'mode|m=i' => \$mode) or warn "Using default parameters: $0 --from NAME\n";

chomp($infile,$outfile,$mode);

# in
open(IN,$infile)||die "$!";
my@newin = <IN>;
# TODO: add checks for infile , cleanup, adde xtra option for leaving commented lines in ,
# out
open(OUT,">",$outfile)|| die "$!";
#print ND "coordinates\tstrand\tsampleid\tunique_counts\tscore\tscore\tRefseqID\n";


# makeover: get the names of the columns into a array, resort them according to target array, print line indices accordingly
# target: we do not need different modes, infile will be automatically converted into one that works

# missing: Found_in_self_BSSSI_molecules#
# for the type: we need only the Type one as we triple it, for later printing we call them type type1 and type2
          # indices 1  2   3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
my@target_colnames=("SmapEntryID","QryContigID","RefcontigID1","RefcontigID2","QryStartPos","QryEndPos","RefStartPos","RefEndPos","Confidence","Type","Type","Type","XmapID1","XmapID2","LinkID","QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx","Zygosity","Genotype","GenotypeGroup","RawConfidence","RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter","Sample","Algorithm","Size","Present_in_%_of_BNG_control_samples","Present_in_%_of_BNG_control_samples_with_the_same_enzyme","Fail_BSPQI_assembly_chimeric_score","Fail_BSSSI_assembly_chimeric_score","OverlapGenes","NearestNonOverlapGene","NearestNonOverlapGeneDistance","Found_in_parents_assemblies","Found_in_parents_molecules","Found_in_self_BSPQI_molecules","Mother_molecule_count","Father_molecule_count","Found_in_self_molecules");



foreach my $line (@newin){
    # ignore lines with "#"
    if(!($line=~/#/)){
        chomp $line;
        my@parts=split(/\t/,$line);	# split line to find coordinates of gene
        # main functions here




    }
    elsif($line=~/#h /){
      # header with actual column names included
      my@header_in_file=split("\t",$line);

      # check h ere what header we found and which we did not find. we should be able to find all?

    }

}
