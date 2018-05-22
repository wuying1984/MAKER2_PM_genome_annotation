#!/usr/bin/perl
use warnings;
use strict;

#############################################################################
####################For get high quality gene list###########################
####################with standard start codon (M)############################
#############################################################################

#usage example perl 3get_hc_M.pl HC_gene.txt proteins.fa HC_gene_M.txt

my $sp_file = $ARGV[0]; #A gene list, first column is gene id
open R, $sp_file;
my @sp_lines = <R>;
close R;
chomp @sp_lines;
my %name2sp;
foreach my $line (@sp_lines){
        if ($line =~ /2nd/){ #need to change the symble here
                my ($gene,$d,$cut)=split /,/,$line;
                $gene=~ s/-RA//;
                $name2sp{$gene}="$d\t$cut";
        }
}

my $in_fasta= $ARGV[1]; #input the proteins 
open R, $in_fasta;
my @lines = <R>;
close R;
chomp @lines;
my %name2pro;
my $name = '';

foreach my $line (@lines){
        if ($line =~ />(.*) .* .* .* .*/){
                $name=$1;
                $name=~ s/-RA//;
#       print $name, "\n";
        }
        else {
                $name2pro{$name} = $line;
        }
}

my $out_fasta =$ARGV[2]; #output file name

open W, ">$out_fasta";
foreach my $gene (keys %name2sp){
        print $gene if !$name2pro{$gene};
        my $pro = $name2pro{$gene};
#       delete $name2sp{$gene} if $pro !~ /^M/;
        if ($pro =~ /^M/){
                my $len = length $pro;
                my (undef,$cut)=split /\t/,$name2sp{$gene};
                my $nosp= substr ($pro,$cut,($len-$cut));
                print W "$gene\n";

        }
}


#my %name2spnotm;
#foreach my $line (@tm_lines){
#       if ($line !~ /#/){
#               my @eles = split /\s+/,$line;
#               my $gene = $eles[0];
#               $gene =~ s/-RA//;
#               my $helix = $eles[4];
#               $helix =~ s/PredHel\=//;
#               print "$gene,$helix\n";
#               if ($helix ==0 && $name2sp{$gene} ){
#                       $name2spnotm{$gene}=$name2sp{$gene};
#               }
 #
#       }
#}

#open W, ">all.spnotmhmm.list";
#print W "gene\tdscore\tcut\n";
#foreach my $gene (keys %name2spnotm){
#       print W $gene, "\t", $name2spnotm{$gene}, "\n";
#}

#close W;
#

#open W, ">UCSC1.sp.list";
#print W "gene\tdscore\tcut\n";
#foreach my $gene (keys %name2sp){
#       print W $gene, "\t", $name2sp{$gene}, "\n";
#}
