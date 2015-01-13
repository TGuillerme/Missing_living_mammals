#!/usr/bin/perl -w

#this script converts a tree ouput from TNT to NEXUS (1 tree per file)
# Carlos PeÃ±a 2011-03-16

use strict;
use Bio::TreeIO;
        my $usage = "script.pl INFILE OUTFILE\n";
        my $infile = shift or die $usage;
        my $outfile = shift or die $usage;

# convert manually to newick format
open ( my $file, "<", $infile) or die $!;
my @data = <$file>;
my $tree = $data[1];
$tree =~ s/\s{1}\)/\)/g;
$tree =~ s/\s{1}/,/g;
$tree =~ s/\)\(/\),\(/g;
$tree =~ s/;//g;
$tree =~ s/,$//g;
close($file);
open ( $file, ">", $infile) or die $!;
print $file $tree;
close($file);

my ($filein,$fileout) = @ARGV;
my ($format,$oformat) = qw(newick nexus);
my $in = Bio::TreeIO->new(-file => $infile, -format => $format);
my $out= Bio::TreeIO->new(-format => $oformat, -file => ">$outfile");

while( my $t = $in->next_tree ) { 
        $out->write_tree($t);
}
