#!/usr/bin/perl -w
use strict;
use lib '/ifs5/PC_PA_AP/USER/lixf/Package';
use GFF;
unless (@ARGV == 5) {
	die <<End;

Description:
    This script is used to add location information to each ortholog pair.

Usage:
    perl $0 xyz.ort1 ref_gff hit_gff ref_species hit_species >xyz.ort2

End
}
my $in_file = shift;
my $gff_file1 = shift;
my $gff_file2 = shift;
my $specise1 = shift;
my $specise2 = shift;
my %species = ($specise1=>1, $specise2=>1);

my (%chrGene1, %geneIndex1);
getChrGene($gff_file1, \%chrGene1, \%geneIndex1);
my (%chrGene2, %geneIndex2);
getChrGene($gff_file2, \%chrGene2, \%geneIndex2);

my %result;
open IN, $in_file;
while (<IN>) {
	my ($g1, $g2, $level) = (split /\s+/)[0,1,-1];
	my ($chr1, $bg1, $ed1, $strand1, $index1) = @{$geneIndex1{$g1}};
	die "$g2" unless $geneIndex2{$g2};
	my ($chr2, $bg2, $ed2, $strand2, $index2) = @{$geneIndex2{$g2}};
	$result{$chr1}{$index1} = "$g1\t$chr1\t$index1\t$strand1\t$bg1\t$ed1\t$g2\t$chr2\t$index2\t$strand2\t$bg2\t$ed2\t$level\n";
}
close IN;

foreach my $chr1 (keys %chrGene1) {
	foreach my $p (@{$chrGene1{$chr1}}) {
		my ($g1, $bg1, $ed1, $strand1) = @$p;
		my $index1 = $geneIndex1{$g1}->[-1];
		next if $result{$chr1}{$index1};
		$result{$chr1}{$index1} = "$g1\t$chr1\t$index1\t$strand1\t$bg1\t$ed1\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
}

foreach my $chr (sort keys %result) {
	foreach my $index (sort {$a <=> $b} keys %{$result{$chr}}) {
		print $result{$chr}{$index};
	}
}
