#!/usr/bin/perl -w
use strict;
die "Usage: <IN> <gff dir>\n" unless @ARGV == 2;

my %chrGeneNum;
my $gff_dir = $ARGV[1];
opendir DH, $gff_dir;
while (my $file = readdir DH) {
	next unless $file =~ /\.gff$/;
	my $sp = (split /\./, $file)[0];
	my $in_file = "$gff_dir/$file";
	open IN, $in_file;
	while (<IN>) {
		my @info = split /\s+/;
		next unless $info[2] eq "mRNA";
		$chrGeneNum{$sp}{$info[0]} ++;
	}
	close IN;
}
closedir DH;

## restore the whole table
my @data;
open IN, $ARGV[0];
@data = <IN>;
close IN;
## creat psudo lines
my @tmp = split /\s+/, $data[0];
my $col_n = @tmp;
@tmp = ();
for (my $i = 0; $i < $col_n; $i ++) {
	push @tmp, 10e10;
}
my $tmp = join "\t", @tmp;
## add three pseudo lines to the head and the tail
for (my $i = 0; $i < 5; $i ++) {
	unshift @data, "$tmp\n";
	push @data, "$tmp\n";
}

for (my $i = 5; $i < @data - 5;  $i ++) {
	my @info = split /\s+/, $data[$i];	
	next if $info[7] eq "NA";
	#print "$data[$i]";
	my $up_check = 0;
	my $up_i;
	for (my $j = $i - 1; $j >= 0; $j --) {
		my @up = split /\s+/, $data[$j];
		next if $up[7] eq "NA";
		if ($up[1] eq $info[1] && $up[7] eq $info[7] && abs($up[8] - $info[8]) < 5) {
			$up_check = 1;
			$up_i = $j;
			last;
		} else {
			last;
		}
	}

	my $down_check = 0;
	my $down_i;
	for (my $j = $i + 1; $j < @data; $j ++) {
		my @down = split /\s+/, $data[$j];
		next if $down[7] eq "NA";
		if ($down[1] eq $info[1] && $down[7] eq $info[7] && abs($down[8] - $info[8]) < 5) {
			$down_check = 1;
			$down_i = $j;
			last;
		} else {
			last;
		}
	}
	
	my $chr_gene_num_check = 0;
	for (my $j = 0; $j < @info-1; $j += 6) {
		my $sp = (split /_/, $info[$j])[0];
		my $chr = $info[$j+1];
		$chr_gene_num_check = 1 if $chrGeneNum{$sp}{$chr} == 1;
	}
	
	if ($up_check && $down_check) {
		#print $data[$up_i], "#$data[$i]", $data[$down_i], "\n";
		print "$data[$i]";
	} elsif ($up_check) {
		#print $data[$up_i], "#$data[$i]", "\n";
		print "$data[$i]";
	} elsif ($down_check) {
		#print "#$data[$i]", $data[$down_i], "\n";
		print "$data[$i]";
	} elsif ($chr_gene_num_check) {
		print "$data[$i]";
	}
}
