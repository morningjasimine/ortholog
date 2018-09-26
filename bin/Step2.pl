#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
use File::Basename;
die "Usage: <1_blastp dir> <config> <gff dir> <outdir>\n" unless @ARGV == 5;
my $blastp_dir = shift;
my $contig = shift;
my $gff_dir = shift;
my $outdir = shift;

foreach my $p (\$blastp_dir, \$contig, \$gff_dir, \$outdir) {
	$$p = abs_path($$p);
}
my $ort1 = "$Bin/../util/ort1.pl";
my $ort2 = "$Bin/../util/ort2.pl";
my $link_ort = "$Bin/../util/link_ort.pl";
my $syntenic_ort = "$Bin/../util/syntenic_ort.pl";
foreach my $script ($ort1, $ort2, $link_ort, $syntenic_ort) {
	die "$script not exist!" unless -e $script;
}
my $step2_out = "$outdir/2_ortholog";
mkdir $step2_out unless -e $step2_out;
chdir $step2_out;
############################################################################################################

## get gff file
my %gff;
opendir DH, $gff_dir;
while (my $file = readdir DH) {
	next unless $file =~ /\.gff$/;	
	my $in_file = "$gff_dir/$file";
	my $sp = (split /\./, $file)[0];
	$gff{$sp} = $in_file;
}
closedir DH;

## read contig and get cutoff
my %cutoff; 
open IN, $contig;
$/ = "#--END--";
while (<IN>) {
	print STDERR "test:$_";
	next unless /alignment\s+cutoff/;
	s/^\s+|\s+$//g;
	my @lines = split /\n/;
	foreach my $line (@lines) {
		next if $line =~ /^#/;
		my @info = split /\s+/, $line;
		$cutoff{$info[0]}{$info[1]} = [$info[2], $info[3]];
		$cutoff{$info[1]}{$info[0]} = [$info[2], $info[3]];
	}
	last;
}
$/ = "\n";
close IN;


my @ort_files;
open SH, ">step2.run.sh";
print SH "date\n";
opendir DH, $blastp_dir;
while (my $file = readdir DH) {
	next unless $file =~ /\.idAdd$/;
	my $in_file = "$blastp_dir/$file";
	my $sp_pair = (split /\./, $file)[0];
	my ($ref_sp, $tar_sp) = (split /_/, $sp_pair);
	my ($alignRate, $identity) = @{$cutoff{$ref_sp}{$tar_sp}};
	my $basename = basename($in_file);
	print SH <<cmd;
perl $ort1 $in_file $ref_sp $tar_sp >$basename.rbh;
awk '\$3>=$alignRate && \$4>=$alignRate && \$5>=$identity' $basename.rbh >$basename.rbh.filter
perl $ort2 $basename.rbh.filter $gff{$ref_sp} $gff{$tar_sp} $ref_sp $tar_sp >$basename.rbh.filter.ort
perl $syntenic_ort $basename.rbh.filter.ort $gff_dir >$basename.rbh.filter.ort.syntenic

cmd
	push @ort_files, "$basename.rbh.filter.ort";
}
closedir DH;
my $n = @ort_files;
$n ++;
print SH "perl $link_ort @ort_files $contig >${n}species.ort\n";
print SH "date\n";
close SH;

system("sh step2.run.sh");
