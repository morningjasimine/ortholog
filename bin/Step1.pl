#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use FindBin qw($Bin $Script);
die "Usage: <m8 file> <pep file> <config> <outdir>\n" unless @ARGV == 4;
my $m8_file = shift;
my $pep_file = shift;
my $config = shift;
my $outdir = shift;
mkdir $outdir unless -e $outdir;

foreach my $p (\$m8_file, \$pep_file, \$config, \$outdir) {
	$$p = abs_path($$p);
}
my $select_m8 = "$Bin/../util/select_m8.pl";
my $solar = "$Bin/../util/solar.pl";
my $solar_add_realLen = "$Bin/../util/solar_add_realLen.pl";
my $solar_add_identity = "$Bin/../util/solar_add_identity.pl"; 
foreach my $script ($select_m8, $solar, $solar_add_realLen, $solar_add_identity) {
	die "$script not exist!" unless -e $script;
}
my $step1_out = "$outdir/1_blastp";
mkdir $step1_out unless -e $step1_out;
chdir $step1_out;
############################################################################################################

my @cmd;
push @cmd, "date";
push @cmd, "perl $select_m8 $m8_file $config $step1_out";
push @cmd, "for i in *.m8;do $solar -a prot2prot -f m8 \$i >\$i.solar;done";
push @cmd, "for i in *.solar;do perl $solar_add_realLen \$i $pep_file >\$i.cor;done";
push @cmd, "for i in *.m8;do perl $solar_add_identity --solar \$i.solar.cor --m8 \$i >\$i.solar.cor.idAdd;done";
push @cmd, "rm *.solar *.cor";
push @cmd, "date";

my $sh_file = "$step1_out/step1.run.sh";
open SH, ">$sh_file";
print SH "$_\n" foreach @cmd;
close SH;

system("sh $sh_file");
