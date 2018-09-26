#!/usr/bin/perl -w
use strict;
die "Usage: <CR.mutipleCopy.pos> <IR.pep>\n" unless @ARGV == 2;

my %pep;
open IN, $ARGV[1];
$/ = ">";
<IN>;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//g;
	/(.+)\n/;
	my $id = (split /\s+/, $1)[0];
	s/.+\n//;
	s/U/X/ig;
	$pep{$id} = $_;
}
$/ = "\n";
close IN;

open IN, $ARGV[0];
while (<IN>) {
	my $id = (split /\s+/)[0];
	print ">$id\n$pep{$id}\n";
}
close IN;
