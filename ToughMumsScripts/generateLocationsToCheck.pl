#!/usr/bin/perl

# Usage: perl ./generateLocationsToCheck.pl singleInputAnnotTabFile FileWithLocations

use strict;
use warnings;

my %locationsInFile;
my $chr = '';

my $bedOutput = 0;
if( scalar(@ARGV) > 2 && $ARGV[2] eq "bedfile") {
	$bedOutput = 1;
}

open (FILE, $ARGV[0]);
while (my $line = <FILE>) {
	if (index($line, "Position") != -1) {
		next;
	}
	chomp($line);
	my @lineParts = split("\t", $line);
	my @posParts = split(":", $lineParts[0]);
	if (index($posParts[0], "chr") != -1) {
		$posParts[0] =~ tr/chr//d;
		$chr = 'True';
	}
	$locationsInFile{$posParts[0]}{$posParts[1]} = 1;
}
close(FILE);

open (FILE, $ARGV[1]);
while (my $line = <FILE>) {
	if (index($line, "Chromosome") != -1) {
		next;
	}
	chomp($line);
	my @lineParts = split("\t", $line);
	if (not exists $locationsInFile{$lineParts[0]}{$lineParts[1]}) {
		if ($chr) {
			if($bedOutput) {
				print "chr$lineParts[0]\t".($lineParts[1]-1)."\t".($lineParts[1]-1)."\n";
			}
			else {
				print "chr$lineParts[0]:$lineParts[1]\n";
			}
		}
		else {
			print "$lineParts[0]:$lineParts[1]\n";
		}
	}
}
close(FILE); 
