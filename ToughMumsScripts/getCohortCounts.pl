#!/usr/bin/perl

# Author : Katie Stouffer

# Usage: perl ./getCohortCounts.pl fileWithFileNames
# Note: fileWithFileNames should have full filepath of each one of files to be considered in the cohort; one per line
# Assume total allele count is 2 x number of individuals
# Looks at all snps now (assuming possible errors in variant calling that might not label common snps as common)

use strict;
use warnings;

open my $handle, '<', $ARGV[0];
chomp(my @vcfAnnotTabFilePaths = <$handle>);
close $handle;

my %muts; # Should be $mut{$chrom}{$pos}{$allele} = count
my $posCol = 0;
my $genotypeCol = 5;
my $changeCol = 1;
my $numAlleles = @vcfAnnotTabFilePaths;
my %refAlleles;
$numAlleles = $numAlleles + $numAlleles;

for my $vcfAnnotTabPath (@vcfAnnotTabFilePaths) {
	open (VCFAnnot, $vcfAnnotTabPath);
	my %indMuts;
        print STDERR "Looking into ".$vcfAnnotTabPath."\n";
	while (my $line = <VCFAnnot>) {
		if (index($line, "Position") != -1) {
			next;
		}
		my @lineParts = split("\t", $line);
		$lineParts[$posCol] =~tr/chr//d;
		my @posParts = split(":", $lineParts[$posCol]);
		my $chrom = $posParts[0];
		my $loc = $posParts[1];
		my @changeParts = split(">", $lineParts[$changeCol]);
		my $obsAllele = $changeParts[1];
		my $refAllele = $changeParts[0];
		if (not exists $indMuts{$chrom}{$loc}{$obsAllele}{$lineParts[$genotypeCol]}) {
			$indMuts{$chrom}{$loc}{$obsAllele}{$lineParts[$genotypeCol]} = 0;
			$refAlleles{$chrom}{$loc}=$refAllele
		}
	}
	close(VCFAnnot);
	foreach my $chrom (sort keys %indMuts) {
		foreach my $loc (sort keys %{ $indMuts{$chrom} }) {
			foreach my $obsAllele (sort keys %{ $indMuts{$chrom}{$loc} }) {
				if (exists $muts{$chrom}{$loc}{$obsAllele}) {
					if (exists $indMuts{$chrom}{$loc}{$obsAllele}{"HET"}) {
						$muts{$chrom}{$loc}{$obsAllele} += 1;
					}
					elsif (exists $indMuts{$chrom}{$loc}{$obsAllele}{"HOMO"}) {
						$muts{$chrom}{$loc}{$obsAllele} += 2;
					}
				}
				else {
					if (exists $indMuts{$chrom}{$loc}{$obsAllele}{"HET"}) {
						$muts{$chrom}{$loc}{$obsAllele} = 1;
					}
					elsif (exists $indMuts{$chrom}{$loc}{$obsAllele}{"HOMO"}) {
						$muts{$chrom}{$loc}{$obsAllele} = 2;
					}
				}
			}
		}
	}
}

my @header = qw( Chromosome Position Observed_Allele Allele_Count Assumed_Total_Alleles Reference_Allele );

print join("\t", @header) . "\n";

foreach my $chrom (sort keys %muts) {
	foreach my $loc (sort keys %{ $muts{$chrom} }) {
		foreach my $obsAllele (sort keys %{ $muts{$chrom}{$loc} }) {
			my @line = (
			$chrom, $loc,
			$obsAllele, $muts{$chrom}{$loc}{$obsAllele},
			$numAlleles, $refAlleles{$chrom}{$loc} );
			print join("\t", @line) . "\n";
		}
	}
}

