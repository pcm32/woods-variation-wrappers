#!/usr/bin/perl

# Author : Katie Stouffer

# Usage: perl ./caseControlCheck.pl --filterOut=outputFromFilterRunWithAffected --family=numberOfFamilyInFilterOutput [optional args]

# Dominant Check: file with filenames of family members who are not affected; check to ensure that mutation in affected is not found in non-affected family members
# Use parameters: --dom --unAffected=fileNameWithUnaffectedFiles (one per line)

# Recessive Check
# Consanguineous: checks to ensure that mother and father each only have 1 copy (i.e. HET) for the mutation in affected individuals; unaffected can only have 1 copy if any
# Use parameters: --rec --con --mother=filenameOfMother --father=filenameOfFather (if only one file is specified, then only that one is checked for having 1 copy) --unAffected=fileNameWithUnaffected

# Non-consanguineous: checks to ensure that mother and father each only have 1 copy of 1 of the mutations in the affected individuals (i.e. no parent has the same two het mutations or homo for one of the mutations).  If non-affected sibs are sequenced, checks to ensure that no non-affected sib has both HET mutations or homo for one of the mutations.
# Use parameters: --rec --mother=filenameOfMother --father=filenameOfFather --unAffected=fileNameWithUnaffected (or sibs unaffected)

use strict;
use warnings;
use Getopt::Long;

my $filterOut = '';
my $mother = '';
my $father = '';
my $unAffected = '';
my $rec = '';
my $dom = '';
my $con = '';
my $family = 0;

usage() if (@ARGV < 3 or ! GetOptions('filterOut=s' => \$filterOut,
	    'unAffected:s' => \$unAffected,
	    'mother:s' => \$mother,
	    'father:s' => \$father,
	    'rec' => \$rec,
	    'dom' => \$dom,
	    'con' => \$con,
	    'family=i' => \$family));
sub usage
{
	print "usage: ./caseControlCheck.pl --family=numberOfFamilyInFilterOutput --filterOut=outputFromFilterWithAffected [--dom] [--rec] [--con] [--unAffected=fileNamesOfUnaffected(sibs or non-relations)] [--mother=fileNameOfMother] [--father=fileNameOfFather]\n";
	exit;
}

my @unAffectedFileNames;
if ($unAffected) {
	open my $handle, '<', $unAffected;
	chomp(@unAffectedFileNames = <$handle>);
	close $handle;
}

# Store per gene all of the mutations in that gene
my %affectedLocs;
my %familyCount;
open (FILE, $filterOut);
while (my $line = <FILE>) {
	# assume first line begins with genename
	chomp($line);
	if ($line=~/^\s*$/) {
		next;
	}
	if (index($line, "-----") != -1) {
		next;
	}
	$line =~ tr/://d;
	my $gene = $line;
	$line = <FILE>; # first dashed line
	$line = <FILE>; # second dashed line
	# Get family number
	$line = <FILE>;
	chomp($line);
	$line =~ tr/://d;
	$line =~ tr/Family//d;
	$line =~ s/^\s+|\s+$//g;
	my $exitLoop = 0;
	while ($line != $family) {
		while (index($line, "-----") == -1) {
			$line = <FILE>;
		}
		$line = <FILE>; # blank line
		$line = <FILE>; # next family or another dash
		if (index($line, "Family") != -1) {
			chomp($line);
			$line =~ tr/://d;
			$line =~ tr/Family//d;
			$line =~ s/^\s+|\s+$//g;
		}
		else {
			$exitLoop = 1;
			last;
		}
	}
	if ($exitLoop) {
		$line = <FILE>; # dashed line
		next;
	}
	# Found family with genes
	else {
		$line = <FILE>; # blank line
		$line = <FILE>; # location
		while (index($line, "*") != -1) {
			chomp($line);
			$line =~ tr/*//d;
			$line =~ s/^\s+|\s+$//g;
			my @locParts = split(":",$line); # [0] is chromosome, [1] is base pair number
			$line = <FILE>; # file name of first individual
			$line = <FILE>; # mutation in Individual (use to get base change)
		 	my @lineParts = split("\t", $line); # [0] is location chrom:number, [1] is mutation refAllele>obsAlle
			#             Gene   Chrom         Base pair     Mutation
			$affectedLocs{$gene}{$locParts[0]}{$locParts[1]}{$lineParts[1]} = 1;
			while (index($line, "*") == -1) {
				$line = <FILE>;
				if (index($line, "-----") != -1) {
					last;
				}
			}
		}
		$line = <FILE>; # blank line
		# Skip down to 2 dashed lines
		while ((index($line, ":") != -1) || (($line=~/^\s*$/) || (index($line, "/") != -1))) {
			$line = <FILE>;
			if (index($line, "-----") != -1) {
				$line = <FILE>;
			}
		}
	}
}
close(FILE);

my %locsToEliminate; #chr, pos, baseChange
my $changeCol = 1;
my $geneCol = 6;
my $Genotype = 5;

# Dominant Case
if ($dom) {
	my $numLocsElim = 0;
	if ($mother) {
		push @unAffectedFileNames, $mother;
	}
	if ($father) {
		push @unAffectedFileNames, $father;
	}
	for my $file (@unAffectedFileNames) {
		open (FILE, $file);
		while (my $line = <FILE>) {
			if (index($line, "Position") != -1) {
				next;
			}
			chomp($line);
			my @lineParts = split("\t", $line);
			my @posParts = split(":", $lineParts[0]);
			$posParts[0] =~ tr/chr//d; 
			if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
				$locsToEliminate{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
			}
		}
		close(FILE);
	}
	foreach my $g (sort keys %locsToEliminate) {
		foreach my $c (sort keys %{ $locsToEliminate{$g}}) {
			foreach my $p (sort keys %{ $locsToEliminate{$g}{$c}}) {
				foreach my $chge (sort keys %{ $locsToEliminate{$g}{$c}{$p}}) {
					$numLocsElim++;
				}
			}
		}
	}
}
elsif ($rec && $con) {
	if ($unAffected) {
		for my $file (@unAffectedFileNames) {
			open (FILE, $file);
			while (my $line = <FILE>) {
				if (index($line, "Position") != -1) {
					next;
				}
				chomp($line);
				my @lineParts = split("\t", $line);
				my @posParts = split(":", $lineParts[0]);
				$posParts[0] =~ tr/chr//d;
				if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
					if ($lineParts[$Genotype] eq "HOMO") {
						$locsToEliminate{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
					}
				}
			}
			close(FILE);
		}
	}
	my %affectedLocsHET;
	my @parents;
	if ($mother) {
		push @parents, $mother;
	}
	if ($father) {
		push @parents, $father;
	}
	if ($mother || $father) {
		for my $p (@parents) {
			open (FILE, $p);
			while (my $line = <FILE>) {
				if (index($line, "Position") != -1) {
					next;
				}
				chomp($line);
				my @lineParts = split("\t", $line);
				my @posParts = split(":", $lineParts[0]);
				$posParts[0] =~ tr/chr//d;
				if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
					if ($lineParts[$Genotype] eq "HET") {
						$affectedLocsHET{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}{$p} = 1;
					}
				}
			}
			close(FILE);
		}
	}
	my $numParents = @parents;
	# add locs to eliminate all those that are not HET
	foreach my $g (sort keys %affectedLocs) {
		foreach my $c (sort keys %{ $affectedLocs{$g}}) {
			foreach my $p (sort keys %{ $affectedLocs{$g}{$c}}) {
				foreach my $chge (sort keys %{ $affectedLocs{$g}{$c}{$p}}) {
					foreach my $par (@parents) {
						if (not exists $affectedLocsHET{$g}{$c}{$p}{$chge}{$par}) {
							$locsToEliminate{$g}{$c}{$p}{$chge} = 1;
						}	
					}
				}
			}
		}
	}
}

# Assume recessive case is left
else {
	# Make Dictionary of total number of mutations per gene
	my %numberOfMuts;
	foreach my $g (sort keys %affectedLocs) {
		my $count = 0;
		foreach my $c (sort keys %{ $affectedLocs{$g}}) {
			foreach my $p (sort keys %{ $affectedLocs{$g}{$c}}) {
				foreach my $chge (sort keys %{ $affectedLocs{$g}{$c}{$p}}) {
					$count++;
				}
			}
		}
		$numberOfMuts{$g} = $count;
	}
	
	# Check for any of mutations being HOMO in unaffected individuals or individuals having all of mutations associated with gene
	if ($unAffected || ($mother || $father)) {
		if ($mother) {
			push @unAffectedFileNames, $mother;
		}
		if ($father) {
			push @unAffectedFileNames, $father;
		}
		for my $file (@unAffectedFileNames) {
			open (FILE, $file);
			my %indLocs;
			while (my $line = <FILE>) {
				if (index($line, "Position") != -1) {
					next;
				}
				chomp($line);
				my @lineParts = split("\t", $line);
				my @posParts = split(":", $lineParts[0]);
				$posParts[0] =~ tr/chr//d;
				if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
					if ($lineParts[$Genotype] eq "HOMO") {
						$locsToEliminate{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
					}
					$indLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
				}
			}
			close(FILE);
			
			# Count up number of mutations per gene in person and if equals number in total dictionary, eliminate all of the mutations in gene
			foreach my $g (sort keys %indLocs) {
				my $count = 0;
				foreach my $c (sort keys %{ $indLocs{$g}}) {
					foreach my $p (sort keys %{ $indLocs{$g}{$c}}) {
						foreach my $chge (sort keys %{ $indLocs{$g}{$c}{$p}}) {
							$count++;
						}
					}
				}
				if ($count == $numberOfMuts{$g}) {
					foreach my $c (sort keys %{ $indLocs{$g}}) {
						foreach my $p (sort keys %{ $indLocs{$g}{$c}}) {
							foreach my $chge (sort keys %{ $indLocs{$g}{$c}{$p}}) {
								$locsToEliminate{$g}{$c}{$p}{$chge} = 1;
							}	
						}
					}
				}
			}
		}
	}
	# If have Mother and Father, eliminate mutations that are both in mother and father or that are in neither mother and father
	if ($mother && $father) {
		my %motherLocs;
		my %fatherLocs;
		# Record locations that are HET per gene
		open (FILE, $mother);
		while (my $line = <FILE>) {
			if (index($line, "Position") != -1) {
				next;
			}
			chomp($line);
			my @lineParts = split("\t", $line);
			my @posParts = split(":", $lineParts[0]);
			$posParts[0] =~ tr/chr//d;
			if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
				if ($lineParts[$Genotype] eq "HET") {
					$motherLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
				}
			}
		}
		close(FILE);
		open (FILE, $father);
		while (my $line = <FILE>) {
			if (index($line, "Position") != -1) {
				next;
			}
			chomp($line);
			my @lineParts = split("\t", $line);
			my @posParts = split(":", $lineParts[0]);
			$posParts[0] =~ tr/chr//d;
			if (exists $affectedLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]}) {
				if ($lineParts[$Genotype] eq "HET") {
					$fatherLocs{$lineParts[$geneCol]}{$posParts[0]}{$posParts[1]}{$lineParts[$changeCol]} = 1;
				}
			}
		}
		close(FILE);
		# for each gene, check each pair of mutations and if pair splits between parents, add both locs to locstokeep
		my %locsToKeep;
		foreach my $g (sort keys %affectedLocs) {
			my %muts;
			foreach my $c (sort keys %{ $affectedLocs{$g}}) {
				foreach my $p (sort keys %{ $affectedLocs{$g}{$c}}) {
					foreach my $chge (sort keys %{ $affectedLocs{$g}{$c}{$p}}) {
						my $s = "$c:$p:$chge";
						$muts{$s} = 1;
					}
				}
			}
			my @mutsArray = keys %muts;
			my $numMuts = @mutsArray;
			for (my $i = 0; $i < $numMuts; $i++) {
				my @mut1Parts = split(":", $mutsArray[$i]);
				my $motherNotFather1 = ((exists $motherLocs{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]}) && (not exists $fatherLocs{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]}));
				my $fatherNotMother1 = ((not exists $motherLocs{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]}) && (exists $fatherLocs{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]}));
				if (not ($motherNotFather1 || $fatherNotMother1)) {
					next;
				}
				for (my $j = $i+1; $j < $numMuts; $j++) {
					my @mut2Parts = split(":", $mutsArray[$j]);
					my $motherNotFather2 = ((exists $motherLocs{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]}) && (not exists $fatherLocs{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]}));
					my $fatherNotMother2 = ((not exists $motherLocs{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]}) && (exists $fatherLocs{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]}));
					if (not ($motherNotFather2 || $fatherNotMother2)) {
						next;
					}
					elsif ($motherNotFather1 && $fatherNotMother2) {
						$locsToKeep{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]} = 1;
						$locsToKeep{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]} = 1;
					}
					elsif ($motherNotFather2 && $fatherNotMother1) {
						$locsToKeep{$g}{$mut1Parts[0]}{$mut1Parts[1]}{$mut1Parts[2]} = 1;
						$locsToKeep{$g}{$mut2Parts[0]}{$mut2Parts[1]}{$mut2Parts[2]} = 1;
					}
				}
			}	
		}
		# Add all affected locs not in locs to keep to locs to remove
		foreach my $g (sort keys %affectedLocs) {
			foreach my $c (sort keys %{ $affectedLocs{$g}}) {
				foreach my $p (sort keys %{ $affectedLocs{$g}{$c}}) {
					foreach my $chge (sort keys %{ $affectedLocs{$g}{$c}{$p}}) {
						if (not exists $locsToKeep{$g}{$c}{$p}{$chge}) {
							$locsToEliminate{$g}{$c}{$p}{$chge} = 1;
						}	
					}
				}
			}
		}
	}
}
# Eliminate all Locs from AffectedLocs based on locsToEliminate
my %newAffectedLocs;
foreach my $g (sort keys %affectedLocs) {
	my $countOfMuts = 0;
	my %tempLocs;
	foreach my $c (sort keys %{ $affectedLocs{$g}}) {
		foreach my $p (sort keys %{ $affectedLocs{$g}{$c}}) {
			foreach my $chge (sort keys %{ $affectedLocs{$g}{$c}{$p}}) {
				if (not exists $locsToEliminate{$g}{$c}{$p}{$chge}) {
					$tempLocs{$c}{$p}{$chge} = 1;
					$countOfMuts++;
				}
			}
		}
	}
	if (((not $rec) || ($rec && $con)) || (($rec && (not $con)) && ($countOfMuts > 1))) {
		foreach my $c (sort keys %tempLocs) {
			foreach my $p (sort keys %{ $tempLocs{$c}}) {
				foreach my $chge (sort keys %{ $tempLocs{$c}{$p}}) {
					$newAffectedLocs{$g}{$c}{$p}{$chge} = 1;
				}
			}
		}
	}
}

# Read in original file and skip down to genes and families that kept, if location is in newAffectedLocs with base change per individual, then write the line

open (FILE, $filterOut);
while (my $line = <FILE>) {
	my $gene;
	# assume first line is Gene
	if ($line=~/^\s*$/) {
		print $line;
		next;
	}
	chomp($line);
	$line =~ tr/://d;
	$gene = $line;
	# if Gene was not altered, then print all of lines to next gene
	if (not exists $affectedLocs{$line}) {
		print "$line:\n";
		$line = <FILE>;
		print $line;
		$line = <FILE>;
		print $line;
		$line = <FILE>;
		while (index($line, "-----") == -1) {
			print $line;
			$line = <FILE>;
			if (index($line, "-----") != -1) {
				print $line;
				$line = <FILE>;
			}
		}
		print $line;
	}
	else {
		my $printFamily = '';
		if (exists $newAffectedLocs{$line}) {
			$printFamily = 'True';
		}
		my @geneHeader;
		push @geneHeader, "$line:\n";
		$line = <FILE>; # dashed line 1
		push @geneHeader, $line;
		$line = <FILE>; # dashed line 2
		push @geneHeader, $line;
		$line = <FILE>; # Family 1
		while (index($line, "Family") == 0) { 
			chomp($line);
			$line =~ tr/://d;
			$line =~ tr/Family//d;
			$line =~ s/^\s+|\s+$//g;
			if ($line == $family) {
				if (not $printFamily) {
					while (index($line, "-----") == -1) {
						$line = <FILE>;
					}
					$line = <FILE>; # blank line
					$line = <FILE>; # next family or dashed line
				}
				else {
					push @geneHeader, "Family $line:\n";
					$line = <FILE>;
					while (index($line, "-----") == -1) {
						if ($line=~/^\s*$/) {
							push @geneHeader, "$line";
							$line = <FILE>;
						}
						elsif (index($line, "*") != -1) {
							chomp($line);
							my @locParts = split(":", $line);
							$locParts[0] =~ tr/*//d;
							$locParts[1] =~ tr/*//d;
							$locParts[0] =~ s/^\s+|\s+$//g;
							$locParts[1] =~ s/^\s+|\s+$//g;
							if (exists $newAffectedLocs{$gene}{$locParts[0]}{$locParts[1]}) {
								# print lines corresponding to this location and base changes
								push @geneHeader, "$line\n";
								$line = <FILE>;
								while ((index($line, "*") == -1) && (index($line, "-----") == -1)) {
									if ($line=~/^\s*$/) {
										push @geneHeader, $line;
									}
									elsif (index($line, "/") != -1) {
										push @geneHeader, $line;
									}
									else {
										my @subLineParts = split("\t", $line);
										if (exists $newAffectedLocs{$gene}{$locParts[0]}{$locParts[1]}{$subLineParts[1]}) {
											push @geneHeader, $line;
										}
									}
									$line = <FILE>;
								}
							}
							else {
								$line = <FILE>;
								while ((index($line, "*") == -1) && (index($line, "-----") == -1)) {
									$line = <FILE>;
								}
							}
						}
					}
					push @geneHeader, "$line";
					$line = <FILE>; # blank line
					push @geneHeader, "$line";
					$line = <FILE>; # next family or dashed line
				}
			}
			else {
				push @geneHeader, "Family $line:\n";
				while (index($line, "-----") == -1) {
					$line = <FILE>;
					push @geneHeader, $line;
				}
				$line = <FILE>; # blank line
				push @geneHeader, $line;
				$line = <FILE>; # next family or dashed line
			}
		}
		# Dashed line ending families for this gene
		if (@geneHeader <= 3) {
			$line = <FILE>; # dashed line
			$line = <FILE>; # blank line
		}
		else {
			push @geneHeader, $line; # dashed line
			$line = <FILE>; # dashed line 2
			push @geneHeader, $line;
			$line = <FILE>; # blank line
			push @geneHeader, $line;
			foreach my $gh (@geneHeader) {
				print $gh;
			}
		}
	}
}
close(FILE);	
