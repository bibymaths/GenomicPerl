#!/usr/bin/perl -I . -I ../perllib
################################################################
### Copyright (c) 2002 Rex A. Dwyer.
### Distributed with the book
### Genomic Perl: From Bioinformatics Basics to Working Code
### copyright (c) 2002 Cambridge University Press
### under terms described therein.
################################################################

#########################################
## This program implements a branch & bound method to find the
## optimal multiple alignment of a number of sequences.
#########################################

use strict;
use Util;

my $g=-2; ## gap penalty

#########################################
sub planMultipleAlignment {
##  fills in portions of two dynamic programming tables, %M and %how, 
##  which then tell the score of the best alignment through "relevant"
##  entries of the full table.  Relevance is determined by invoking
##  the subroutine &$relevant, and may be heuristic or exact.
##  RETURNS: best score, and reference to @how array.
#########################################
    my ($relevant,   ## reference to a boolean subroutine.
	@ss)         ## list of strings to be aligned.
	= @_;   
    my %M;     ## holds best alignment scores.
    my %how;  ## holds best predecessor of entries of %M.

    my $numDirections = (1<<@ss)-1;
    my @pendingLayers = (); foreach (@ss) { push @pendingLayers, []; }
    my $origin = join(',' , ((0) x (scalar @ss)));
    my @goal = map { length($_) } @ss;
    my $goalLayer=0; foreach (@goal) { $goalLayer += $_; }
    my $currentLayer = [$origin];
    
    $M{$origin} = 0;
    my $count = 0;
    for (my $layer=0; $layer<=$goalLayer; $layer++) {
	print "LAYER $layer:";
	while (my $vec = pop @$currentLayer) {  ### something in current layer.
	    my @vec = split(',', $vec);
	    my $score = $M{$vec};
	    
	    my $inRange = 1;
	    foreach (0..$#vec) {
		($inRange=0, last) if $vec[$_] > $goal[$_];
	    }
	    if ($inRange && &$relevant($score, \@vec)) {
#		print "Relevant: $vec $M{$vec}\n";
		foreach my $direction (1..$numDirections) {
		    my @succ = @vec;
		    my @column = (("-") x @succ);
		    my $layerDelta = 0;
		    for (my $i=0; $direction; $i++, $direction>>=1) {
			if ($direction & 1) {
			    $succ[$i]++;
			    $column[$i] = substr($ss[$i],$succ[$i]-1,1);
			    $layerDelta++;
			}
		    }
		    my $succ = join(',',@succ);
		    my $nuscore = $score + scoreColumn(@column);
		    if (!defined($M{$succ})) {
			($M{$succ},$how{$succ}) = ($nuscore, $vec);
#			print "push $succ onto list $layerDelta\n";
			push @{$pendingLayers[$layerDelta-1]}, $succ;
		    } elsif ($nuscore >= $M{$succ}) {
			($M{$succ},$how{$succ}) = ($nuscore, $vec);
		    }
		}
	    } else { 
#		print "Not relevant: $vec $M{$vec} $how{$vec}\n";
		delete $M{$vec};
	        delete $how{$vec};
	    }
	}
	push @pendingLayers, [];
	$currentLayer = shift @pendingLayers;
    }
    return ($M{join(',',@goal)}, \%how);
}

#########################################
sub scoreColumn {
##  given a list of symbols in a column,
##  RETURNS: the sum-of-pairs score for the column.
#########################################
    my @col = @_;  ## the symbols in the column.
    my ($gaps,$aas,$score) = (0,0,0);
    while (@col) {
	my $aa = shift @col;
	($gaps++, next) if $aa eq "-";
	$aas++;
	foreach my $aa1 (@col) {
	    next if $aa1 eq "-";
#	    $score += $blosum62{$aa,$aa1};
	    $score += ($aa eq $aa1) ? +1 : -1;
	}
    }
    return $score + ($g * $gaps * $aas);
}

#########################################
sub prependColumnToAlignment {
##  adds a new column to the left end of an existing alignment.
##  RETURNS: nothing.
#########################################
    my ($A,   ## alignment (reference to list of strings)
	@col) ## list of symbols to be added to each string of alignment
	= @_;
    foreach (@$A) {$_ = (shift @col).$_};
}


#########################################
sub reconstructAlignment {
##  constructs a multiple alignment based on individual decisions recorded
##  in dynamic programming table.
##  RETURNS: reference to list of sequences with gap symbols added.
#########################################
    my ($how,  ## table recording how each interim score was derived.
	@ss)   ## list of sequences.
	= @_;
    my @result = (("") x @ss);
    my $origin = join(',', map {"0"} @ss);
    my $current = join(',', map {length($_)} @ss);
    my @gaps = (("-")x@ss);
    while ($current ne $origin) {
	my @current = split(',',$current);
	my $previous = $$how{$current};
	my @previous = split(',',$previous);
	my @column = @gaps;
	for (my $i=0; $i<@ss; $i++) {
	    if ($current[$i] != $previous[$i]) {
		$column[$i] = substr($ss[$i],$previous[$i],1);
	    }
	}
	prependColumnToAlignment(\@result, @column);
	$current = $previous;
    }
    return \@result;
}


#########################################
sub tunnelAlign {
##  implements a simple branch-and-bound heuristic. 
##  Fills only table entries lying near a direct route from source to goal.
##  RETURNS: best score and best alignment found
#########################################
    my ($width,    ## width of tunnel; score non-decreasing function of width
	@ss) = @_; ## list of sequences (strings)
    my @goal = map { length($_) } @ss;
    my $goalSq = 0;
    foreach (@goal) { $goalSq += $_*$_; }

    sub tunnelRelevant {
	my ($score,$vecref) = @_;
	my ($i,$veclen2, $proj, $sslen2) = (0,0,0,0);
	foreach $i (0..@ss-1) {
	    $veclen2 += $$vecref[$i]*$$vecref[$i];
	    $proj += $$vecref[$i]*$goal[$i];
	}
	return (($veclen2 - ($proj*$proj / $goalSq)) < $width*$width);
    }

    my ($score, $how) = planMultipleAlignment(\&tunnelRelevant, @ss);
    return ($score,reconstructAlignment($how,@ss));
}

#########################################
sub computeSuffixPairSimilarities {
##  uses dynamic programming to fill a table of optimal alignment scores
##  for all possible suffixes of all possible pairs of input sequences.
##  RETURNS: a reference to the table of scores.
#########################################
    my(@ss) = @_;   ## list of sequences (strings)
    my @c;
    foreach my $p (0..$#ss) {
	my $s = $ss[$p];
	my $m = length $s;
	foreach my $q ($p+1..$#ss) {
	    my $t = $ss[$q];
	    my $n = length $t;
	    $c[$p][$q][$m][$n] = 0;
	    for (my $i=$m; $i>=0; $i--) { $c[$p][$q][$i][$n] = $g * ($m-$i); }
	    for (my $j=$n; $j>=0; $j--) { $c[$p][$q][$m][$j] = $g * ($n-$j); }
	    for (my $i=$m-1; $i>=0; $i--) {
		for (my $j=$n-1; $j>=0; $j--) {
		    my $match = scoreColumn(substr($s,$i,1),substr($t,$j,1));
		    $c[$p][$q][$i][$j] = max($c[$p][$q][$i+1][$j] + $g,
					     $c[$p][$q][$i][$j+1] + $g,
					     $c[$p][$q][$i+1][$j+1] + $match);
		}
	    }
	}
    }
    return \@c;
}

#########################################
sub branchBoundAlign {
##  implements an exact branch-and-bound method.
##  Fills only table entries that cannot be proved irrelevant by
##  considering pairwise alignment scores.
##  RETURNS: best score and best alignment found.
#########################################
    my ($lowerBound,   ## score of best alignment know so far.
	@ss) = @_;     ## list of sequences to align (strings).

    my $C = computeSuffixPairSimilarities(@ss);

    sub branchBoundRelevant {
	my ($score,$vecref) = @_;
	my $K = @ss;

	my $upperBound = $score;
	foreach my $p (0..$#ss) {
	    foreach my $q ($p+1..$#ss) {
		$upperBound += $$C[$p][$q][$$vecref[$p]][$$vecref[$q]];
	    }
	}
#	print join(',',@$vecref), " return($lowerBound <= $upperBound)\n";
	return ($lowerBound <= $upperBound);
    }

    my ($score, $how) = planMultipleAlignment(\&branchBoundRelevant, @ss);
    return ($score,reconstructAlignment($how,@ss));
}

#########################################
sub scoreMultipleAlignment {
##  sums up the sum-of-pairs scores of all columns of an alignment.
##  RETURNS: the sum-of-pairs score of the entire alignment.
#########################################
    my ($alignment) = @_;  ## reference to list of strings
    my $score;
    foreach my $i (0..length($$alignment[0])-1) {
	$score += scoreColumn(map {substr($_,$i,1)} @$alignment);
    }
    return $score;
}

#########################################
sub printAlignment {
##  prints a multiple alignment with a header.
##  RETURNS: nothing.
#########################################
    my ($title,      ## heading for alignment
	$alignment)  ## reference to list of strings.
	= @_;
    print "\n***** $title\n";
    foreach (@$alignment) { print "$_\n"; }
    my $score = scoreMultipleAlignment($alignment);
    print "score $score\n\n";
}    

#########################################
##  MAIN PROGRAM
#########################################
my @ss;
while (my $s = <DATA>) {
    chomp($s);
    push @ss, $s if $s;  ## tolerate blank lines
}

my ($tunnelScore, $tunnelAlignment) = tunnelAlign(2,@ss);
printAlignment("Tunnel Strategy", $tunnelAlignment);

my ($lcScore, $lcAlignment) = branchBoundAlign($tunnelScore, @ss);
printAlignment("Branch & Bound Strategy", $lcAlignment);



__END__
GVLTDVQVALVKSSFEEFNANIPKNTHRFFTLVLEIAPGAKDLFSFLKGSS
SPLTADEASLVQSSWKAVSHNEVEILAAVFAAYPDIQNKFSQFAGK
VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS
VHLSGGEKSAVTNLWGKVNINELGGEALGRLLVVYPWTQRFFEAFGDLS
VLSAADKTNVKGVFSKIGGHAEEYGAETLERMFIAYPQTKTYFPHFDLS
