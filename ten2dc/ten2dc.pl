#!/usr/bin/perl -w

use strict;
@ARGV == 6 or die "
Usage: perl $0 Mrr Mtt Mff Mrt Mrf Mtf
	(the exponent is not taken into account)\n
";

my ($mrr, $mtt, $mff, $mrt, $mrf, $mtf) = @ARGV;

my $line1 = `./ten2axe $mrr $mtt $mff $mrt $mrf $mtf`;
my @axis = split /\s+/,$line1;
print STDERR "scalar_moment= $axis[14] iso= $axis[12] clvd= $axis[13]\n";

my $line2=`./axe2dc $axis[10] $axis[11] $axis[4] $axis[5]`;

print STDERR "strike1 dip1 rake1 strike2 dip2 rake2:\n";
print STDERR "\t$line2";
