
$dd=tojulian(2015,12,31);

print $dd;



open(MY, '>myfile.txt') or die ("Could not open file\n");
print MY "The date is $dd\n";
close(MY);


#========================================================================
# tojulian-perl routine convert from calendar date to julian date
#
# 22-Sep-1996 T.Hutchinson hutch@hummock.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  A perl subroutine that returns the julian day (# days into year).
#
#    input: date in yyyy,mm,dd format
#
#    return: date in yyyy,ddd format
#

sub tojulian {

    #my($year) = substr($_[0], 0, 4);
    #my($month) = substr($_[0], 4, 2);
    #my($day) = substr($_[0], 6, 2);
    my ( $year, $month, $day ) = @_;
    my($juldate,$month_count);

    my(@daysinmonth) = 
        ("31","28","31","30","31","30","31","31","30","31","30","31");
    my($julday) = 0;

    for ($month_count = 1; $month_count < $month; $month_count++) {
        $julday += $daysinmonth[$month_count-1];
    }
    
    $julday += $day;

    if (($month > 2) && (isLeapyear($year))) { $julday++; }

    $juldate=sprintf("%d,%03d", $year, $julday);
    return $juldate;
}

#========================================================================
# isLeapyear-perl routine to determine if year is a leap year
#
# 22-Sep-1996 T.Hutchinson hutch@hummock.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  Check to see if year is a leap year.
#
#    input: year (4 digits)
#
#    return: True ("1") for leap year, False ("0") for non-leap year.
#

sub isLeapyear {
    
    my($year) = $_[0];
    
    if (($year % 4 == 0) && ((($year % 100) != 0) || (($year % 400) == 0))) {
	return 1;
    } else {
	return 0;
    }
}

#========================================================================
# fromjulian-perl routine to convert from julian date to calendar date
#
# 22-Sep-1996 T.Hutchinson hutch@hummock.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  A perl subroutine that returns the calendar day.
#
#    input: date in yyyyddd format
#
#    return: date in yyyymmdd format
#

sub fromjulian {
    my($julday) = $_[0];
    my($day) = substr($_[0],4,3);
    my($year) = substr($_[0],0,4);
    my($month);
    my($date);

# a trick to remove leading zeros from $day
    $day += 0;

    my(@daysinmonth) = 
        ("31","28","31","30","31","30","31","31","30","31","30","31");

    if ((isLeapyear($year)) && ($day >= 60)) { $daysinmonth[1] = 29; }
    
    for ($month = 0;$day > $daysinmonth[$month]; $month++) {
	$day -= $daysinmonth[$month];
    }

    $month++;
    $date = sprintf("%4d%02d%02d", $year, $month, $day);
    return $date;
}

#========================================================================
# cal_plus-perl routine to add or subtract from the calendar date
#
# 22-Sep-1996 T.Hutchinson hutch@hummock.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  Get the next day of the year.
#
#    input: first: date in yyyymmdd format
#           second: number of days to add (subtract if < 0) to date
#
#    return: next day in yyyymmdd format.
#

sub cal_plus {

    my($yyyyddd) = tojulian($_[0]);
    my($diff) = $_[1];
    $yyyyddd = julian_plus($yyyyddd,$_[1]);
    my $outdate = fromjulian($yyyyddd);
    return $outdate;
}

#========================================================================
# julian_plus-perl routine to add or subtract from the julian date
#
# 22-Sep-1996 T.Hutchinson hutch@hummock.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  Get the next day of the year.
#
#    input: $juldate: in yyyyddd format
#           $diff: number of days to add (subtract if < 0) to $juldate
#
#    return: next day in yyyyddd format.
#
sub julian_plus
{
    
    my($juldate, $diff) = ($_[0], $_[1]);
    my $day = substr($juldate,4,3);
    my $year = substr($juldate,0,4);
    if ($diff > 0) {
	while ($diff > 0) {
	    my $days_this_year = &isLeapyear($year) ? 366 : 365;
	    if ($diff > 364) {
		$juldate = &julian_plus($juldate, 364);
		$day = substr($juldate,4,3);
		$year = substr($juldate,0,4);
		$diff -= 364;
	    } else {
		$day += $diff;
		$diff = 0;
		if ($day > $days_this_year) {
		    $year++;
		    $day -= $days_this_year;
		}
	    }
	}
    } elsif ($diff < 0) {
	while ($diff < 0) {
	    my $days_previous_year = &isLeapyear($year-1) ? 366 : 365;
	    if ($diff < -364) {
		$juldate = &julian_plus($juldate, -364);
		$day = substr($juldate,4,3);
		$year = substr($juldate,0,4);
		$diff += 364;
	    } else {
		$day += $diff;
		$diff = 0;
		if ($day < 1) {
		    $year--;
		    $day += $days_previous_year;
		}
	    }
	}
    }
    $juldate = sprintf("%d%03d", $year, $day);
    return $juldate;
}

#========================================================================
# date2_to_date - converts yymmdd format to yyyymmdd format
#
# 11-May-1998 T.Haran tharan@colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#========================================================================
#
# $Header: /home/haran/navdir/src/scripts/date.pl,v 1.28 2003/08/05 16:55:34 haran Exp $
#
#  Convert a date in yymmdd format to yyyymmdd format
#
#    input: $date2: in yymmdd format
#
#    return: $date in yyyymmdd format.
#
sub date2_to_date
{
    
    my($date2) = @_;

    my $date = (substr($date2, 0, 2) >= 70) ?
	"19" . "$date2" : "20" . "$date2";
    return $date;
}



#==============================================================================
# date_ok - Checks to see if date is in valid yyyymmdd format.
#
# 21-Aug-1997 Terry Haran haran@kryos.colorado.edu 303-492-1847
# National Snow & Ice Data Center, University of Colorado, Boulder
# Boulder, CO  80309-0449
#==============================================================================
#
#
#  Checks to see if date is in valid format.
#
#    input:  date in yyyymmdd format
#            string identifying date
#
#    return: $true if date is ok.
#            $false otherwise.
#

sub date_ok {

    my $date   = $_[0];
    my $string = $_[1];

    my $yyyy = substr($date,0,4);
    my $mm = substr($date,4,2);
    my $dd = substr($date,6,2);
    my $success = $false;
    if (($yyyy < 1970) || ($yyyy >= 2070)) {
	print "ERROR: $string year is out of range, year: $yyyy.\n";
    } elsif (($mm < 1) || ($mm > 12)) {
	print "ERROR: $string month is out of range, month: $mm.\n";
    } elsif (($dd < 1) || ($dd > 31)) {
	print "ERROR: $string day is out of range, day: $dd.\n";
    } else {
	my $juldate = &tojulian($date);
	my $date_check = &fromjulian($juldate);
	if ($date ne $date_check) {
	    print "ERROR: $string is an invalid date: $date.\n";
	} else {
	    $success = $true;
	}
    }
    return $success;
}

