#! /usr/bin/perl
#########################################################################
# A program decompress current directory seed into daily sac file
# and storage seed into coresspond directory
#
# usage:
#       put this program into directory that contain seed file then run
#       it, and here we assume seed file name Do NOT end with "pl".
# input :
#               current directory seed file
# output :
#           correspond directory that contain seed and sac file
#########################################################################


$pwd=qx{pwd};
chomp($pwd);

#-----------------------------------------------------------------------
#   note:
#   Because we need handle seed file that start year 2002
#   day of 91 to 121.so here we define $year 2002,$day 91 and
#   $days 31.
#   When you need this program to handle others date seed file
#   you need modify those value to you need.
#   $statlen is the station name characters length in seed file
# 
#-----------------------------------------------------------------------
$year=2016;
$beginday=0;
#$day=91;
$days=365;
#$statlen=4;


chomp($year);
chomp($day);
chomp($days);

#-----------------------------------------------------------------------
# First we create directory that contain correspond seed file
#-----------------------------------------------------------------------
@filenames=qx{ls};
foreach $fn(@filenames){
   chomp($fn);
   if(-e $fn && substr($fn,-2,2) ne "pl"){
	#The seed file's name contain information about
	#network name and station name
	@sf=split(/\./,$fn);
	$statname=$sf[1];
	chomp($statname);
	mkdir($statname);
	#Copy seed file to correspond directory
	rename($fn,$statname.'/'.$fn);
    }
}
@statime=("0:00","0:30","1:00","1:30","2:00","2:30","3:00","3:30","4:00","4:30","5:00","5:30","6:00","6:30"
          ,"7:00","7:30","8:00","8:30","9:00","9:30","10:00","10:30","11:00","11:30","12:00","12:30","13:00"
          ,"13:30","14:00","14:30","15:00","15:30","16:00","16:30","17:00","17:30","18:00","18:30"
          ,"19:00","19:30","20:00","20:30","21:00","21:30","22:00","22:30","23:00");
@endtime=("1:00","1:30","2:00","2:30","3:00","3:30","4:00","4:30","5:00","5:30","6:00","6:30","7:00","7:30"
          ,"8:00","8:30","9:00","9:30","10:00","10:30","11:00","11:30","12:00","12:30","13:00","13:30","14:00"
          ,"14:30","15:00","15:30","16:00","16:30","17:00","17:30","18:00","18:30","19:00","19:30"
          ,"20:00","20:30","21:00","21:30","22:00","22:30","23:00","23:30","23:59");
#------------------------------------------------------------------------
# then we enter every directory that contain seed file and decompress
# each seed in daily lenght sac data.
#------------------------------------------------------------------------
@dirnames=qx{ls};
foreach $fn(@dirnames){
   chomp($fn);
   print "%%%%%%%%%%%%%%%%%%%$fn\n";
   if(-d $fn){
	chdir($fn);
	$temp=qx{pwd};
	print "+++++++++++++++++$temp";
	$seedfile=qx{ls};
	chomp($seedfile);
#	print "*************************$seedfile\n";
	$day=$beginday;
	open(RDSEED, "|rdseed\n");
	# cycling decompress seed file into daily lenght sac file
	for($i=0;$i<$days;$i++){
            for($j=0;$j<47;$j++){
               print "###########################$seedfile";
               #Input  File (/dev/nrst0) or 'Quit' to Exit: YLDINX20020401-20020501.4844
               print RDSEED "$seedfile\n";
               #Output File (stdout)    :
                
               # print "---------------------$seedfile \n"; 
               print RDSEED "\n";
               #Volume #  [(1)-N]       :
               print RDSEED "\n";
               #Options [acCsSpRtde]    : d
               print RDSEED "d\n";
               #Summary file (None)     :
               print RDSEED "\n";
               #Station List (ALL)      :
               print RDSEED "\n";
               #Channel List (ALL)      :
               print RDSEED "\n";
               #Network List (ALL)      :
               print RDSEED "\n";
               #Loc Ids (ALL ["--" for spaces]) :
               print RDSEED "\n";
               #Output Format [(1=SAC), 2=AH, 3=CSS, 4=mini seed, 5=seed, 6=sac ascii, 7=SEGY] :
               print RDSEED "1\n";
               #Output file names include endtime? [Y/(N)](New in 5.0 version rdseed)
               print RDSEED "\n";
               #Output poles & zeroes ? [Y/(N)]
               print RDSEED "\n";
               #Check Reversal [(0=No), 1=Dip.Azimuth, 2=Gain, 3=Both]:
               print RDSEED "\n";
               #Select Data Type [(E=Everything), D=Data of Undetermined State, R=Raw waveform Data, Q=QC'd data] :
               print RDSEED "\n";
               #Start Time(s) YYYY,DDD,HH:MM:SS.FFFF :
               print RDSEED "$year,$day,$statime[$j]\n";
               print "$year,$day,$statime[$j]\n";
               #print RDSEED $day,;
               #print RDSEED $statime[$j];
               #print RDSEED "\n";
               #End Time(s)   YYYY,DDD,HH:MM:SS.FFFF :
               print RDSEED "$year,$day,$endtime[$j]\n";
               #print RDSEED $day,;
               #print RDSEED $endtime[$j];
               #print RDSEED "\n";
               #Sample Buffer Length [2000000]:
               print RDSEED "\n";
               #Extract Responses [Y/(N)]     :
               if($i==0){
                   print RDSEED "y\n";
               }else{
                   print RDSEED "\n";
               }
            }
            $day++;
	}
	#Input  File (/dev/nrst0) or 'Quit' to Exit: quit
	print RDSEED "quit\n";
	close(RDSEED);
	# For the cautiou of failure of genetating resopnse file
        qx{rdseed -f $seedfile -R};
	chdir($pwd);	
    }
}

#-----------------------------------------------------------------
# next we enter in the generated directory and delete rdseed log file
# and rename seismeter response file to RESP
#-----------------------------------------------------------------


$pwd = qx{pwd};
chomp($pwd);
@dirlist = qx{ls};
foreach $dirname(@dirlist){
    chomp($dirname);
    if(-e $dirname && -d $dirname){
      chdir($dirname);
      @filelist=qx{ls};
      foreach $fl(@filelist){
         chomp($fl);
         if(substr($fl,-3,3) eq "SAC"){
            print "$fl :is a SAC file ,Jump it\n";
            next;
         }elsif(substr($fl,0,4) eq "RESP" && substr($fl,-3,2) eq "BH"){
            $respfn=substr($fl,0,4);
            print " Rename Seismeter response file $fl to $respfn\n";
            qx{mv $fl $respfn};
         }else{
            print "Delete other file\n";
            qx{rm $fl}
         }
         
      }
      chdir($pwd);  
    }
}
