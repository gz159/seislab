#! /usr/bin/perl -w
#*********************************************************************************
# A program used to select the gathered disperion data
# the selected disperion data should in the range of c-2*sigma and c+2*sigma
# where "c" is the averaging of phase vlocity and sigma is the standard deviation of 
# phase velocity 
#
#	by guozhi 20101126
#*********************************************************************************


$pwd=qx{pwd};
chomp($pwd);

# Create a directory named "shearV" to keep the inversed shear wave velocity
qx{rm selected -R}if(-e "selected" && -d "selected");
mkdir("selected");


print STDOUT "Begin select the disperion \n";

foreach $period (5..50){
    # First we calculate the averaging phase velocity of each period
    $fr="paths".$period;
    open(LINE,"$fr") || die("can't open input $fr file");
    $filenum=0;
    $sumvel=0;
    while($line = <LINE>){
	chomp($line);
        @data=split(' ',$line);
        $sumvel=$sumvel+$data[4];
        $filenum=$filenum+1;
    }
    close(LINE);
    $averagevel=$sumvel/$filenum;
    # Then we calculate the standard deviation of phase velocity
    open(LINE,"$fr") || die("can't open input $fr file");
    $dvel=0;
    while($line = <LINE>){
        chomp($line);
        @data=split(' ',$line);
        $dvel=$dvel+($data[4]-$averagevel)**2;
    }
    close(LINE);
    $stdvel=sqrt($dvel/$filenum);
    # Get the allowing range of the phase velocity
    $top=$averagevel+2*$stdvel;
    $bottom=$averagevel-2*$stdvel;

    open(LINE,"$fr") || die("can't open input $fr file");
    chdir("selected");
    unlink($fr)if(-e $fr && -d $fr);
    open(FN,">$fr") || die("can't creat file $fr");
    $i=0;
    while($line = <LINE>){
	chomp($line);
        $i++;
        @data=split(' ',$line);
        if($data[4] >= $bottom && $data[4] <= $top){
             print FN "$data[0] $data[1] $data[2] $data[3] $data[4] $data[5]\n";
        } else {
             print STDOUT "Line $i of $fr have been delected\n";
             print STDOUT "v=$data[4],agev=$averagevel,stdv=$stdvel,V range $bottom - $top\n";
        }
    }
    close(FN);
    close(LINE);
    chdir($pwd);
    print STDOUT "****Finish process file: $fr\n";

}

print STDOUT "--------------------------\n";
print STDOUT "Finish all selection\n";


