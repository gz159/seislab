#!/use/bin/env perl


open(IN, "< ../tmpout.txt");
@catlog= <IN>;
chomp(@catlog);
close(IN);

open(IN, "< ../station.dat");
@sta= <IN>;
chomp(@sta);
close(IN);

open(OUT,"> TravDist.txt") or die("Can not build file outfile");

$i=0;
while($i<@catlog){
	@testchars=split(' ',$catlog[$i]);
	$ftc=$testchars[0];
	chomp($ftc);
	if($ftc eq "#"){  
             $evtlat=$testchars[7];
             chomp($evtlat);
             $evtlon=$testchars[8];
			 chomp($evtlon);
             #print "$evtlat $evtlon\n";
        }else{
             $evstan=$testchars[0];
			 chomp($evstan);
             $evsttrv=$testchars[1];
			 chomp($evsttrv);
             $evstphn=$testchars[3];
			 chomp($evstphn);
             $j=0;
             while($j<@sta){
                  @stainfo=split(' ',$sta[$j]);
                  $stnam=$stainfo[0];
				  chomp($stnam);
           
                  if($evstan eq $stnam){
                      $stlat=$stainfo[1];
					  chomp($stlat);
                      $stlon=$stainfo[2];
					  chomp($stlon);
                      $dist=qx(udelaz -ELAT $evtlat -ELON $evtlon -SLAT $stlat -SLON $stlon -DELKM);
                      chomp($dist);
                      print OUT "$dist $evsttrv $evstphn\n";
                  }
                  $j++;
             }
        }
        $i++;
}
close(OUT);
