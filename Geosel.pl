#!/bin/perl
#######################################################################
#  This perl script use to select EQ event in a given region coordinates file
#
#                       guozhi.china(at)gmail.com
#                                        20220210
#######################################################################



#@TienShanpolygon = (39.47,75.42,42.62,70.62,43.82,81.25,41.79,85.39,41.24,80.11);

@periods=("01","02","03","04","05","06");

foreach $period(@periods){

    chomp($period);
    # At first we read in the selection area coordinate file
    $clippath="Haiyuan_subregions/hy".$period.".txt";
    print "Processing file: $clippath\n";
    open(IN,$clippath) || die("Can't open region coordinates file:$clippath\n");
    #empty array contain area coordiantes information
    @clippathpolygon=();
    while($line=<IN>){
       chomp(@line);
       @data=split(' ',$line);
       $a=$data[0];
       $b=$data[1];
       #   $data[0] longitude    $data[1]  latitude
       print $data[0],$data[1];
       push(@clippathpolygon,$data[0],$data[1]);
    }
    close(IN);  
    
    # create file to keep selected EQ
    $deviashearv="EQsel".$period.".dat";
    chomp($deviashearv);
    unlink($deviashearv)if(-e $deviashearv);    
    open(FNR,">$deviashearv") || die("can't creat file:$deviashearv\n");
    
    # read in file containes all of the EQ
    # loop over every EQ to selcet
    $inputf="tmpout.txt";
    print "Processing input EQ file: $inputf\n";
    open(IN,$inputf) || die("Can't open EQ file:$inputf\n");
    while($line=<IN>){
       chomp(@line);
       @data=split(' ',$line);

       if($data[0] eq "#"){
            # if the line start with #, then this line contain EQ information

            $lon=$data[8];
            $lat=$data[7];
            if(point_in_polygon($lon,$lat,@clippathpolygon)){
                $seleq=1;
                print FNR "$line";
            }else{
                $seleq=0;
            }
        
       }else{
            # this line only contain seismic station and travel time information
            if($seleq){
                print FNR "$line";
            }
         
       }

    }# end while loop
    close(IN); 
    close(FNR);

}

# point_in_polygon ( $x, $y, @xy )
#
#    Point ($x,$y), polygon ($x0, $y0, $x1, $y1, ...) in @xy.
#    Returns 1 for strictly interior points, 0 for strictly exterior
#    points. For the boundary points the situation is more complex and
#    beyond the scope of this book.  The boundary points are
#    exact, however: if a plane is divided into several polygons, any
#    given point belongs to exactly one polygon.
#
#    Derived from the comp.graphics.algorithms FAQ,
#    courtesy of Wm. Randolph Franklin.
#
sub point_in_polygon {
    my ( $x, $y, @xy ) = @_;

    my $n = @xy / 2;                      # Number of points in polygon.
    my @i = map { 2 * $_ } 0 .. (@xy/2);  # The even indices of @xy.
    my @x = map { $xy[ $_ ]     } @i;     # Even indices: x-coordinates.
    my @y = map { $xy[ $_ + 1 ] } @i;     # Odd indices: y-coordinates.

    my ( $i, $j );                        # Indices.

    my $side = 0;                         # 0 = outside, 1 = inside.

    for ( $i = 0, $j = $n - 1 ; $i < $n; $j = $i++ ) {
        if (
            (

             # If the y is between the (y-) borders ...
             ( ( $y[ $i ] <= $y ) && ( $y < $y[ $j ] ) ) ||
             ( ( $y[ $j ] <= $y ) && ( $y < $y[ $i ] ) )
            )
            and
            # ...the (x,y) to infinity line crosses the edge
            # from the ith point to the jth point...
            ($x
             <
             ( $x[ $j ] - $x[ $i ] ) *
             ( $y - $y[ $i ] ) / ( $y[ $j ] - $y[ $i ] ) + $x[ $i ] )) {
          $side = not $side; # Jump the fence.
      }
    }

    return $side ? 1 : 0;
}

#@polygon = ( 1, 1,  3, 5,  6, 2,  9, 6,  10, 0,  4,2,  5, -2);
#print "( 3, 4): ", point_in_polygon( 3, 4, @polygon ), "\n";
#print "( 3, 1): ", point_in_polygon( 3, 1, @polygon ), "\n";