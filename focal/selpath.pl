# At first we read in the clip path file
$clippath="selpath.txt";
open(IN,$clippath) || die("Can't open clip path file:$clippath\n");
while($line=<IN>){
   chomp(@line);
   @data=split(', ',$line);
   #   $data[0] longitude    $data[1]  latitude
   push(@clippathpolygon,$data[0],$data[1]);
}
close(IN);

 $pf = "gmtsol";
 chomp($pf);
 open(LINE,$pf) || die("can't open globalCMT file: $pf");
 @lines=<LINE>;

 $fline=qx{gawk 'END{print NR}' $pf};


 $deviashearv="selglobalcmt";
 chomp($deviashearv);
 unlink($deviashearv)if(-e $deviashearv);    
 open(FNR,">$deviashearv") || die("can't creat file:$deviashearv\n");

 for ($i=0;$i<$fline;$i+=1){
    chomp($i);
    $line=$lines[$i];
    @data=split(' ',$line);
    $lon=$data[0];
    $lat=$data[1];
    if(point_in_polygon($lon,$lat,@clippathpolygon)){
        print FNR "$data[0] $data[1] $data[2] $data[3] $data[4] $data[5] $data[6] $data[7] $data[8] $data[9]\n";
    }else{
      print "$data[0] $data[1] $data[2] $data[3] $data[4] $data[5] $data[6] $data[7] $data[8] $data[9]\n";
    }

 }
close(FNR);

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
#print "( 3,-2): ", point_in_polygon( 3,-2, @polygon ), "\n";
#print "( 5, 4): ", point_in_polygon( 5, 4, @polygon ), "\n";
#print "( 5, 1): ", point_in_polygon( 5, 1, @polygon ), "\n";
#print "( 5,-2): ", point_in_polygon( 5,-2, @polygon ), "\n";