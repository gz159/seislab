gmt begin CN-focalmechanism jpg
    gmt set MAP_GRID_PEN_PRIMARY 0.25p,gray,2_2
    gmt coast -JM105/35/10c -R70/138/20/54 -Ba10f5g10 -G244/243/239 -S167/194/223
    gmt basemap -Lg85/17.5+c17.5+w800k+f+u --FONT_ANNOT_PRIMARY=1p
    gmt grdimage @earth_relief_06m -I+d
    gmt plot CN-border-La.gmt -W0.1p
    gmt meca selglobalcmt -Gred -Sz0.3 -V
gmt end show