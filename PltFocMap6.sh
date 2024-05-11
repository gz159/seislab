#!/bin/bash

#python GenEQ.py
##################################################
# A script to plot himalayan range topography
#
##################################################

gmt begin plot_focal_mechanism ps,jpg
    gmt set FONT_ANNOT_PRIMARY 8 FONT_LABEL 10
    gmt set MAP_FRAME_TYPE plain MAP_FRAME_PEN thicker MAP_DEGREE_SYMBOL degree FORMAT_GEO_MAP=ddd:mm:ssF
    gmt set MAP_GRID_CROSS_SIZE_PRIMARY 0


    #p=-Jm102.8/24.2/2.4i
    tr=-R117/123/20/26
    R=-R119/123/21.4/25.8
    J=-Jm120/23/1.5i
    B1=-Bf0.5a1
    B2=-BsWeN
    #p=-JOa76/40/10/4i

    #gmt grdraster 36 $tr -I0.5m -GNinger30s.grd
    gmt grdcut ${HOME}/Bin/China_relief_30s.grd $R -GNinger30s.grd
    # using Aster GDEM topography file
    # grdcut Sichuanchangnin.nc $r -GNinger30so.grd
    gmt grdsample Ninger30s.grd -I0.4k -GNinger30s.grd
    gmt grdgradient Ninger30s.grd -Nt1 -A45 -GNinger30s_i.grd

    ########################## Draw aftershocks relocation ######################### 
    #gmt grd2cpt Ninger30s.grd -CDEM_screen.cpt -Z 

    gmt basemap $R $J $B1 $B2
    gmt makecpt -Cterra -T-10000/10000/200 -Z -D
    gmt grdimage Ninger30s.grd -INinger30s_i.grd
    gmt coast $B1 $B2 -I1 -Df -I1/0.25p,blue -C -W0.25p -A250 
    gmt plot -W0.6p,black,-  chinafault.xy


gmt makecpt -Cseis -T0/40/2

gmt meca -Sa1.5 -A+s0.2c -C << EOF
# lon    lat    depth  strike   dip    rake    mag     lon    lat   title 
121.74  23.81   31.5   228.2    28.5   121.4   7.43   122.31   22.7   202404022358_Mw7.43
121.65  24.1    14.0   264.6    41.2   153.8   6.39   121.45   25.3   202404030011_Mw6.39
121.51  23.76   10.0   217.6    20.2    97.6   5.66   120.11   23.46  202404221050_Mw5.66
121.49  23.85    8.0   230.1    32.9    89.2   5.80   120.0    24.65  202404221411_Mw5.80
121.68  23.7     9.0   205.6    22.6    85.9   6.17   121.58   22.2   202404221826_Mw6.17
121.53  23.81   10.0   216.1    32.7    77.4   6.05   120.03   24.03  202404221832_Mw6.05
121.58  23.69   17.0   224.5    22.1   104.7   5.67   120.18   22.89  202404222049_Mw5.67
121.59  23.87   13.0   216.8    16.1   109.1   6.0    120.7    25.1   202404230004_Mw5.82
121.77  24.08   33.0   208.2    85.4  -118.9   6.0    122.49   24.22  202404261821_Mw5.65
121.77  24.21   27.0   263.8    28.0   162.6   5.28   122.1    25.41  202404261849_Mw5.28
121.62  23.71   20.0   252.0    40.8   128.7   5.46   120.32   22.34  202405060945_Mw5.46
121.89  24.22   11.0    99.5    73.7   13.1    5.77   122.47   24.81  202405100745_Mw5.77
EOF

gmt colorbar -DjBL+w5c/0.5c+ml+o0.8c/0.4c -Bx+l"Focal Depth" -By+lkm -F+gwhite

    # plot legend

    #gmt basemap -R0/2/0/1.5 -Jx7c/1c -B -Y-2.5c 

    #gmt meca -Sa0.7 -C0.5pP2p -W1p,red -Gred << EOF 
    # x  y   depth strike dip rake mag   x    y     title
    #0.2  0.6  9.5   0      90  0   6.0  0.2  0.6   Strike-slip
    #0.4  0.6  9.5   0      45  90  6.0  0.4  0.6   Thrust
    #0.6  0.6  9.5   0      45 -90  6.0  0.6  0.6   Normal
    #0.8  0.6  9.5  30      60  46  6.0  0.8  0.6   Oblique
    #EOF
    #gmt plot -Sa0.5c -Gred -W1.0p,red << EOF
    #1.0 0.6
    #EOF
    #gmt plot -M -W1.0p,black,-  << EOF
    #1.2 0.6
    #1.5 0.6
    #EOF

gmt end show

rm *.grd
