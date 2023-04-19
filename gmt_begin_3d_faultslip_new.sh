gmt begin 3d_faultslip png

gmt set MAP_FRAME_TYPE=plain

#gmt basemap -R101/101.5/37.7/38/-20/0  -JX6i/4i -JZ2i -Bxagf -Byagf  -Bzagf  -Bz+l"Depth(km)"  -BWeSnZ -p160/20  
gmt makecpt  -Cjet -T0/3.5/0.2 
gmt plot3d -R101/101.5/37.7/38/-20/0  -JX6i/4i -JZ2i -Bxagf -Byafg  -Bzafg  -Bzafg+l"Depth(km)" -BWeSnZ -p160/20 -W0p,gray lj_slip_3dgrid.gmtline -C -L    
gmt end show