#!/usr/bin/env -S bash -e
gmt begin CN-border-JM ps,jpg
    gmt set MAP_GRID_PEN_PRIMARY 0.25p,gray,2_2
    # 绘制中国地图
    gmt coast -JM105/35/15c -R70/138/13/56 -Ba10f5g10 -G244/243/239 -S167/194/223
    gmt basemap -Lg85/17.5+c17.5+w800k+f+u --FONT_ANNOT_PRIMARY=4p
    gmt plot CN-border-La.gmt -W0.1p
    gmt plot selpath.txt -W1p,red
    cat CEICselesta.dat | gawk '{print $6,$5}' | gmt plot -Gblue -St0.1c -W0.1p,white

    # 绘制南海区域
    gmt inset begin -DjRB+w1.8c/2.2c -F+p0.5p
        gmt coast -JM? -R105/123/3/24 -G244/243/239 -S167/194/223 -Df
        gmt plot CN-border-La.gmt -W0.1p
        cat CEICselesta.dat | gawk '{print $6,$5}' | gmt plot -Gblue -St0.1c -W0.1p,white
    gmt inset end
gmt end show
