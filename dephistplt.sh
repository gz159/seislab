cat tmpout.txt | gawk '{if($1=="#"){print $10} }' | pshistogram -R0/40/0/20 -Ba10f5:"Depth(km)":/a5f5:"Frequency"::,%:WSne -JX8c/5c -S -Z1 -W1 -C -L0.5p -V -K > tmpp.ps
cat hypoDD.reloc | gawk '{print $4}' | pshistogram -R -B  -JX -S -Z1 -W1 -C -L0.5p,red -V -O >> tmpp.ps
