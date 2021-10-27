set log y
set xlabel "q(k)"
set ylabel "log(I(q))
set style line 1 lc rgb '#1E90FF' # dark blue
set style line 2 lc rgb '#DC143C' # red
set style line 3 lc rgb '#808080' # grey


plot 'exp.dat' title 'Experimental Data'
replot 'theo.dat' title 'SAXS-factor'
replot 'pdbmodel_exp.fit' u 1:4 title 'FoXS'


set term jpeg size 1024,720 enhanced font 'Helvetica,20' linewidth 2
set output "comp.jpeg"
replot
quit