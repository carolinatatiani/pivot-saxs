set log y
set xlabel "q(k)"
set ylabel "log(I(q))
set style line 1 lt rgb '#d11f60' # Strong Pink
set style line 2 lt rgb '#1f60d1' # Strong Blue
set style line 3 lt rgb '#1fd190' # Strong Cyan
set style line 4 lt rgb '#d1901f' # Strong Orange

set term jpeg size 1024,720 enhanced font 'Helvetica,20' linewidth 3

set output 'theo_compare.jpeg'

plot 'exp.dat' w yerrorbars ls 4 title 'PROT',  'theo.dat' ls 1 w lines title 'SAXS-Factor (Chi=XXXX)', 'pdbmodel00.fit' u 1:4 ls 2 w lines title 'Crysol (Chi=YYYY)', 'pdbmodel_exp.fit' u 1:4 ls 3 w lines title 'FoXS (Chi=ZZZZ)

set xlabel "q(k)*q(k)"
set xrange [:0.016]
set autoscale y
set output 'theo_guinier.jpeg' 

plot 'exp.dat' u ($1*$1):2:3 w yerrorbars ls 4 title 'PROT',  'theo.dat' u ($1*$1):2 ls 1 w lines title 'SAXS-Factor (Chi=XXXX)', 'pdbmodel00.fit' u ($1*$1):4 ls 2 w lines title 'Crysol (Chi=YYYY)', 'pdbmodel_exp.fit' u ($1*$1):4 ls 3 w lines title 'FoXS (Chi=ZZZZ)

quit