set xlabel "Residue number"
set ylabel "Residue number"

set style line 1 lt rgb '#d11f60' # Strong Pink
set style line 2 lt rgb '#1f60d1' # Strong Blue
set style line 3 lt rgb '#1fd190' # Strong Cyan
set style line 4 lt rgb '#d1901f' # Strong Orange

set term jpeg size 1024,720 enhanced font 'Helvetica,20' linewidth 3

set output 'contact_map.jpeg'

plot 'pdbmodel-md.mq' u 2:4 ls 2 notitle

quit
