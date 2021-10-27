set log y
set xlabel "q(k)"
set ylabel "log(I(q))
set style line 1 lt rgb '#d11f60' # Strong Pink
set style line 2 lt rgb '#1f60d1' # Strong Blue
set style line 3 lt rgb '#1fd190' # Strong Cyan
set style line 4 lt rgb '#d1901f' # Strong Orange

set term jpeg size 1024,720 enhanced font 'Helvetica,20' linewidth 2

set output "RUUN_all4.jpeg"
plot '../exp.dat' w yerrorbars ls 4 title 'PROT', 'conformationLAST.dat' ls 3 title 'Last conformation (chi=LASCHI)', 'conformationCHIM.dat' ls 1 title 'Smallest Chi (chi=CHIII)','conformationENEER.dat' ls 2 title 'Smallest HS_Energy (chi=ENECHI)'

set output "RUUN_last.jpeg"
plot '../exp.dat' w yerrorbars ls 4 title 'PROT', 'conformationLAST.dat' ls 3 title 'Last conformation (chi=LASCHI)'


set output "RUUN_chi.jpeg"
plot '../exp.dat' w yerrorbars ls 4 title 'PROT', 'conformationCHIM.dat' ls 1 title 'Smallest Chi (chi=CHIII)'


set output "RUUN_energy.jpeg"
plot '../exp.dat' w yerrorbars ls 4 title 'PROT', 'conformationENEER.dat' ls 2 title 'Smallest Energy_{HS} (chi=ENECHI)'

quit