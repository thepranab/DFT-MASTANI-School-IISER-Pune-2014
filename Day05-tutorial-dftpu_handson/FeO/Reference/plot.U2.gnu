set term post eps enhanced color 24
set output "FeO.pdos.U2.eps"
set xlabel "E (eV)"
set ylabel "DOS (states/eV/atom)"
set xrange [-10.0:5.0]
set yrange [0.0:4.0]
set key top left
set style line 1 lt 2 lc 3 lw 3.0
set style line 2 lt 1 lc 1 lw 3.0
set style line 3 lt 1 lc 2 lw 3.0
set style line 4 lt 4 lc 7 lw 2.0
set arrow from 0.0,0.0 to 0.0,4.0 nohead ls 4
plot "./results_FeO_U2/feo.pdos_atm#3(Fe1)_wfc#2(d)" u (-10.1559+$1):2 w l ls 1 t 'Fe-d (maj)',\
     "./results_FeO_U2/feo.pdos_atm#3(Fe1)_wfc#2(d)" u (-10.1559+$1):3 w l ls 2 t 'Fe-d (min)',\
      "./results_FeO_U2/feo.pdos_atm#1(O)_wfc#2(p)" u (-10.1559+$1):($2+$3) w l ls 3 t 'O-p'
