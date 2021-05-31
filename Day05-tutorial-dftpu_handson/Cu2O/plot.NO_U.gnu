set term post eps enhanced color 24
set output "Cu2O.pdos.gga.eps"
set xlabel "E (eV)"
set ylabel "DOS (states/eV/atom)"
set xrange [-2:15]
set yrange [0:8.5]
set style line 1 lt 2 lc 3 lw 3.0
set style line 2 lt 1 lc 1 lw 3.0
set style line 3 lt 1 lc 2 lw 3.0
plot "./results_Cu2O_GGA/cu2o.pdos_atm#1(Cu)_wfc#1(d)" u 1:2 w l ls 1 t 'Cu-d',\
     "./results_Cu2O_GGA/cu2o.pdos_atm#1(Cu)_wfc#2(s)" u 1:2 w l ls 2 t 'Cu-s',\
      "./results_Cu2O_GGA/cu2o.pdos_atm#5(O)_wfc#2(p)" u 1:2 w l ls 3 t 'O-p'
