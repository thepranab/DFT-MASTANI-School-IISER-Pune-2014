#!/bin/bash
#PBS -N lfpo
#PBS -o lfpo.o
#PBS -e lfpo.e
#PBS -l select=2:ncpus=16
#PBS -l walltime=6:00:00
##PBS -q P_theos

cd $PBS_O_WORKDIR

module load intelmpi/4.1.3
module load intel/14.0.1
module load fftw/3.3.3-intel_intelmpi-13.0.1_4.1.0

PWdir='/scratch/cococcio/Bin'
tempdir1='/scratch/cococcio/FePO4/temp1gga'
tempdir2='/scratch/cococcio/FePO4/temp2gga'

if [ ! -d $tempdir1 ]; then
   mkdir  $tempdir1
fi

if [ ! -d $tempdir2 ]; then
   mkdir  $tempdir2
fi

resdir='/scratch/cococcio/FePO4/results_fepo4_Ucalc_gga'

if [ ! -d $resdir ]; then
   mkdir  $resdir
fi

pseudodir='/scratch/cococcio/pseudi'

U=1.d-10

ecut=65.0
ecutrho=`echo "$ecut * 8" |bc -l`

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=0.6,
starting_magnetization(3)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
!startingwfc='file'
!startingpot='file'
!diago_thr_init = 2.0E-11
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
Fe3 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.00000000     0.00000000     0.00000000
Fe2      0.49999864     0.00000000    -0.41393362
Fe3     -0.04946574     0.50000000    -0.49995407
Fe3      0.45052136     0.50000000    -0.91403681
P       -0.17960241     0.00000000    -0.55113785
P        0.63011874     0.50000000    -0.36301502
P        0.13012077     0.50000000    -0.05093371
P        0.32039682     0.00000000    -0.86279486
O       -0.15607465     0.00000000    -0.24285946
O        0.60672077     0.50000000    -0.67134041
O        0.10676155     0.50000000    -0.74258636
O        0.34392697     0.00000000    -0.17106974
O        0.17079443     0.00000000    -0.79852584
O        0.27968840     0.50000000    -0.11554598
O       -0.22030405     0.50000000    -0.29851658
O        0.67079501     0.00000000    -0.61541718
O       -0.10551538    -0.20278462    -0.69922038
O        0.55588345     0.70261140    -0.21495383
O        0.05586523     0.70261110    -0.19895061
O        0.39448098    -0.20279149    -0.71474814
O        0.55588345     0.29738860    -0.21495383
O       -0.10551538     0.20278462    -0.69922038
O        0.39448098     0.20279149    -0.71474814
O        0.05586523     0.29738890    -0.19895061
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.Fe1.unp.out

ethr=`grep ethr $resdir/fepo4.Fe1.unp.out |tail -1 |awk '{print $3}'`

rm -rf $tempdir2/*
cp -r $tempdir1/* $tempdir2/

for alp in 0.0 -0.04 0.04 -0.08 0.08 #0.0 -0.02 0.02 -0.05 0.05
do

rm -rf $tempdir1/*
cp -r $tempdir2/* $tempdir1/

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=0.6,
starting_magnetization(3)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
hub_pot_fix = .true.,
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
Hubbard_alpha(1) = $alp
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
startingwfc='file'
startingpot='file'
diago_thr_init = $ethr
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
Fe3 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.00000000     0.00000000     0.00000000
Fe2      0.49999864     0.00000000    -0.41393362
Fe3     -0.04946574     0.50000000    -0.49995407
Fe3      0.45052136     0.50000000    -0.91403681
P       -0.17960241     0.00000000    -0.55113785
P        0.63011874     0.50000000    -0.36301502
P        0.13012077     0.50000000    -0.05093371
P        0.32039682     0.00000000    -0.86279486
O       -0.15607465     0.00000000    -0.24285946
O        0.60672077     0.50000000    -0.67134041
O        0.10676155     0.50000000    -0.74258636
O        0.34392697     0.00000000    -0.17106974
O        0.17079443     0.00000000    -0.79852584
O        0.27968840     0.50000000    -0.11554598
O       -0.22030405     0.50000000    -0.29851658
O        0.67079501     0.00000000    -0.61541718
O       -0.10551538    -0.20278462    -0.69922038
O        0.55588345     0.70261140    -0.21495383
O        0.05586523     0.70261110    -0.19895061
O        0.39448098    -0.20279149    -0.71474814
O        0.55588345     0.29738860    -0.21495383
O       -0.10551538     0.20278462    -0.69922038
O        0.39448098     0.20279149    -0.71474814
O        0.05586523     0.29738890    -0.19895061
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.Fe1.$alp.out

done 

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
!startingwfc='file'
!startingpot='file'
!diago_thr_init = 2.0E-11
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P1  31.0 P.pbe-n-rrkjus_psl.0.1.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.17960241     0.00000000     0.55113785
Fe1      0.67960106     0.00000000     0.13720423
Fe2      0.13013667     0.50000000     0.05118378
Fe2      0.63012377     0.50000000    -0.36289895
P1       0.00000000     0.00000000     0.00000000
P        0.80972115     0.50000000     0.18812283
P        0.30972318     0.50000000     0.50020414
P        0.49999923     0.00000000    -0.31165701
O        0.02352776     0.00000000     0.30827839
O        0.78632318     0.50000000    -0.12020255
O        0.28636397     0.50000000    -0.19144851
O        0.52352938     0.00000000     0.38006811
O        0.35039685     0.00000000    -0.24738798
O        0.45929082     0.50000000     0.43559187
O       -0.04070163     0.50000000     0.25262128
O        0.85039742     0.00000000    -0.06427932
O        0.07408703    -0.20278462    -0.14808253
O        0.73548586     0.70261140     0.33618402
O        0.23546764     0.70261110     0.35218725
O        0.57408340    -0.20279149    -0.16361029
O        0.73548586     0.29738860     0.33618402
O        0.07408703     0.20278462    -0.14808253
O        0.57408340     0.20279149    -0.16361029
O        0.23546764     0.29738890     0.35218725
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.P1.unp.out

ethr=`grep ethr $resdir/fepo4.P1.unp.out |tail -1 |awk '{print $3}'`

rm -rf $tempdir2/*
cp -r $tempdir1/* $tempdir2/

for alp in 0.0 -0.04 0.04 -0.08 0.08 #0.0 -0.02 0.02 -0.05 0.05
do

rm -rf $tempdir1/*
cp -r $tempdir2/* $tempdir1/

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
hub_pot_fix = .true.,
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
Hubbard_alpha(3) = $alp
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
startingwfc='file'
startingpot='file'
diago_thr_init = $ethr
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P1  31.0 P.pbe-n-rrkjus_psl.0.1.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.17960241     0.00000000     0.55113785
Fe1      0.67960106     0.00000000     0.13720423
Fe2      0.13013667     0.50000000     0.05118378
Fe2      0.63012377     0.50000000    -0.36289895
P1       0.00000000     0.00000000     0.00000000
P        0.80972115     0.50000000     0.18812283
P        0.30972318     0.50000000     0.50020414
P        0.49999923     0.00000000    -0.31165701
O        0.02352776     0.00000000     0.30827839
O        0.78632318     0.50000000    -0.12020255
O        0.28636397     0.50000000    -0.19144851
O        0.52352938     0.00000000     0.38006811
O        0.35039685     0.00000000    -0.24738798
O        0.45929082     0.50000000     0.43559187
O       -0.04070163     0.50000000     0.25262128
O        0.85039742     0.00000000    -0.06427932
O        0.07408703    -0.20278462    -0.14808253
O        0.73548586     0.70261140     0.33618402
O        0.23546764     0.70261110     0.35218725
O        0.57408340    -0.20279149    -0.16361029
O        0.73548586     0.29738860     0.33618402
O        0.07408703     0.20278462    -0.14808253
O        0.57408340     0.20279149    -0.16361029
O        0.23546764     0.29738890     0.35218725
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.P1.$alp.out

done 

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
!startingwfc='file'
!startingpot='file'
!diago_thr_init = 2.0E-11
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.15607465     0.00000000     0.24285946
Fe1      0.65607330     0.00000000    -0.17107416
Fe2      0.10660891     0.50000000    -0.25709461
Fe2      0.60659601     0.50000000    -0.67117735
P       -0.02352776     0.00000000    -0.30827839
P        0.78619339     0.50000000    -0.12015556
P        0.28619542     0.50000000     0.19192575
P        0.47647147     0.00000000    -0.61993540
O1       0.00000000     0.00000000     0.00000000
O        0.76279542     0.50000000    -0.42848094
O        0.26283621     0.50000000    -0.49972690
O        0.50000162     0.00000000     0.07178972
O        0.32686909     0.00000000    -0.55566638
O        0.43576306     0.50000000     0.12731348
O       -0.06422939     0.50000000    -0.05565712
O        0.82686966     0.00000000    -0.37255772
O        0.05055927    -0.20278462    -0.45636092
O        0.71195810     0.70261140     0.02790563
O        0.21193988     0.70261110     0.04390885
O        0.55055564    -0.20279149    -0.47188868
O        0.71195810     0.29738860     0.02790563
O        0.05055927     0.20278462    -0.45636092
O        0.55055564     0.20279149    -0.47188868
O        0.21193988     0.29738890     0.04390885
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O1.unp.out

ethr=`grep ethr $resdir/fepo4.O1.unp.out |tail -1 |awk '{print $3}'`

rm -rf $tempdir2/*
cp -r $tempdir1/* $tempdir2/

for alp in 0.0 -0.04 0.04 -0.08 0.08 #0.0 -0.02 0.02 -0.05 0.05
do

rm -rf $tempdir1/*
cp -r $tempdir2/* $tempdir1/

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
hub_pot_fix = .true.,
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
Hubbard_alpha(4) = $alp
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
startingwfc='file'
startingpot='file'
diago_thr_init = $ethr
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.15607465     0.00000000     0.24285946
Fe1      0.65607330     0.00000000    -0.17107416
Fe2      0.10660891     0.50000000    -0.25709461
Fe2      0.60659601     0.50000000    -0.67117735
P       -0.02352776     0.00000000    -0.30827839
P        0.78619339     0.50000000    -0.12015556
P        0.28619542     0.50000000     0.19192575
P        0.47647147     0.00000000    -0.61993540
O1       0.00000000     0.00000000     0.00000000
O        0.76279542     0.50000000    -0.42848094
O        0.26283621     0.50000000    -0.49972690
O        0.50000162     0.00000000     0.07178972
O        0.32686909     0.00000000    -0.55566638
O        0.43576306     0.50000000     0.12731348
O       -0.06422939     0.50000000    -0.05565712
O        0.82686966     0.00000000    -0.37255772
O        0.05055927    -0.20278462    -0.45636092
O        0.71195810     0.70261140     0.02790563
O        0.21193988     0.70261110     0.04390885
O        0.55055564    -0.20279149    -0.47188868
O        0.71195810     0.29738860     0.02790563
O        0.05055927     0.20278462    -0.45636092
O        0.55055564     0.20279149    -0.47188868
O        0.21193988     0.29738890     0.04390885
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O1.$alp.out

done 

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
!startingwfc='file'
!startingpot='file'
!diago_thr_init = 2.0E-11
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1     -0.17079443     0.00000000     0.79852584
Fe1      0.32920421     0.00000000     0.38459221
Fe2     -0.22026018     0.50000000     0.29857176
Fe2      0.27972693     0.50000000    -0.11551097
P       -0.35039685     0.00000000     0.24738798
P        0.45932430     0.50000000     0.43551082
P       -0.04067367     0.50000000     0.74759213
P        0.14960238     0.00000000    -0.06426902
O       -0.32686909     0.00000000     0.55566638
O        0.43592634     0.50000000     0.12718543
O       -0.06403288     0.50000000     0.05593948
O        0.17313254     0.00000000     0.62745609
O1       0.00000000     0.00000000     0.00000000
O        0.10889397     0.50000000     0.68297986
O       -0.39109848     0.50000000     0.50000926
O        0.50000058     0.00000000     0.18310866
O       -0.27630981    -0.20278462     0.09930545
O        0.38508902     0.70261140     0.58357201
O       -0.11492921     0.70261110     0.59957523
O        0.22368655    -0.20279149     0.08377769
O        0.38508902     0.29738860     0.58357201
O       -0.27630981     0.20278462     0.09930545
O        0.22368655     0.20279149     0.08377769
O       -0.11492921     0.29738890     0.59957523
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O2.unp.out

ethr=`grep ethr $resdir/fepo4.O2.unp.out |tail -1 |awk '{print $3}'`

rm -rf $tempdir2/*
cp -r $tempdir1/* $tempdir2/

for alp in 0.0 -0.04 0.04 -0.08 0.08 #0.0 -0.02 0.02 -0.05 0.05
do

rm -rf $tempdir1/*
cp -r $tempdir2/* $tempdir1/

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
hub_pot_fix = .true.,
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
Hubbard_alpha(4) = $alp
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
startingwfc='file'
startingpot='file'
diago_thr_init = $ethr
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1     -0.17079443     0.00000000     0.79852584
Fe1      0.32920421     0.00000000     0.38459221
Fe2     -0.22026018     0.50000000     0.29857176
Fe2      0.27972693     0.50000000    -0.11551097
P       -0.35039685     0.00000000     0.24738798
P        0.45932430     0.50000000     0.43551082
P       -0.04067367     0.50000000     0.74759213
P        0.14960238     0.00000000    -0.06426902
O       -0.32686909     0.00000000     0.55566638
O        0.43592634     0.50000000     0.12718543
O       -0.06403288     0.50000000     0.05593948
O        0.17313254     0.00000000     0.62745609
O1       0.00000000     0.00000000     0.00000000
O        0.10889397     0.50000000     0.68297986
O       -0.39109848     0.50000000     0.50000926
O        0.50000058     0.00000000     0.18310866
O       -0.27630981    -0.20278462     0.09930545
O        0.38508902     0.70261140     0.58357201
O       -0.11492921     0.70261110     0.59957523
O        0.22368655    -0.20279149     0.08377769
O        0.38508902     0.29738860     0.58357201
O       -0.27630981     0.20278462     0.09930545
O        0.22368655     0.20279149     0.08377769
O       -0.11492921     0.29738890     0.59957523
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O2.$alp.out

done 

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
!startingwfc='file'
!startingpot='file'
!diago_thr_init = 2.0E-11
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.10551538     0.20278462     0.69922038
Fe1      0.60551402     0.20278462     0.28528676
Fe2      0.05604964     0.70278462     0.19926631
Fe2      0.55603674     0.70278462    -0.21481642
P       -0.07408703     0.20278462     0.14808253
P        0.73563412     0.70278462     0.33620536
P        0.23563615     0.70278462     0.64828667
P        0.42591220     0.20278462    -0.16357448
O       -0.05055927     0.20278462     0.45636092
O        0.71223615     0.70278462     0.02787998
O        0.21227693     0.70278462    -0.04336598
O        0.44944235     0.20278462     0.52815064
O        0.27630981     0.20278462    -0.09930545
O        0.38520378     0.70278462     0.58367440
O       -0.11478866     0.70278462     0.40070381
O        0.77631039     0.20278462     0.08380320
O1       0.00000000     0.00000000     0.00000000
O        0.66139883     0.90539602     0.48426655
O        0.16138061     0.90539572     0.50026978
O        0.49999636    -0.00000687    -0.01552776
O        0.66139883     0.50017322     0.48426655
O        0.00000000     0.40556925     0.00000000
O        0.49999636     0.40557612    -0.01552776
O        0.16138061     0.50017352     0.50026978
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O3.unp.out

ethr=`grep ethr $resdir/fepo4.O3.unp.out |tail -1 |awk '{print $3}'`

rm -rf $tempdir2/*
cp -r $tempdir1/* $tempdir2/

for alp in 0.0 -0.04 0.04 -0.08 0.08 #0.0 -0.02 0.02 -0.05 0.05
do

rm -rf $tempdir1/*
cp -r $tempdir2/* $tempdir1/

rm -rf fepo4.in
cat > fepo4.in << EOF
&control
calculation='scf'
restart_mode='from_scratch',
prefix='FePO4',
pseudo_dir = '$pseudodir',
outdir='$tempdir1',
tprnfor = .false.
tstress = .false.
verbosity = 'high'
etot_conv_thr = 1.d-5,
forc_conv_thr = 1.d-4,
/
&system
ibrav=0,
celldm(1) = 18.868676564
nat=24, ntyp=5,
nspin = 2
starting_magnetization(1)=0.6,
starting_magnetization(2)=-0.6,
ecutwfc = $ecut,
ecutrho = $ecutrho,
occupations='smearing',
smearing='gauss',
degauss=0.01
hub_pot_fix = .true.,
lda_plus_u = .true.
lda_plus_u_kind = 0,
U_projection_type = 'ortho-atomic',
Hubbard_U(1) = $U
Hubbard_U(2) = $U
Hubbard_U(3) = $U
Hubbard_U(4) = $U
Hubbard_U(5) = $U
Hubbard_U(6) = $U
Hubbard_alpha(4) = $alp
/
&electrons
conv_thr = 1.0e-8
mixing_beta = 0.25
mixing_ndim = 20
mixing_mode = 'local-TF'
startingwfc='file'
startingpot='file'
diago_thr_init = $ethr
electron_maxstep = 200
/
&ions
upscale = 2.d0,
trust_radius_ini = 0.4
/
&cell
press = 0.0
!cell_dofree = 'z'
/
ATOMIC_SPECIES
Fe1 56.0 Fe.pbe-nd-rrkjus.UPF
Fe2 56.0 Fe.pbe-nd-rrkjus.UPF
P   31.0 P.pbe-n-rrkjus_psl.0.1.UPF
O1  16.0 O.pbe-rrkjus.UPF
O   16.0 O.pbe-rrkjus.UPF
ATOMIC_POSITIONS {crystal}
Fe1      0.10551538     0.20278462     0.69922038
Fe1      0.60551402     0.20278462     0.28528676
Fe2      0.05604964     0.70278462     0.19926631
Fe2      0.55603674     0.70278462    -0.21481642
P       -0.07408703     0.20278462     0.14808253
P        0.73563412     0.70278462     0.33620536
P        0.23563615     0.70278462     0.64828667
P        0.42591220     0.20278462    -0.16357448
O       -0.05055927     0.20278462     0.45636092
O        0.71223615     0.70278462     0.02787998
O        0.21227693     0.70278462    -0.04336598
O        0.44944235     0.20278462     0.52815064
O        0.27630981     0.20278462    -0.09930545
O        0.38520378     0.70278462     0.58367440
O       -0.11478866     0.70278462     0.40070381
O        0.77631039     0.20278462     0.08380320
O1       0.00000000     0.00000000     0.00000000
O        0.66139883     0.90539602     0.48426655
O        0.16138061     0.90539572     0.50026978
O        0.49999636    -0.00000687    -0.01552776
O        0.66139883     0.50017322     0.48426655
O        0.00000000     0.40556925     0.00000000
O        0.49999636     0.40557612    -0.01552776
O        0.16138061     0.50017352     0.50026978
K_POINTS {automatic}
2 4 4 0 0 0
CELL_PARAMETERS (alat)
   1.000000000   0.000000000   0.000000000
   0.000000000   0.592611671   0.000000000
   0.000000000   0.000000000   0.489214880
EOF


mpirun -genv I_MPI_FABRICS shm:tmi -np 32 $PWdir/pw.x -npool 2  < fepo4.in > $resdir/fepo4.O3.$alp.out

done 
