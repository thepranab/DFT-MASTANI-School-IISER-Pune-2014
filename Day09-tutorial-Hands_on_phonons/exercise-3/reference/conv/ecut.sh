#!/bin/sh
NAME='ecut'
RUNDIR='/home/sharmila/Documents/MASTANI-2014/phonons/exercise-3'
for a in 45 50 55 60    
#for a in 10 15 20 25 30 35 40    
  do
  cat > $RUNDIR/$NAME.$a.in <<EOF
&control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='alas',
    tstress = .true.
    tprnfor = .true.
    pseudo_dir = '../../pseudo',
    outdir='tmp'
 /
 &system
    ibrav=  2,  nat=  2, ntyp= 2,celldm(1)=10.50,
    ecutwfc =$a,
    occupations='fixed',
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.5
    conv_thr =  1.0d-10
 /
 &ions
 /
 &cell
 /

ATOMIC_SPECIES
 Al  26.98  Al.pz-vbc.UPF
 As  74.92  As.pz-bhs.UPF 

ATOMIC_POSITIONS
 Al 0.00 0.00 0.00
 As 0.25 0.25 0.25

K_POINTS AUTOMATIC
4 4 4 1 1 1
EOF
pw.x < $RUNDIR/$NAME.$a.in > $RUNDIR/$NAME.$a.out
done

