#!/bin/bash

if [ ! -f r.x ]; then
   ./comp_resp_mat.j
fi

cp ../results_NiO/ni*.out ./
cp ../results_NiO/o*.out ./

rm dn*dat
./grepnalfa_nio_r16

echo " &input_mat "> resp_mat0.in
echo "   ntyp = 2 ">> resp_mat0.in #
echo "   na(1) = 8 ">> resp_mat0.in
echo "   na(2) = 8 ">> resp_mat0.in
echo "   nalfa = 3 ">> resp_mat0.in
echo "   magn = .true. ">> resp_mat0.in
echo "   filepos = 'pos_nio_r16' ">> resp_mat0.in
echo "   back = 'no' ">> resp_mat0.in
echo "   filednda = 'file.nio.r16' ">> resp_mat0.in

rm Ur16.out
for n1 in 1 2 3 4 5
do

rm resp_mat.in
cp resp_mat0.in resp_mat.in

echo "   n1 = $n1 ">> resp_mat.in
echo "   n2 = $n1 ">> resp_mat.in
echo "   n3 = $n1 ">> resp_mat.in
echo " &end ">> resp_mat.in

./r.x < resp_mat.in

nat=`grep 'number of atoms' Umat.out |tail -1|awk '{print $9}'`
U=`grep U1 Umat.out |awk '{print $5}'`
echo $nat $U >> Ur16.out

mv Umat.out Umat_nio.$n1.$n1.$n1.out

done
