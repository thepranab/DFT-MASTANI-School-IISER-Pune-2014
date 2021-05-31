#!/bin/bash
#####################

if [ ! -f r.x ]; then
   ./comp_resp_mat.j
fi

#cp ../results_Cu2O/cu*.out ./
#cp ../results_Cu2O/o*.out ./

#rm -f dn*dat
#./grepnalfa_cu2o

echo " &input_mat "> resp_mat0.in
echo "   ntyp = 2 ">> resp_mat0.in 
echo "   na(1) = 16">> resp_mat0.in # Cu-d states
echo "   na(2) = 8 ">> resp_mat0.in # O-p states
echo "   nalfa = 3 ">> resp_mat0.in
echo "   magn = .false. ">> resp_mat0.in
echo "   filepos = 'pos_cu2o' ">> resp_mat0.in
echo "   back = 'neutral' ">> resp_mat0.in
echo "   filednda = 'file.cu2o.24' ">> resp_mat0.in

rm -f Ur24.out
for n1 in 1 2 3 4  
do

rm -f resp_mat.in
cp resp_mat0.in resp_mat.in

echo "   n1 = $n1 ">> resp_mat.in
echo "   n2 = $n1 ">> resp_mat.in
echo "   n3 = $n1 ">> resp_mat.in
echo " &end ">> resp_mat.in

./r.x < resp_mat.in

nat=`grep 'number of atoms' Umat.out |tail -1|awk '{print $9}'`
U=`grep U1 Umat.out |awk '{print $5}'`
echo $nat $U >> Ur24.out

mv Umat.out Umat_cuo.$n1.$n1.$n1.out

done
