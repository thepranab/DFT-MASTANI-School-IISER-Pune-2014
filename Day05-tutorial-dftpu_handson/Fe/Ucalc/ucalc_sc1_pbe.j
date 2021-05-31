if [ ! -f r.x ]; then
   sh comp_resp_mat.j
fi

echo " &input_mat "> resp_mat0.in
echo "   ntyp = 1 ">> resp_mat0.in
echo "   na(1) = 2 ">> resp_mat0.in
echo "   nalfa = 5 ">> resp_mat0.in
echo "   filepos = 'pos_sc1' ">> resp_mat0.in
echo "   back = 'no' ">> resp_mat0.in
echo "   filednda = 'file_sc1' ">> resp_mat0.in

rm -f Usc1_pbe.out
for n1 in 1 2 3 4 5 6
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
echo $nat $U >> Usc1_pbe.out

done

mv Umat.out Umat_sc1_pbe.out

