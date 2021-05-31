/home/sharmila/codes/espresso-4.3.2/bin/pw.x < asy.scf.in > asy.scf.out
cp -r tmp-asy tmpx
cp -r tmp-asy tmpy
cp -r tmp-asy tmpz
/home/sharmila/codes/espresso-4.3.2/bin/pw.x < asy.bscfx.in > asy.bscfx.out
/home/sharmila/codes/espresso-4.3.2/bin/pw.x < asy.bscfy.in > asy.bscfy.out
/home/sharmila/codes/espresso-4.3.2/bin/pw.x < asy.bscfz.in > asy.bscfz.out

#/usr/local/apps/espresso-5.1/bin/pw.x < asy.scf.in > asy.scf.out
#cp -r tmp-asy tmpx
#cp -r tmp-asy tmpy
#cp -r tmp-asy tmpz
#/usr/local/apps/espresso-5.1/bin/pw.x < asy.bscfx.in > asy.bscfx.out
#/usr/local/apps/espresso-5.1/bin/pw.x < asy.bscfy.in > asy.bscfy.out
#/usr/local/apps/espresso-5.1/bin/pw.x < asy.bscfz.in > asy.bscfz.out
