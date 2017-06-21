mm
mv rad3d.x pdfv_s6_3

sed -i -e 's/s6/s5/g' params.f90
sed -i -e 's/_6.dat/_5.dat/g' params.f90

mm
mv rad3d.x pdfv_s5_3

sed -i -e 's/s5/s4/g' params.f90
sed -i -e 's/_5.dat/_4.dat/g' params.f90

mm
mv rad3d.x pdfv_s4_3

sed -i -e 's/s4/s3/g' params.f90
sed -i -e 's/_4.dat/_3.dat/g' params.f90

mm
mv rad3d.x pdfv_s3_3

sed -i -e 's/s3/s2/g' params.f90
sed -i -e 's/_3.dat/_2.dat/g' params.f90

mm
mv rad3d.x pdfv_s2_3

sed -i -e 's/s2/s1/g' params.f90
sed -i -e 's/_2.dat/_1.dat/g' params.f90

mm
mv rad3d.x pdfv_s1_3

sed -i -e 's/s1/s6/g' params.f90
sed -i -e 's/_1.dat/_6.dat/g' params.f90

sed -i -e 's/nt_last=700/nt_last=650/g' params.f90

mm
mv rad3d.x pdfv_s6_2

sed -i -e 's/s6/s5/g' params.f90
sed -i -e 's/_6.dat/_5.dat/g' params.f90

mm
mv rad3d.x pdfv_s5_2

sed -i -e 's/s5/s4/g' params.f90
sed -i -e 's/_5.dat/_4.dat/g' params.f90

mm
mv rad3d.x pdfv_s4_2

sed -i -e 's/s4/s3/g' params.f90
sed -i -e 's/_4.dat/_3.dat/g' params.f90

mm
mv rad3d.x pdfv_s3_2

sed -i -e 's/s3/s2/g' params.f90
sed -i -e 's/_3.dat/_2.dat/g' params.f90

mm
mv rad3d.x pdfv_s2_2

sed -i -e 's/s2/s1/g' params.f90
sed -i -e 's/_2.dat/_1.dat/g' params.f90

mm
mv rad3d.x pdfv_s1_2

sed -i -e 's/s1/s6/g' params.f90
sed -i -e 's/_1.dat/_6.dat/g' params.f90

sed -i -e 's/nt_last=650/nt_last=600/g' params.f90

mm
mv rad3d.x pdfv_s6

sed -i -e 's/s6/s5/g' params.f90
sed -i -e 's/_6.dat/_5.dat/g' params.f90

mm
mv rad3d.x pdfv_s5

sed -i -e 's/s5/s4/g' params.f90
sed -i -e 's/_5.dat/_4.dat/g' params.f90

mm
mv rad3d.x pdfv_s4

sed -i -e 's/s4/s3/g' params.f90
sed -i -e 's/_4.dat/_3.dat/g' params.f90

mm
mv rad3d.x pdfv_s3

sed -i -e 's/s3/s2/g' params.f90
sed -i -e 's/_3.dat/_2.dat/g' params.f90

mm
mv rad3d.x pdfv_s2

sed -i -e 's/s2/s1/g' params.f90
sed -i -e 's/_2.dat/_1.dat/g' params.f90

mm
mv rad3d.x pdfv_s1

sed -i -e 's/s1/s6/g' params.f90
sed -i -e 's/_1.dat/_6.dat/g' params.f90
sed -i -e 's/nt_last=600/nt_last=700/g' params.f90
