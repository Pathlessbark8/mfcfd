#/bin/bash
cp -r solution temp
sed -i 1d temp/sol-* 
cat temp/sol-* > temp.dat
sort -k1 -g temp.dat > sol.dat
rm -rf temp temp.dat
