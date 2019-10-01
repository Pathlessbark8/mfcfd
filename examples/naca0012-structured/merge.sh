#/bin/bash
cp -r solution temp
cat temp/sol-* > temp.dat
sort -k1 -g temp.dat > sol.dat
rm -rf temp temp.dat
