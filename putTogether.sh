#author: Jennifer Havens
#dependancy for Snakefile.linearCount
#adds thrid col (the read counts) to growing raw read count expression matrix
#before running need to have expression matrix with window location names in first col
#first argument is cellName second argument is output matix name

printf "$1 \n" > temp
awk '{print $3}' $1 >> temp
paste $2 temp > temp2
mv temp2 $2