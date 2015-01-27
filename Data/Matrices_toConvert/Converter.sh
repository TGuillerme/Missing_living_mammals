for f in *.nex
do
prefix=$(basename $f .nex)
echo $prefix

sed 's/#NEXUS/#NEXUS\
\
Begin data;\
    Dimensions ntax=@ nchar=§;\
    Format datatype=standard gap=- missing=?;\
    Matrix\
/g' ${prefix}.nex |
sed '$s/$/\
    ;\
End;/' >  ${prefix}.nex2

variables=$(grep [0-9]sp  ${prefix}.nex)
species=$(echo $variables | sed 's/sp.*//')
character=$(echo $variables | sed 's/.*sp//' | sed 's/c//')

mv ${prefix}.nex ${prefix}.bkp

sed 's/'"$variables"'//g'  ${prefix}.nex2 | sed 's/@/'"$species"'/' | sed 's/§/'"$character"'/' >  ${prefix}.nex

rm  ${prefix}.nex2

done