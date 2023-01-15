IMPFILE=$1
AKAFILE=$2
OUTFILE=$3

cp ${IMPFILE} ${OUTFILE}

awk -F '\t' '{print $1"\t"$2}' ${AKAFILE} | while read GENE; do 

 OLDNAME=$(echo "${GENE}" | cut -f1)
 NEWNAME=$(echo "${GENE}" | cut -f2)

 echo 'Replacing ' $OLDNAME ' with' $NEWNAME

 sed -i "s/${OLDNAME}/${NEWNAME}/g" ${OUTFILE}
 sed -i '/NA/d' ${OUTFILE}

done
