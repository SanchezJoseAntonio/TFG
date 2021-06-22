#Jerarquia: tRNA > sncRNA > repeats > pseudogene > lncRNA > protein coding
#Circanno devuelve: chr coordinate coordinate peak:information:type reads strand ...
#chr1    2235865 2235925 G10.1:ENST00000378536.4_1|SKI|SKI-001|protein_coding|cds:overlap        1597    +       2235865 2235925     0       1       60,     0,
#Recibe dos argumentos, el primero el archivo bed que se le pasa y el segundo el numero de replica que es. 
if [ ! -d ./DATA ]
    then
        echo "El directorio ./DATA no existe, lo voy a crear..."
        mkdir ./DATA
    fi    
#Ejecucion de circAnno con la jerarquia. Se eliminan con las segundas lineas aquellos que ya han coincidido.
#Parto del archivo de alineamientos y voy eliminando de él aquellas líneas que no son intergénicos. Utilizo el nuevo archivo, que 
#Solo tiene intergénicos para la siguiente iteración de la jerarquía.
#circanno          bed12                    bed                 resultado.txt
#resultado.txt                                                                          bed         bed
#tRNA
#Copio el archivo de alineamientos
#En estas variables guardo la jerarquía (salvo intergenicos)
FIRST="tRNA"
REST="sncRNA repeats pseudogene lncRNA CDS utr3 utr5 "
LAST="intron"
cp $1 Copy$1
circAnno/bin/circAnno -i 0 ./${FIRST}.bed12 ./$1 | grep -v "intergenic" > ./DATA/${FIRST}$2.txt
#Con esta linea cojo los identificadores de cada pico que no han sido intergenicos y los elimino del archivo de alineamientos
cat ./DATA/${FIRST}$2.txt | cut -f 4 | cut -d ":" -f 1,3  | cut -d ":" -f 1 | sed 's/[.]/\\\\./' \
| tr "\n" "|" | sed 's/|/\\\\|/g' | sed 's/\\\\|$//'| xargs -I{} sed -i '/{}/d' $1 
cp  $1 ./DATA/Rep$21.bed #Copia los que se hayan determinado como intergenicos para el siguiente nivel
rm $1
#Copio el archivo de alineamientos para que tenga el nombre anterior y elimino el que tiene el nombre copia
cp Copy$1 $1
rm Copy$1

NUMERO=1
for SUBGROUP in $REST
do
    circAnno/bin/circAnno -i 0 ./${SUBGROUP}.bed12 ./DATA/Rep$2${NUMERO}.bed | grep -v "intergenic" > ./DATA/${SUBGROUP}$2.txt
    cat ./DATA/${SUBGROUP}$2.txt | cut -f 4 | cut -d ":" -f 1,3  | cut -d ":" -f 1 | sed 's/[.]/\\\\./' \
    | tr "\n" "|" | sed 's/|/\\\\|/g' | sed 's/\\\\|$//'| xargs -I{} sed -i '/{}/d'  ./DATA/Rep$2${NUMERO}.bed
    cp ./DATA/Rep$2${NUMERO}.bed ./DATA/Rep$2$(($NUMERO + 1)).bed
    rm ./DATA/Rep$2${NUMERO}.bed
    NUMERO=$(($NUMERO + 1))
done

#protein coding/intron
circAnno/bin/circAnno -i 0 ./${LAST}.bed12 ./DATA/Rep$2${NUMERO}.bed > ./DATA/tmp.txt
#En este último caso, divido el archivo en intergénicos o no intergénicos. 
grep -v "intergenic" ./DATA/tmp.txt > ./DATA/${LAST}$2.txt
grep "intergenic" ./DATA/tmp.txt > ./DATA/Intergenicos$2.txt
rm ./DATA/Rep$2${NUMERO}.bed ./DATA/tmp.txt

#Jerarquia: tRNA > sncRNA > repeats > pseudogene > lncRNA > protein coding
#Hago la suma de las reads de cada tipo

for HIERARCHY in $FIRST $REST $LAST Intergenicos
do
    awk -v FS="\t" -v HIERARCHY="$HIERARCHY" 'BEGIN{SUM=0};{SUM+=$5};END{print HIERARCHY "\t" SUM}' ./DATA/${HIERARCHY}$2.txt > Resultado${HIERARCHY}.txt
done

#Lo dejo todo en un unico archivo de resultados
cat Resultado* > ./DATA/Resultados$2.txt
rm Resultado*
