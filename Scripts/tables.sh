mkdir -p ./Docs/Tables/

for f in $(ls ./Data/Exported/*.csv); do
    rm ./Docs/Tables/$(basename ${f}).tex
    python3 ./tools/tably/tably.py ${f} -o ./Docs/Tables/$(basename ${f}).tex -s ";"
    sed -i '' \
        -e 's/\[htb\]/\[H\]/g' \
        -e 's/ *\\bottomrule//g' \
        -e 's/ *\\toprule//g' \
        -e 's/ *\\midrule//g' \
        -e 's/ccc/|c|c|c|/g' \
        -e 's/\\\\/\\\\\\hline/g' \
        -e 's/labels/\\hline labels/g' \
        -e 's/\\\_/ /g' \
        ./Docs/Tables/$(basename ${f}).tex
done