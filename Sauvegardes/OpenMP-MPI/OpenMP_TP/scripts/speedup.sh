SCRIPTDIR=$(cd "$( dirname "$0" )" && pwd)

fileList=$(grep -l CPU *.res)

for file in $fileList
do
    python $SCRIPTDIR/speedup.py $file
    echo
done
