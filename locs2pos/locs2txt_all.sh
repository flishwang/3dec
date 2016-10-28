echo "Usage:"
echo "  ./locs2txt_all.sh dirA dirB"
echo "      Convert all locs files in dir1 into _pos.txt and output them into dir2"
input=${1:-.}
output=${2:-.}
mkdir -p $output
for i in `ls $input/*.locs`;do 
t=`basename $i .locs`
tt="./locs2txt $input/${t}.locs $output/${t}_pos.txt"
echo $tt
$tt
done