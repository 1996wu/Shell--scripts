#! /bin/bash 
function replace_coordinate {
    awk 'BEGIN{var="'"$1"'"; var_1="'"$2"'"}
    {
        if ( (var + 2)* var_1 +3  <= NR && NR<= (var + 2) * (var_1 + 1) )
        print $0 
    }' simulation.xyz  > new_coordinate
    awk 'BEGIN {
        if (getline < "new_coordinate" > 0)
        {l=1}
        close("new_coordinate")
        while (getline line < "new_coordinate" >0 )
        {
            value[l] = line
            l++;
        }
    }
    {
        l = 1 
        while (l <=NR){
            if ($0 ~ /[a-zA-Z]{1,2}(\s*(-?[0-9]+\.[0-9]*)){3}/){
                i++
            }
            gsub(/[a-zA-Z]{1,2}(\s*(-?[0-9]+\.[0-9]*)){3}/,value[i])
            print $0
            l++
            next
        } 
    }' "$3" > "$4"
    sed -i 's/^M//g' "$4"
    rm new_coordinate
}
#  $1 natom  $2 ntime $3 old_file $4 new_file

if [ -f tmp.gjf  ] && [  -f simulation.xyz ];then
    true
else
    echo "The file 'tmp.gjf' or 'simulation.xyz' do not exit"
    exit 1
fi 
l=0 
natom=$(head  -n1  simulation.xyz)
ntime=$(grep "" -c simulation.xyz | awk '{print (int ($0 / ("'"${natom}"'" + 2)) ) }')
while [ "${ntime}" -gt "${l}" ]
do 
    replace_coordinate "$natom"  "${l}" tmp.gjf  "${l}".gjf
    g16 "${l}".gjf  
    l=$((l+1))
done 
