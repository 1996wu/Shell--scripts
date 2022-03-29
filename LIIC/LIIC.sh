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


template="tmp.inp"
suffix=${template##*.}


natom=$(grep -E "[0-9]*" -o -m 1 simulation.xyz)
ntime=$(grep "" -c simulation.xyz | awk '{print (int ($0 / ("'"${natom}"'" + 2)) ) }')
N_tmp=$(grep -E -c "[a-zA-Z]{1,2}(\s*(-?[0-9]+\.[0-9]*)){3}" "${template}")

if [ -f "${template}" ] && [  -f simulation.xyz ];then
    true
else
    echo "The file '${template}' or 'simulation.xyz' do not exit"
    exit 1
fi 


if [ "${N_tmp}" -ne "${natom}" ];then
  echo "The natom($N_tmp) of template is error"
  exit 1
fi

k=0
l=0 # read l time geom from 'simulation.xyz'
while [ "${ntime}" -gt "${l}" ]
do 
	new_file="${k}.${suffix}" 
    replace_coordinate "$natom" "${l}" "${template}" "${new_file}"
    #g16 "${new_file}" # running software
    sh bdf_template.sh "${new_file}" > "${k}".log
    if [ "$((l+1))" -gt "${ntime}"  ];then 
        cp "${k}.scforb" "$((k+1)).scforb"
    fi 
	l=$((l+1))
	k=$((k+1))
done 
