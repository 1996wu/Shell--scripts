#/bin/bash 
N=3 #n states
l=$((N+2)) 
A=$(echo ${PWD##*/})
eV=27.2114

for i in $(seq 0 ${N}) 
do
  if [ "${i}" == "0" ];then
    echo -e  "Liic_path \t S${i}  "  >> tmp_1
  else 
    echo  "S${i} " >>tmp_1
  fi
done 
cat tmp_1 | xargs > tmp_2


for inf in $* 
do
filename=$(echo ${inf%.*}) 
suffix=$(echo ${inf##*.}) 
if [ "${suffix}" == "out" -o "${suffix}" == "log" ] 2> /dev/null ;then  
  for i in $(seq 0 ${N}) 
  do
  if  [ "${i}" == "0" ] ; then 
    grep "SCF Done" ${inf} | awk  '{print $5}' >tmp_3
  else 
    grep "Excited State   ${i}" ${inf} | awk  '{print $5}' >> tmp_3 
  fi 
  done 
  cat tmp_3 | xargs | awk -v var=${eV} '{
  for (i=1; i<=NF; i=i+1)
    if ( i==1) 
      printf("%14.5f", $1*var)
    else
      printf ("%14.5f", $i+$1*var)} 
  END{print "\t"}'  | sed "1s/^/&${filename}   /" >>tmp_4 
else 
  echo "suffix(${suffix}) is error" 
fi
done

#min=$(awk '{print $2}' tmp_4 |sort -n -b | head -n1 )
#paste $(echo "$(seq 0  ${N})" |xargs) | sort -n -b  >tmp_3 
min=-89742.95384872193 #!!!!!! Unit eV 


echo "${min}"
sort -n -b tmp_4  | awk  -v  var=${min} '{ 
  for (i=1; i<=NF; i=i+1) 
  if (i==1 )  
    printf ("%s",$i)
  else
    printf("%14.5f ",$i-var)} 
  END{print "\n"} ' | xargs -n ${l}  > tmp_5 
cat tmp_2 tmp_5 | column -t > Energy_${A}

rm tmp_* [0-9] 2> /dev/null
 
