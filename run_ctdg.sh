readarray pfams < $1
CPUS=$2
for acc in ${pfams[@]}
do
    test -d CTDG_out/$acc && echo $acc Exists || python ../ctdg_finder/ctdg_finder.py -c $CPUS -n $acc -p $acc -S 1000 -d ../ctdg_db/mammals_18
done

