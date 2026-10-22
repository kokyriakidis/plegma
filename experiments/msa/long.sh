cd ~/plegma-bench/prof
B=~/plegma-bench/msa/bin; K=~/karan/build/cli/karan
run() { name=$1; L=$2; shift 2; (ulimit -v 67108864; /usr/bin/time -f "%e s %M KB" timeout 3600 "$@" > out.$name.$L.fa 2> err.$name.$L; echo "$name $L exit=$? $(tail -1 err.$name.$L) ed=$($K ref_$L.fa out.$name.$L.fa 2>/dev/null | tail -1 | cut -f4)"); }
for f in "10000 50 0.1" "100000 10 0.05" "1000000 10 0.01"; do set -- $f; R=reads_$1_$2_$3.fa
  run plegma_S $1 ~/plegma-dev/build/plegma -S $R
  run theseus $1 $B/theseus_msa -m 0 -x 6 -o 4 -e 3 -t 2 -s $R -f /dev/stdout
  run abpoa $1 $B/abpoa $R
done
