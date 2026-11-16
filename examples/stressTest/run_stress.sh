#!/bin/bash
# R4_Trans_Full stress test (stress_trans.f03): wall time and peak memory (VmHWM)
# per class and size, averaged over the runs, one run at a time.
# The program is built by Build_Tests (make check, or the target below).
# Usage: ./run_stress.sh [-M max] [-m min] [-s step] [-n runs] [-t type] [-p threads] [-v]
#   defaults: -M 60 -m 10 -s 10 -n 5 -t all; -v prints the progress
#   -p: OpenMP threads (OMP_NUM_THREADS); default: the environment value, or all cores
#   type: all, real, complex, or any single class of stress_trans (Full, FullC, S2,
#   S2C, S2P, S4, Klein, KleinC, R8, R4, R2, R8C, R8LR; f = Full, fc = FullC).
#   Output: results_avg.txt
#   all/real/complex: rel.time and rel.mem against Full (real classes) or FullC
#   (complex classes) at the same n; their runs are interleaved (round r runs the
#   reference and every class once), so a change of machine load hits all of them
Max_size=60
Min_size=40
Step_size=10
run_number=5
run_type='all'
verbose=0
threads_num="${OMP_NUM_THREADS:-}"
# ordered slowest to fastest by the expected n^5 work (theory; ties by measured n=40)
real_classes='S2 R2 Klein R8LR R4 R8'   # rel. to Full:  0.75 0.75 0.50 0.50 0.375 0.30
complex_classes='S2C S2P KleinC R8C S4' # rel. to FullC: 0.75 0.625 0.50 0.50 0.50
# expected rel.time from the n^5 GEMM work of R4_Trans_Full (Full = 4n^5), shown as
# "(th.)" next to the measured rel.time
declare -A theory=([S2]=0.750 [R2]=0.750 [Klein]=0.500 [R8LR]=0.500 [R4]=0.375 [R8]=0.300 \
  [S2C]=0.750 [S2P]=0.625 [KleinC]=0.500 [R8C]=0.500 [S4]=0.500)

# Parse command-line options (getopts takes single-letter options only)
while getopts ":M:m:s:n:t:p:v" opt; do
  case $opt in
    'M') Max_size="${OPTARG}"
      ;;
    'm') Min_size="${OPTARG}"
      ;;
    's') Step_size="${OPTARG}"
      ;;
    'n') run_number="${OPTARG}"
      ;;
    't') run_type="${OPTARG}"
      ;;
    'p') threads_num="${OPTARG}"
      ;;
    'v') verbose=1
      ;;
    ':') echo "option -${OPTARG} needs a value"; exit 1
      ;;
    *) echo "unknown option -${OPTARG}"; exit 1
      ;;
  esac
done
if [[ -n $threads_num ]]; then
  if ! [[ $threads_num =~ ^[1-9][0-9]*$ ]]; then echo "-p needs a positive integer"; exit 1; fi
  export OMP_NUM_THREADS=$threads_num
  threads_label=$threads_num
else
  unset OMP_NUM_THREADS
  threads_label="$(nproc) (all cores)"
fi
cd "$(dirname "$0")"
trap 'rm -f runs_*.tmp' EXIT
echo  START $run_number runs for n = `seq $Min_size $Step_size $Max_size` using $run_type, threads $threads_label  > results_avg.txt
printf '\n| %-8s | %6s | %10s | %8s | %6s | %10s | %7s | %14s | %8s | %-5s |\n' \
  class n 'mean (s)' sd 'sd %' 'min (s)' 'peak MB' 'rel.time (th.)' rel.mem ref >> results_avg.txt
printf '|----------+--------+------------+----------+--------+------------+---------+----------------+----------+-------|\n' >> results_avg.txt
# run_one n class file: one run of stress_trans, its result line appended to file;
# stops the script when stress_trans rejects the class or fails to pack it
run_one(){ local line
  line=$(./stress_trans $1 $2 2>&1 | tail -1)
  if [[ $line != *'->'* ]]; then echo "stress_trans $1 $2 failed:" >&2; echo "$line" >&2; exit 1; fi
  echo "$line" >> $3; }
# stats n class file             : no reference (rel.time, rel.mem = -)
# stats n class file ref         : this row is the reference (1.000, ref = (ref))
# stats n class file reftime refmem refName : ratios against the reference mean and peak MB
# appends the row to results_avg.txt and prints "mean peakMB" on stdout
stats(){ local n=$1 c=$2 f=$3; shift 3
  awk -v c=$c -v n=$n -v reftime="$1" -v refmem="${2:-0}" -v refname="$3" -v th="${theory[$c]:-}" \
    '{t[NR]=$3; s+=$3; if(min==""||$3<min) min=$3; m=$NF}
     END{mean=s/NR; for(i=1;i<=NR;i++) v+=(t[i]-mean)^2; sd=(NR>1)?sqrt(v/(NR-1)):0;
     if(reftime=="ref"){tr="1.000"; mr="1.000"; rf="(ref)"}
     else if(reftime+0>0&&refmem+0>0){tr=sprintf("%.3f",mean/reftime); if(th!="") tr=tr" ("th")"; mr=sprintf("%.3f",m/refmem); rf=refname}
     else {tr="-"; mr="-"; rf="-"}
     sp=(mean>0)?sprintf("%.1f%%",100*sd/mean):"-";
     printf "| %-8s | %6d | %10.4f | %8.4f | %6s | %10.4f | %7d | %14s | %8s | %-5s |\n", c, n, mean, sd, sp, min, m, tr, mr, rf >> "results_avg.txt";
     printf "%.6f %d\n", mean, m}' $f; }
# avg n class [ref | reftime refmem refName]: run_number runs of one class, back to back
avg(){ : > runs.tmp
  for r in `seq $run_number`; do run_one $1 $2 runs.tmp; done
  stats $1 $2 runs.tmp "${@:3}"; }
# group n refClass classes...: interleaved runs (round r runs the reference and every
# class once), so a change of machine load hits all of them; then the reference row
# and each class against it
group(){ local n=$1 ref=$2; shift 2
  local c r
  for c in $ref "$@"; do : > runs_$c.tmp; done
  if [[ $verbose == 1 ]]; then printf " $ref+$#:"; fi
  for r in `seq $run_number`; do
    for c in $ref "$@"; do run_one $n $c runs_$c.tmp; done
    if [[ $verbose == 1 ]]; then printf " $r"; fi
  done
  read reftime refmem < <(stats $n $ref runs_$ref.tmp ref)
  for c in "$@"; do stats $n $c runs_$c.tmp $reftime $refmem $ref > /dev/null; done
  for c in $ref "$@"; do rm -f runs_$c.tmp; done; }
sep(){ printf '|----------+--------+------------+----------+--------+------------+---------+----------------+----------+-------|\n' >> results_avg.txt; }
case $run_type in
  'all')
    for n in `seq $Min_size $Step_size $Max_size`; do
      if [[ $verbose == 1 ]]; then printf " $n"; fi
      group $n Full $real_classes
      group $n FullC $complex_classes
      if [[ $verbose == 1 ]]; then printf "\n"; fi
      sep
    done
    ;;
  'real')
    for n in `seq $Min_size $Step_size $Max_size`; do
      if [[ $verbose == 1 ]]; then printf " $n"; fi
      group $n Full $real_classes
      if [[ $verbose == 1 ]]; then printf "\n"; fi
      sep
    done
    ;;
  'complex')
    for n in `seq $Min_size $Step_size $Max_size`; do
      if [[ $verbose == 1 ]]; then printf " $n"; fi
      group $n FullC $complex_classes
      if [[ $verbose == 1 ]]; then printf "\n"; fi
      sep
    done
    ;;
  'f')
    for n in `seq $Min_size $Step_size $Max_size`;
    do avg $n 'Full' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi;
    ;;
  'S2')
    for n in `seq $Min_size $Step_size $Max_size`;
      do avg $n 'S2' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi;
    ;;
  'S4')
    for n in `seq $Min_size $Step_size $Max_size`;
      do avg $n 'S4' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi;
    ;;
  'R8')
    for n in `seq $Min_size $Step_size $Max_size`;
      do avg $n 'R8' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi
    ;;
  'fc')
    for n in `seq $Min_size $Step_size $Max_size`;
      do avg $n 'FullC' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi
    ;;
  'Klein')
    for n in `seq $Min_size $Step_size $Max_size`;
      do avg $n 'Klein' > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi;
    done
    if [[ $verbose == 1 ]]; then printf "\n";  fi
    ;;
  *)
    for n in `seq $Min_size $Step_size $Max_size`; do
      avg $n $run_type > /dev/null
      if [[ $verbose == 1 ]]; then printf " $n"; fi
    done
    if [[ $verbose == 1 ]]; then printf "\n"; fi
    ;;
esac
printf '\nDONE\n' >> results_avg.txt
