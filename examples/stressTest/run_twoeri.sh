#!/bin/bash
# TwoERI_Trans stress test (stress_twoeri.f03): wall time and peak memory (VmHWM)
# per case and size, averaged over the runs, one run at a time.
# The program is built by Build_Tests (make check, or the target below).
# Usage: ./run_twoeri.sh [-M max] [-m min] [-s step] [-n runs] [-t type] [-p threads] [-v]
#   defaults: -M 60 -m 10 -s 10 -n 5 -t nogen; -v prints the progress
#   -p: OpenMP threads (OMP_NUM_THREADS); default: the environment value, or all cores
#   -o: output file (default results_twoeri.txt)
#   type: all (every case), nogen (all but genC and genC_full, the slowest), real,
#   complex, or any single case of stress_twoeri (spaceR, spaceC, spaceC_LR, spinR,
#   spinC, genC, spaceR_full, spaceC_full, spinR_full, spinC_full, genC_full).
#   all/nogen/real/complex: rel.time and rel.mem against the *_full case of the same
#   MO type (same MOs, AO ERI stored without symmetry, same transformed blocks) at the
#   same n; their runs are interleaved (round r runs the reference and every case
#   once), so a change of machine load hits all of them.
#   The final aaaa column is the symmetry class of the aaaa block of the output MO
#   ERIs (after TwoERI_Trans); the AO input is R_8_1234, or StorFull for *_full.
Max_size=60
Min_size=10
Step_size=10
run_number=5
run_type='nogen'
verbose=0
threads_num="${OMP_NUM_THREADS:-}"
# group definitions: "reference:cases" (cases compared with that reference)
space_real='spaceR_full:spaceR'
spin_real='spinR_full:spinR'
space_complex='spaceC_full:spaceC spaceC_LR'
spin_complex='spinC_full:spinC'
gen_complex='genC_full:genC'
# expected rel.time from the n^5 GEMM work of R4_Trans_Full summed over the transformed
# terms of each case (Full term = 1; R8 0.30, R4 0.375, S2P 0.625, R8C/R8LR 0.50,
# InFold only 0.75), divided by the terms of its *_full reference; shown as "(th.)".
# spinR (0.30+0.30+0.375)/3, spinC (0.50+0.50+0.625)/3; genC 8.75/14 (7 blocks, one
# R4_Two_Piece_Trans_Full call each with 2 bra and 2 ket pieces = 2 terms without folds;
# aaaa/bbbb/abab 1.0 each, aabb 1.25, abba/aaab/abbb 1.5 each)
declare -A theory=([spaceR]=0.300 [spaceC]=0.500 [spaceC_LR]=0.500 \
  [spinR]=0.325 [spinC]=0.542 [genC]=0.625)
out_file='results_twoeri.txt'

# Parse command-line options (getopts takes single-letter options only)
while getopts ":M:m:s:n:t:p:o:v" opt; do
  case $opt in
    'M') Max_size="${OPTARG}" ;;
    'm') Min_size="${OPTARG}" ;;
    's') Step_size="${OPTARG}" ;;
    'n') run_number="${OPTARG}" ;;
    't') run_type="${OPTARG}" ;;
    'p') threads_num="${OPTARG}" ;;
    'o') out_file="${OPTARG}" ;;
    'v') verbose=1 ;;
    ':') echo "option -${OPTARG} needs a value"; exit 1 ;;
    *) echo "unknown option -${OPTARG}"; exit 1 ;;
  esac
done
# groups: the "reference:cases" groups to run; empty for a single case (no reference)
case $run_type in
  'all')     groups=("$space_real" "$spin_real" "$space_complex" "$spin_complex" "$gen_complex") ;;
  'nogen')   groups=("$space_real" "$spin_real" "$space_complex" "$spin_complex") ;;
  'real')    groups=("$space_real" "$spin_real") ;;
  'complex') groups=("$space_complex" "$spin_complex" "$gen_complex") ;;
  'spaceR_full'|'spaceC_full'|'spinR_full'|'spinC_full'|'genC_full'| \
  'spaceR'|'spaceC'|'spaceC_LR'|'spinR'|'spinC'|'genC') groups=() ;;
  *) echo "unknown type $run_type"; exit 1 ;;
esac
if [[ -n $threads_num ]]; then
  if ! [[ $threads_num =~ ^[1-9][0-9]*$ ]]; then echo "-p needs a positive integer"; exit 1; fi
  export OMP_NUM_THREADS=$threads_num
  threads_label=$threads_num
else
  unset OMP_NUM_THREADS
  threads_label="$(nproc) (all cores)"
fi
cd "$(dirname "$0")"
trap 'rm -f runs_twoeri_*.tmp' EXIT
sep(){ printf '|-------------+--------+------------+----------+--------+------------+---------+----------------+----------+-------------+-------------|\n' >> $out_file; }
echo START $run_number runs for n = `seq $Min_size $Step_size $Max_size` using $run_type, threads $threads_label > $out_file
printf '\n| %-11s | %6s | %10s | %8s | %6s | %10s | %7s | %14s | %8s | %-11s | %-11s |\n' \
  case n 'mean (s)' sd 'sd %' 'min (s)' 'peak MB' 'rel.time (th.)' rel.mem ref 'final aaaa' >> $out_file
sep
# run_one n case file: one run of stress_twoeri, its result line appended to file;
# stops the script when stress_twoeri rejects the case or fails
run_one(){ local line
  line=$(./stress_twoeri $1 $2 2>&1 | tail -1)
  if [[ $line != *'aaaa='* ]]; then echo "stress_twoeri $1 $2 failed:" >&2; echo "$line" >&2; exit 1; fi
  echo "$line" >> $3; }
# stats n case file             : no reference (rel.time, rel.mem = -)
# stats n case file ref         : this row is the reference (1.000, ref = (ref))
# stats n case file reftime refmem refName : ratios against the reference mean and peak MB
# appends the row to the output file and prints "mean peakMB" on stdout
stats(){ local n=$1 c=$2 f=$3; shift 3
  awk -v c=$c -v n=$n -v reftime="$1" -v refmem="${2:-0}" -v refname="$3" -v out="$out_file" -v th="${theory[$c]:-}" \
    '{t[NR]=$3; s+=$3; if(min==""||$3<min) min=$3; m=$NF; cl=$7; sub(/^aaaa=/,"",cl)}
     END{mean=s/NR; for(i=1;i<=NR;i++) v+=(t[i]-mean)^2; sd=(NR>1)?sqrt(v/(NR-1)):0;
     if(reftime=="ref"){tr="1.000"; mr="1.000"; rf="(ref)"}
     else if(reftime+0>0&&refmem+0>0){tr=sprintf("%.3f",mean/reftime); if(th!="") tr=tr" ("th")"; mr=sprintf("%.3f",m/refmem); rf=refname}
     else {tr="-"; mr="-"; rf="-"}
     sp=(mean>0)?sprintf("%.1f%%",100*sd/mean):"-";
     printf "| %-11s | %6d | %10.4f | %8.4f | %6s | %10.4f | %7d | %14s | %8s | %-11s | %-11s |\n", c, n, mean, sd, sp, min, m, tr, mr, rf, cl >> out;
     printf "%.6f %d\n", mean, m}' $f; }
# avg n case: run_number runs of one case, back to back, no reference
avg(){ : > runs_twoeri_one.tmp
  for r in `seq $run_number`; do run_one $1 $2 runs_twoeri_one.tmp; done
  stats $1 $2 runs_twoeri_one.tmp > /dev/null; }
# group n refCase cases...: interleaved runs (round r runs the reference and every
# case once); then the reference row and each case against it
group(){ local n=$1 ref=$2; shift 2
  local c r
  for c in $ref "$@"; do : > runs_twoeri_$c.tmp; done
  if [[ $verbose == 1 ]]; then printf " $ref+$#:"; fi
  for r in `seq $run_number`; do
    for c in $ref "$@"; do run_one $n $c runs_twoeri_$c.tmp; done
    if [[ $verbose == 1 ]]; then printf " $r"; fi
  done
  read reftime refmem < <(stats $n $ref runs_twoeri_$ref.tmp ref)
  for c in "$@"; do stats $n $c runs_twoeri_$c.tmp $reftime $refmem $ref > /dev/null; done
  for c in $ref "$@"; do rm -f runs_twoeri_$c.tmp; done; }
for n in `seq $Min_size $Step_size $Max_size`; do
  if [[ $verbose == 1 ]]; then printf " $n"; fi
  for g in "${groups[@]}"; do group $n ${g%%:*} ${g#*:}; done
  if [[ ${#groups[@]} -eq 0 ]]; then avg $n $run_type; else sep; fi
  if [[ $verbose == 1 ]]; then printf "\n"; fi
done
printf '\nDONE\n' >> $out_file
