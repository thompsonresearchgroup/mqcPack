#! /bin/bash
# Quick check of the R4_Trans_Full stress program (make check): each class once at
# n=10; exit status 0 and the expected input -> output class. Timing and values are
# not checked (run_stress.sh measures the timing).
: "${abs_top_builddir:=$(cd ../.. && pwd)}"
cd "${abs_top_builddir}/examples/stressTest"
nfail=0
check(){
  out=$(./stress_trans 10 "$1" 2>&1); rc=$?
  line=$(printf '%s\n' "$out" | tail -1)
  if [ $rc -eq 0 ] && printf '%s\n' "$line" | grep -q -- " $2 ->$3 "; then
    echo "PASS  $1: $2 -> $3"
  else
    echo "FAIL  $1 (rc=$rc): expected $2 -> $3, got: $line"
    nfail=$((nfail+1))
  fi; }
check Full  unknown   S_1_1234
check FullC unknown   S_1_1234
check S2    S_2_13_24 S_2_13_24
check S4    S_4_1234  S_4_1234
check R8    R_8_1234  R_8_1234
check Klein S_4_1234  S_4_1234
check S2C    S_2_13_24      S_2_13_24
check S2P    S_2_12_34      S_2_12_34
check KleinC S_4_1234_plain S_4_1234_plain
check R4     R_4_12_34      R_4_12_34
check R2     R_2_12_3_4     R_2_12_3_4
check R8C    R_8_1234       S_4_1234
check R8LR   R_8_1234       S_2_13_24
echo "Tests completed: $nfail failure(s)"
[ $nfail -eq 0 ]
