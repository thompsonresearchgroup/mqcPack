#! /bin/bash
# Quick check of the R4_Trans_Full and TwoERI_Trans stress programs (make check): each
# class or case once at n=10; exit status 0 and the expected output class. Timing and
# values are not checked (run_stress.sh and run_twoeri.sh measure the timing).
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
# TwoERI_Trans stress program: exit status 0 and the expected class of the aaaa block
check2(){
  out=$(./stress_twoeri 10 "$1" 2>&1); rc=$?
  line=$(printf '%s\n' "$out" | tail -1)
  if [ $rc -eq 0 ] && printf '%s\n' "$line" | grep -q -- "aaaa=$2 "; then
    echo "PASS  TwoERI $1: aaaa=$2"
  else
    echo "FAIL  TwoERI $1 (rc=$rc): expected aaaa=$2, got: $line"
    nfail=$((nfail+1))
  fi; }
check2 spaceR_full S_1_1234
check2 spaceC_full S_1_1234
check2 spinR_full  S_1_1234
check2 spinC_full  S_1_1234
check2 spaceR      R_8_1234
check2 spaceC      S_4_1234
check2 spaceC_LR   S_2_13_24
check2 spinR       R_8_1234
check2 spinC       S_4_1234
check2 genC        S_4_1234
echo "Tests completed: $nfail failure(s)"
[ $nfail -eq 0 ]
