AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
AM_LDFLAGS = $(MQC_LDFLAGS)
AUTOMAKE_OPTIONS = subdir-objects
TESTS_ENVIRONMENT = abs_top_srcdir='$(abs_top_srcdir)' \
  abs_top_builddir='$(abs_top_builddir)' \
  FC='$(FC)' MQC_TEST_FCFLAGS='$(AM_FCFLAGS) $(FCFLAGS)'
CLEANFILES = ../algebra2/outfile ../fullci/outfile ../hartreefock/outfile \
  ../../unitTests/src/outfile unitTest07-numerical.faf \
  ../CreateDataList/outfile

TESTS = ../algebra2/testit.sh unitTest01 unitTest02 unitTestIterative unitTest08 unitTest09 \
  ../../unitTests/src/testit.sh ../../unitTests/src/unitTest_NEF_collapse_expand \
  ../../unitTests/src/unitTest_NEF_set_algebra \
  ../../unitTests/src/unitTest_NEF_r4_full_transform \
  ../../unitTests/src/unitTest_NEF_twoeris_r4_transform \
  ../../unitTests/src/unitTest_NEF_empty_space_transform \
  ../../unitTests/src/unitTest_NEF_merge_unmerge_roundtrip \
  ../../unitTests/src/unitTest_symmetry ../stressTest/testit.sh
check_PROGRAMS = ../algebra2/algebra_fun unitTest01 unitTest02 unitTestIterative unitTest08 \
  unitTest09 ../../unitTests/src/unitTest_twoERIs ../../unitTests/src/unitTest_NEF_collapse_expand \
  ../../unitTests/src/unitTest_NEF_set_algebra \
  ../../unitTests/src/unitTest_NEF_r4_full_transform \
  ../../unitTests/src/unitTest_NEF_twoeris_r4_transform \
  ../../unitTests/src/unitTest_NEF_empty_space_transform \
  ../../unitTests/src/unitTest_NEF_merge_unmerge_roundtrip \
  ../../unitTests/src/unitTest_symmetry ../stressTest/stress_trans \
  ../stressTest/stress_twoeri

if WITH_GAUOPEN
TESTS += ../fullci/testit.sh ../hartreefock/testit.sh \
  unitTest07
check_PROGRAMS += ../fullci/fullci ../hartreefock/hartreefock unitTest07
# Runs Gaussian when GAU_BINARY is set (mqc_install exports it); skips otherwise.
TESTS += test_gaussian_createdatalist.sh
check_PROGRAMS += ../CreateDataList/CheckInput
if GAUOPEN_FRONTIER
TESTS += unitTest03 unitTest04
check_PROGRAMS += unitTest03 unitTest04
else
TESTS += test_public_gauopen_character.sh
check_PROGRAMS += unitTest03
endif
else
TESTS += unitTest05 test_no_gauopen_matrixfile.sh
check_PROGRAMS += unitTest05 unitTest06
endif

___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.F03
___fullci_fullci_SOURCES = ../fullci/fullci.F03
___hartreefock_hartreefock_SOURCES = ../hartreefock/hartreefock.F03
___CreateDataList_CheckInput_SOURCES = ../CreateDataList/CheckInput.F03
unitTest01_SOURCES = ../../unitTests/src/unitTest01.f03
unitTest02_SOURCES = ../../unitTests/src/unitTest02.f03
unitTest03_SOURCES = ../../unitTests/src/unitTest03.f03
unitTest04_SOURCES = ../../unitTests/src/unitTest04.f03
unitTest05_SOURCES = ../../unitTests/src/unitTest05.f03
unitTest06_SOURCES = ../../unitTests/src/unitTest06.f03
unitTest07_SOURCES = ../../unitTests/src/unitTest07.f03
unitTestIterative_SOURCES = ../../unitTests/src/unitTestIterative.f03
unitTest08_SOURCES = ../../unitTests/src/unitTest08.f03
unitTest09_SOURCES = ../../unitTests/src/unitTest09.f03
______unitTests_src_unitTest_twoERIs_SOURCES = \
  ../../unitTests/src/unitTest_twoERIs.f03
______unitTests_src_unitTest_NEF_collapse_expand_SOURCES = \
  ../../unitTests/src/unitTest_NEF_collapse_expand.F03
______unitTests_src_unitTest_NEF_set_algebra_SOURCES = \
  ../../unitTests/src/unitTest_NEF_set_algebra.f03
______unitTests_src_unitTest_NEF_r4_full_transform_SOURCES = \
  ../../unitTests/src/unitTest_NEF_r4_full_transform.f03
______unitTests_src_unitTest_NEF_twoeris_r4_transform_SOURCES = \
  ../../unitTests/src/unitTest_NEF_twoeris_r4_transform.f03
______unitTests_src_unitTest_NEF_empty_space_transform_SOURCES = \
  ../../unitTests/src/unitTest_NEF_empty_space_transform.f03
______unitTests_src_unitTest_NEF_merge_unmerge_roundtrip_SOURCES = \
  ../../unitTests/src/unitTest_NEF_merge_unmerge_roundtrip.f03
______unitTests_src_unitTest_symmetry_SOURCES = \
  ../../unitTests/src/unitTest_symmetry.f03
___stressTest_stress_trans_SOURCES = ../stressTest/stress_trans.f03
___stressTest_stress_twoeri_SOURCES = ../stressTest/stress_twoeri.f03
EXTRA_DIST = test_no_gauopen_matrixfile.sh \
  test_public_gauopen_character.sh ../algebra2/testit.sh \
  ../algebra2/OUTPUT/out ../fullci/testit.sh ../fullci/OUTPUT/out \
  ../hartreefock/testit.sh ../hartreefock/OUTPUT/out \
  ../data/MatrixFile/rhf_h2-sto3g.mat ../../unitTests/src/testit.sh \
  ../../unitTests/src/OUTPUT/out test_gaussian_createdatalist.sh \
  ../CreateDataList/testit_runGau.sh ../CreateDataList/OUTPUT/out_runGau \
  ../data/Gaussian_input/test.com ../data/Gaussian_input/test1071.com \
  ../data/Gaussian_input/test1129.com ../data/Gaussian_input/test1132.com \
  ../stressTest/run_stress.sh \
  ../stressTest/run_twoeri.sh ../stressTest/testit.sh

LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}

BUILD_PROGRAM_OBJECTS = ../algebra2/algebra_fun.$(OBJEXT)
if WITH_GAUOPEN
BUILD_PROGRAM_OBJECTS += ../fullci/fullci.$(OBJEXT) \
  ../hartreefock/hartreefock.$(OBJEXT) unitTest07.$(OBJEXT)
endif

build: $(BUILD_PROGRAM_OBJECTS)

# Gaussian-independent quadrature rules, partitions, and owning-grid contracts.
TESTS += unitTest10 unitTest11 unitTest12 test_quadrature_invalid.sh test_quadrature_private.sh
check_PROGRAMS += unitTest10 unitTest11 unitTest12
unitTest10_SOURCES = ../../unitTests/src/unitTest10.f03
unitTest11_SOURCES = ../../unitTests/src/unitTest11.f03
unitTest12_SOURCES = ../../unitTests/src/unitTest12.f03
EXTRA_DIST += test_quadrature_invalid.sh test_quadrature_private.sh
