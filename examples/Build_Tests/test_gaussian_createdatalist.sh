#! /bin/bash
#
# Run the CreateDataList example against a live Gaussian installation.
#
# mqc_install exports GAU_BINARY when the installer chooses to run the
# Gaussian tests. Without it, or without the binary in PATH, the test is
# skipped (automake exit status 77).
#
if [ -z "${GAU_BINARY}" ]; then
    echo "GAU_BINARY not set; skipping Gaussian CreateDataList test."
    exit 77
fi
gau_path="`command -v ${GAU_BINARY}`"
if [ ! -x "${gau_path}" ]; then
    echo "${GAU_BINARY} not found in PATH; skipping Gaussian CreateDataList test."
    exit 77
fi
#
# The Gaussian inputs call the external helpers rdmat-ext and testunf, which
# ship as C-shell scripts in the Gaussian bsd directory.
#
gau_bsd="`dirname ${gau_path}`/bsd"
if [ -d "${gau_bsd}" ]; then
    export PATH="${PATH}:${gau_bsd}"
fi
for helper in rdmat-ext testunf; do
    if ! command -v ${helper} > /dev/null 2>&1; then
        echo "Gaussian helper ${helper} not found in PATH or ${gau_bsd}."
        exit 1
    fi
done
if ! command -v csh > /dev/null 2>&1; then
    echo "csh is required to run the Gaussian helpers rdmat-ext and testunf."
    exit 1
fi

bash ../CreateDataList/testit_runGau.sh
