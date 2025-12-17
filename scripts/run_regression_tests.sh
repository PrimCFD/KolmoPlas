#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *regression/integration* tests using build.sh

usage() {
  cat <<'EOF'
Usage: scripts/run_regression_tests.sh [options]

Configure and build (optional) with MPI enabled and run MPI regression tests
via CTest (label: regression), using scripts/mpi_env.sh for consistent MPI
launcher setup on laptops and clusters.

Options:
  -h, --help           Show this help message and exit.

Environment:
  BUILD_DIR            Build directory (default: build-regression).
  SKIP_BUILD           0|1 to skip the build step (default: 0).
  CTEST_PARALLEL_LEVEL Number of tests to run in parallel (default: CPU count).
  CTEST_TIMEOUT        Per-test timeout in seconds (default: 900).
  REPORT_DIR           Output directory for JUnit XML
                       (default: $BUILD_DIR/test-reports/regression).
  CMAKE_BUILD_TYPE     Build type (default: Release).
  MPI_MODE             auto|emulate|cluster for mpi_env.sh (default: auto).
  CTEST_NAME_REGEX     Optional CTest name regex to filter which tests run
                       (used by kolmoplas examples run).
  EXTRA_CMAKE_ARGS     Extra CMake args, including MPI/test toggles.
  EXTRA_CMAKE_ARGS_USER Extra user-provided CMake args.

Examples:
  scripts/run_regression_tests.sh
  SKIP_BUILD=1 scripts/run_regression_tests.sh
  CTEST_NAME_REGEX='regression::fluids::tgv128::run' scripts/run_regression_tests.sh
EOF
}


# Set MPI env
if [[ -f "scripts/mpi_env.sh" ]]; then
  # Override with MPI_MODE=cluster|emulate when needed.
  # shellcheck source=/dev/null
  source scripts/mpi_env.sh "${MPI_MODE:-auto}"
  # Perf-critical: turn on strict PE by default (ranks × PE must fit cores)
  export MPI_STRICT_PE="${MPI_STRICT_PE:-1}"
  export OMP_NUM_THREADS=1 # limit thread number for small runs (32^3)
fi

BUILD_DIR=${BUILD_DIR:-build-regression}
SKIP_BUILD=${SKIP_BUILD:-0}
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-150000}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/regression"}
LABEL_RE=${LABEL_RE:-regression}
CTEST_NAME_REGEX=${CTEST_NAME_REGEX:-}

PETSC_OPTIONS="-ksp_converged_reason \
-ksp_monitor_short \
-ksp_view \
-pc_view \
-error_output_stdout \
-on_error_abort"


# Uncomment for PETSC linear solve logging
# export PETSC_OPTIONS

if [[ "${SKIP_BUILD}" != "1" ]]; then
  # Suggested PETSc configure for vendored builds (ignored if USE_SYSTEM_PETSC=ON)
 : "${OPT_LEVEL:=3}"
 : "${ARCH_FLAGS:=-march=native}"
 : "${PETSC_CONFIGURE_OPTS:=--with-debugging=0 --with-mpi=1 \
   COPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   CXXOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   FOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} }"
  export PETSC_CONFIGURE_OPTS
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_PERF=OFF -DENABLE_TESTS_REGRESSION=ON" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi

mkdir -p "${REPORT_DIR}"

# Ask ctest to emit JUnit if it supports it
if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-regression.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

# Choose filter mode: by label (default) or by test name regex (if CTEST_NAME_REGEX set)
CTEST_FILTER_OPTS=()
if [[ -n "${CTEST_NAME_REGEX}" ]]; then
  CTEST_FILTER_OPTS=( -R "${CTEST_NAME_REGEX}" )
else
  CTEST_FILTER_OPTS=( --label-regex "${LABEL_RE}" )
fi

ctest -V --test-dir "${BUILD_DIR}" \
      "${CTEST_FILTER_OPTS[@]}" \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}"

# Post-run fallback: if nothing produced JUnit, create a tiny placeholder
shopt -s nullglob
reg_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#reg_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="regression" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-regression.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-regression.xml"
fi
