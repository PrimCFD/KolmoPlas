#!/usr/bin/env bash
set -euo pipefail
# Build (optional) and run *perf* tests using build.sh

usage() {
  cat <<'EOF'
Usage: scripts/run_perf_tests.sh [options]

Build (optional) and run performance tests via CTest (label: perf), writing
JUnit XML under $BUILD_DIR/test-reports/perf/. Sets a laptop/cluster-friendly
MPI environment via scripts/mpi_env.sh when available.

Options:
  -h, --help           Show this help message and exit.

Environment:
  BUILD_DIR            Build directory (default: build-perf).
  SKIP_BUILD           0|1 to skip the build step (default: 0).
  CTEST_PARALLEL_LEVEL Number of tests to run in parallel (default: CPU count).
  CTEST_TIMEOUT        Per-test timeout in seconds (default: 900).
  REPORT_DIR           Output directory for JUnit XML
                       (default: $BUILD_DIR/test-reports/perf).
  MPI_MODE             auto|emulate|cluster mode for mpi_env.sh (default: auto).
  MPI_STRICT_PE        0|1 to enforce ranks×PE ≤ cores (default: 1).
  OMP_NUM_THREADS      Threads per MPI rank (default: 1 here).
  ENABLE_CUDA          ON/OFF/AUTO for CUDA perf runs (default: AUTO logic).
  USE_CUDA_UM          ON/OFF for CUDA unified memory (default: OFF or ON if AUTO chooses CUDA).
  PETSC_CONFIGURE_OPTS Extra PETSc configure flags for vendored PETSc.
  EXTRA_CMAKE_ARGS     Extra CMake args (tests toggles are appended automatically).
  EXTRA_CMAKE_ARGS_USER User-provided CMake args appended after toggles.

Examples:
  scripts/run_perf_tests.sh
  SKIP_BUILD=1 scripts/run_perf_tests.sh
  BUILD_DIR=build-perf CTEST_PARALLEL_LEVEL=4 scripts/run_perf_tests.sh
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
  shift
done


# Set MPI env

if [[ -f "scripts/mpi_env.sh" ]]; then
  # Override with MPI_MODE=cluster|emulate when needed.
  # shellcheck source=/dev/null
  source scripts/mpi_env.sh "${MPI_MODE:-auto}"
  # Perf-critical: turn on strict PE by default (ranks × PE must fit cores)
  export MPI_STRICT_PE="${MPI_STRICT_PE:-1}"
  export OMP_NUM_THREADS=1 # limit thread number/MPI
  export MPI_TEST_NP=8
fi

PETSC_OPTIONS="-ksp_converged_reason \
-ksp_monitor_short"

# Uncomment for PETSC linear solve logging
# export PETSC_OPTIONS

BUILD_DIR=${BUILD_DIR:-build-perf}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/perf"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  # Suggested PETSc configure for vendored builds (ignored if USE_SYSTEM_PETSC=ON)
 : "${OPT_LEVEL:=3}"
 : "${ARCH_FLAGS:=-march=native}"
 : "${PETSC_CONFIGURE_OPTS:=--with-debugging=0 --with-mpi=1 \
   COPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   CXXOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} \
   FOPTFLAGS='-O'${OPT_LEVEL}' '${ARCH_FLAGS} }"
 export PETSC_CONFIGURE_OPTS
  # Auto-detect CUDA unless user forced it
  if [[ -z "${ENABLE_CUDA:-}" || "${ENABLE_CUDA}" == "AUTO" ]]; then
    if command -v nvcc >/dev/null || [[ -d "${CUDAToolkit_ROOT:-/usr/local/cuda}" ]]; then
      ENABLE_CUDA=ON
      : "${USE_CUDA_UM:=ON}"
    else
      ENABLE_CUDA=OFF
    fi
  fi

  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  ENABLE_CUDA="${ENABLE_CUDA}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=OFF -DENABLE_TESTS_PERF=ON -DENABLE_TESTS_REGRESSION=OFF" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi


mkdir -p "${REPORT_DIR}"

if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-perf.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

ctest -V --test-dir "${BUILD_DIR}" \
      -L perf \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}" \

# Post-run fallback: only if no XML was produced by CTest/Catch2
shopt -s nullglob
perf_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#perf_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="perf" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-perf.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-perf.xml"
fi
