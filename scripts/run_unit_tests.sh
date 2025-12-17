#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/run_unit_tests.sh [options]

Build (optional) and run unit tests via CTest, writing JUnit XML under
$BUILD_DIR/test-reports/unit/.

Options:
  -h, --help           Show this help message and exit.

Environment:
  BUILD_DIR            Build directory (default: build-unit).
  SKIP_BUILD           0|1 to skip the build step (default: 0).
  CTEST_PARALLEL_LEVEL Number of tests to run in parallel (default: CPU count).
  CTEST_TIMEOUT        Per-test timeout in seconds (default: 900).
  REPORT_DIR           Output directory for JUnit XML
                       (default: $BUILD_DIR/test-reports/unit).
  CMAKE_BUILD_TYPE     Build type for the unit-test build (default: Debug).
  ENABLE_CUDA          ON/OFF (default: OFF for unit tests).
  USE_CUDA_UM          ON/OFF (default: OFF).
  EXTRA_CMAKE_ARGS     Extra flags appended to CMake configure.
  EXTRA_CMAKE_ARGS_USER Extra user-provided CMake args.

Examples:
  scripts/run_unit_tests.sh
  SKIP_BUILD=1 scripts/run_unit_tests.sh
  BUILD_DIR=build-debug CTEST_PARALLEL_LEVEL=8 scripts/run_unit_tests.sh
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 2 ;;
  esac
  shift
done


BUILD_DIR=${BUILD_DIR:-build-unit}
SKIP_BUILD=${SKIP_BUILD:-0}
CTEST_PARALLEL_LEVEL=${CTEST_PARALLEL_LEVEL:-$(command -v nproc >/dev/null && nproc || echo 2)}
CTEST_TIMEOUT=${CTEST_TIMEOUT:-900}
REPORT_DIR=${REPORT_DIR:-"${BUILD_DIR}/test-reports/unit"}

if [[ "${SKIP_BUILD}" != "1" ]]; then
  BUILD_DIR="${BUILD_DIR}" \
  CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Debug}" \
  ENABLE_CUDA="${ENABLE_CUDA:-OFF}" \
  USE_CUDA_UM="${USE_CUDA_UM:-OFF}" \
  BUILD_TESTS=ON \
  TEST_TOGGLES="-DENABLE_TESTS_UNIT=ON -DENABLE_TESTS_PERF=OFF -DENABLE_TESTS_REGRESSION=OFF" \
  EXTRA_CMAKE_ARGS="${EXTRA_CMAKE_ARGS:+$EXTRA_CMAKE_ARGS }${TEST_TOGGLES}${EXTRA_CMAKE_ARGS_USER:+ ${EXTRA_CMAKE_ARGS_USER}}" \
  scripts/build.sh
else
  [[ -d "${BUILD_DIR}" ]] || { echo "❌ BUILD_DIR ${BUILD_DIR} not found; remove SKIP_BUILD=1"; exit 1; }
fi


mkdir -p "${REPORT_DIR}"

if ctest --help 2>/dev/null | grep -q -- '--output-junit'; then
  CTEST_JUNIT_OPTS=( --output-junit "${REPORT_DIR}/ctest-unit.xml" )
else
  CTEST_JUNIT_OPTS=()
fi

ctest --test-dir "${BUILD_DIR}" \
      -L unit \
      -j "${CTEST_PARALLEL_LEVEL}" \
      --timeout "${CTEST_TIMEOUT}" \
      --no-tests=error \
      --output-on-failure \
      "${CTEST_JUNIT_OPTS[@]}"

# Post-run fallback: only if no XML was produced by CTest/Catch2
shopt -s nullglob
unit_xmls=("${REPORT_DIR}"/*.xml)
if [[ ${#unit_xmls[@]} -eq 0 ]]; then
  printf '<testsuite name="unit" tests="0" failures="0" skipped="0"/>\n' > "${REPORT_DIR}/ctest-unit.xml"
  echo "📄 JUnit (fallback): ${REPORT_DIR}/ctest-unit.xml"
fi
