#!/bin/bash
# Local multi-compiler check: test in both profiles with every compiler found
# on PATH, run the examples and the app, then check llms.txt.
# GitHub Actions (.github/workflows/ci.yml) covers more platforms.

source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
# stdlib has C sources; fpm compiles them with icx for ifx, which may not be installed
command -v icx > /dev/null || export FPM_CC=gcc
FAILED=()

# No nvfortran or LLVM flang: they cannot build stdlib or fplotlib (see llms.txt).
for FC in gfortran ifx; do
    command -v "$FC" > /dev/null || continue
    for PROFILE in debug release; do
        echo "--- $FC $PROFILE"
        # The app writes its plot to build/, so the committed images stay untouched
        if fpm test --compiler "$FC" --profile "$PROFILE" \
            && fpm run --compiler "$FC" --profile "$PROFILE" --example --all > /dev/null \
            && fpm run --compiler "$FC" --profile "$PROFILE" -- -o build/modes.png > /dev/null; then
            echo "✅ $FC $PROFILE"
        else
            FAILED+=("$FC $PROFILE")
        fi
    done
done

python3 tools/check_docs.py || FAILED+=("llms.txt")

if [ ${#FAILED[@]} -eq 0 ]; then
    echo "✅ All checks passed"
else
    printf '❌ %s\n' "${FAILED[@]}"
    exit 1
fi
