#!/usr/bin/env bash
# Regenerate all non-golden fixtures (formatting + IO parity) from Python.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${WITCHI_CONDA_ENV:-witchi-0.2.0}"
cd "$ROOT"
python -I scripts/make_fixtures.py tests/fixtures
rm -rf tests/io
python -I scripts/make_io_cases.py tests/io/inputs
python -I scripts/make_io_fixtures.py tests/io/inputs tests/data --out tests/io/expected 2>&1 | grep -v Warning
mkdir -p tests/intermediates
while IFS=$'\t' read -r name input args; do
    [[ -z "$name" || "$name" == \#* ]] && continue
    # shellcheck disable=SC2086
    python -I scripts/dump_intermediates.py "tests/intermediates/$name.txt" --file "tests/data/$input" $args
done < tests/intermediates.tsv
