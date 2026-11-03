#!/bin/bash

# SPDX-License-Identifier: GPL-3

ROOT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
NPROC=$(nproc)
TMPDIR=$(mktemp -d)

echo "Running clang-tidy on $NPROC threads..."

find "$ROOT_DIR/src/oiseau" -name '*.cpp' -print0 | xargs -0 -P"$NPROC" -I{} bash -c '
  BASENAME=$(basename "{}" .cpp)
  CLANG_FORCE_COLOR=1 clang-tidy -p "'"$ROOT_DIR"'/build" "{}" --header-filter="^'"$ROOT_DIR"'/src/oiseau/" --use-color > "'"$TMPDIR"'/${BASENAME}.log" 2>/dev/null
  echo -n "."
'

echo -e "\nFinished.\n"

error_count=0
for LOG in "$TMPDIR"/*.log; do
  if [[ -s "$LOG" ]]; then
    cat "$LOG"
    error_count=$(( error_count + $(grep -c "error:" "$LOG") ))
  fi
done

if (( error_count > 0 )); then
  echo -e "\033[31m❌ Errors found: $error_count\033[0m"
else
  echo -e "\033[32m✅ All checks passed!\033[0m"
fi

rm -rf "$TMPDIR"

(( error_count > 255 )) && error_count=255
exit $error_count
