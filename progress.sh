#!/usr/bin/env bash
# Pass/fail overview of every exercise. Usage: ./progress.sh
cd "$(dirname "$0")" || exit 1
done_count=0; total=0
for dir in unit*/*/; do
  crate=$(basename "$dir")
  total=$((total + 1))
  out=$(cargo test -q -p "$crate" 2>&1)
  if echo "$out" | grep -q "could not compile"; then
    printf "[!] %-26s  does not compile yet\n" "$crate"
  else
    line=$(echo "$out" | grep -m1 "test result")
    passed=$(echo "$line" | grep -oE '[0-9]+ passed' | grep -oE '[0-9]+')
    failed=$(echo "$line" | grep -oE '[0-9]+ failed' | grep -oE '[0-9]+')
    if [ "${failed:-1}" = "0" ]; then
      printf "[x] %-26s  %s/%s\n" "$crate" "$passed" "$passed"
      done_count=$((done_count + 1))
    else
      printf "[ ] %-26s  %s/%s tests passing\n" "$crate" "${passed:-0}" "$(( ${passed:-0} + ${failed:-0} ))"
    fi
  fi
done
echo
echo "$done_count / $total exercises complete"
