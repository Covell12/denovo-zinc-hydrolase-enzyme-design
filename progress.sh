#!/bin/bash
LOG="${1:-$HOME/enzyme/production.log}"
TOTAL="${2:-32}"
WIDTH=42

while true; do
  done_n=$(grep -c "Finished design" "$LOG" 2>/dev/null)
  [ -z "$done_n" ] && done_n=0
  step=$(grep "Flow-matching step" "$LOG" 2>/dev/null | tail -1 | sed -n 's|.*step: \([0-9]*\)/100.*|\1|p')
  [ -z "$step" ] && step=0

  units=$(( done_n * 100 + step ))
  total=$(( TOTAL * 100 ))
  [ "$units" -gt "$total" ] && units=$total
  pct=$(( units * 100 / total ))
  filled=$(( units * WIDTH / total ))

  bar=$(printf '%*s' "$filled" '' | tr ' ' '=')
  gap=$(printf '%*s' $(( WIDTH - filled )) '' | tr ' ' ' ')

  cur=$(( done_n + 1 ))
  [ "$cur" -gt "$TOTAL" ] && cur=$TOTAL

  left=$(( total - units ))
  eta_min=$(( left * 11 / 60 ))
  eta_h=$(( eta_min / 60 )); eta_m=$(( eta_min % 60 ))

  printf '\r[%s%s] %3d%%  design %2d/%d  step %3d/100  ETA %dh%02dm ' \
    "$bar" "$gap" "$pct" "$cur" "$TOTAL" "$step" "$eta_h" "$eta_m"

  if [ "$done_n" -ge "$TOTAL" ]; then
    printf '\nAll %d designs complete.\n' "$TOTAL"
    break
  fi
  sleep 10
done
