#!/bin/bash
FASTA="${1:-$HOME/enzyme/fasta}"
OUT="${2:-$HOME/enzyme/folded}"
LOG="${3:-$HOME/enzyme/boltz.log}"
WIDTH=42
START=$(date +%s)

TOTAL=$(ls "$FASTA"/*.fasta 2>/dev/null | wc -l)
[ "$TOTAL" -eq 0 ] && { echo "no fasta files in $FASTA"; exit 1; }

while true; do
  done_n=$(find "$OUT" -name "*_model_0.pdb" 2>/dev/null | wc -l)
  [ "$done_n" -gt "$TOTAL" ] && done_n=$TOTAL

  pct=$(( done_n * 100 / TOTAL ))
  filled=$(( done_n * WIDTH / TOTAL ))
  bar=$(printf '%*s' "$filled" '' | tr ' ' '=')
  gap=$(printf '%*s' $(( WIDTH - filled )) '' | tr ' ' ' ')

  now=$(date +%s); elapsed=$(( now - START ))
  if [ "$done_n" -gt 0 ] && [ "$elapsed" -gt 0 ]; then
    per=$(( elapsed / done_n ))
    eta=$(( per * (TOTAL - done_n) ))
    eh=$(( eta / 3600 )); em=$(( (eta % 3600) / 60 ))
    etastr=$(printf '%dh%02dm' "$eh" "$em")
    perstr=$(printf '%ds/ea' "$per")
  else
    etastr="--"; perstr="--"
  fi

  alive=$(pgrep -f "boltz predict" >/dev/null && echo RUN || echo STOPPED)

  printf '\r[%s%s] %3d%%  %2d/%d  %s  ETA %s  %s ' \
    "$bar" "$gap" "$pct" "$done_n" "$TOTAL" "$perstr" "$etastr" "$alive"

  if [ "$done_n" -ge "$TOTAL" ]; then printf '\nAll %d predictions complete.\n' "$TOTAL"; break; fi
  if [ "$alive" = "STOPPED" ]; then
    printf '\nProcess is not running. Last lines of the log:\n'
    tail -15 "$LOG" 2>/dev/null
    break
  fi
  sleep 15
done
