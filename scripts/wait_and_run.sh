#!/bin/bash
# Wait for the standalone TC_1 run (launched outside run_all.sh, orphaned from
# an earlier driver bug) to fully finish, then run the rest strictly one at a
# time. Two concurrent long-running streams against the SMB mount have twice
# triggered "gzcat: couldn't close input: Bad file descriptor" -- so no
# overlap at all this time.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOGFILE="$PROJECT_DIR/results/TC_1.run.log"

while true; do
  if grep -q "\[TC_1\] done\." "$LOGFILE" 2>/dev/null; then
    echo "TC_1 finished cleanly." >&2
    break
  fi
  if ! pgrep -f "run_sample.sh TC_1" > /dev/null 2>&1; then
    echo "WARNING: TC_1 process gone but 'done' marker not found in log -- may have crashed." >&2
    break
  fi
  sleep 30
done

echo "Starting TC_2..TC_6 sequentially." >&2
bash "$PROJECT_DIR/scripts/run_all.sh" TC_2 TC_3 TC_4 TC_5 TC_6
