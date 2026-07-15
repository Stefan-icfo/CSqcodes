"""Compare a repaired qcodes db against its corrupted source.

Usage: python verify_repair.py <corrupted.db> <repaired.db>
Prints meta-table counts and sampled result-table row counts (read-only).
"""
import sqlite3, sys

src, dst = sys.argv[1], sys.argv[2]
o = sqlite3.connect(f"file:{src}?mode=ro", uri=True)
n = sqlite3.connect(f"file:{dst}?mode=ro", uri=True)
print("user_version:", n.execute("PRAGMA user_version").fetchone()[0])
for t in ("runs", "experiments", "layouts", "dependencies"):
    print(t, o.execute(f"SELECT count(*) FROM {t}").fetchone()[0], "->",
          n.execute(f"SELECT count(*) FROM {t}").fetchone()[0])
last = n.execute("SELECT max(run_id) FROM runs").fetchone()[0]
if not last:
    sys.exit("repaired db has NO runs")
for rid in sorted({1, last // 3, 2 * last // 3, last - 1, last} - {0}):
    row = n.execute("SELECT result_table_name FROM runs WHERE run_id=?", (rid,)).fetchone()
    if not row:
        print(f"run {rid}: no runs-row in repaired (lost/pruned)")
        continue
    tbl = row[0]
    cn = n.execute(f'SELECT count(*) FROM "{tbl}"').fetchone()[0]
    try:
        co = o.execute(f'SELECT count(*) FROM "{tbl}"').fetchone()[0]
    except sqlite3.Error as e:
        co = f"unreadable in source ({e})"
    print(f"run {rid} ({tbl}): {co} -> {cn}", "OK" if co == cn else "")
