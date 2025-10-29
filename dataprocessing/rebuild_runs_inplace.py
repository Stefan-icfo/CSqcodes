# rebuild_runs_inplace.py
import sqlite3

DB = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4V21_27_10_25.db"

def main():
    conn = sqlite3.connect(DB)
    cur = conn.cursor()

    print("Opening DB…")
    cur.execute("PRAGMA foreign_keys=OFF;")
    conn.commit()

    # --- 1) Read schema of runs (DDL) ---
    ddl = cur.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='runs'"
    ).fetchone()
    if not ddl or not ddl[0]:
        raise RuntimeError("Could not read DDL for 'runs' table.")
    runs_ddl = ddl[0]

    # --- 2) Get column list for runs ---
    cols_info = cur.execute("PRAGMA table_info(runs)").fetchall()
    if not cols_info:
        raise RuntimeError("PRAGMA table_info(runs) returned nothing.")
    cols = [r[1] for r in cols_info]
    cols_csv = ", ".join(f'"{c}"' for c in cols)
    placeholders = ", ".join("?" for _ in cols)

    # --- 3) Create empty staging table with identical schema ---
    # Trick: build a CREATE TABLE statement for runs_ok by tweaking name
    runs_ok_ddl = runs_ddl.replace('CREATE TABLE "runs"', 'CREATE TABLE "runs_ok"')\
                          .replace('CREATE TABLE runs', 'CREATE TABLE runs_ok')
    # Fall back to SELECT .. WHERE 0; if CREATE fails (rare)
    try:
        cur.execute("DROP TABLE IF EXISTS runs_ok;")
        cur.execute(runs_ok_ddl)
        conn.commit()
    except sqlite3.DatabaseError:
        cur.execute("DROP TABLE IF EXISTS runs_ok;")
        cur.execute(f'CREATE TABLE runs_ok AS SELECT {cols_csv} FROM runs WHERE 0;')
        conn.commit()

    print("Copying rows from runs → runs_ok (skipping corrupt pages)…")
    # --- 4) Chunked copy that narrows window on errors ---
    # Use LIMIT/OFFSET. On error, shrink step; if step==1 and error persists, skip that row.
    offset, step = 0, 20000
    copied = 0
    while True:
        try:
            rows = cur.execute(
                f'SELECT {cols_csv} FROM runs LIMIT ? OFFSET ?;', (step, offset)
            ).fetchall()
        except sqlite3.DatabaseError:
            # Narrow window to hop over the bad page(s)
            if step == 1:
                offset += 1
            else:
                step = max(1, step // 2)
            continue

        if not rows:
            break

        cur.executemany(
            f'INSERT OR IGNORE INTO runs_ok ({cols_csv}) VALUES ({placeholders})',
            rows
        )
        conn.commit()
        copied += len(rows)
        offset += len(rows)
        # Slowly increase step again for speed
        if step < 20000:
            step = min(20000, step * 2)

        if copied % 50000 == 0:
            print(f"  Copied {copied} rows…")

    print(f"Finished copy. Total rows copied: {copied}")

    # --- 5) Swap tables (keep the original as runs_bad for now) ---
    print("Swapping tables…")
    cur.execute('ALTER TABLE runs RENAME TO runs_bad;')
    cur.execute('ALTER TABLE runs_ok RENAME TO runs;')
    conn.commit()

    # --- 6) Drop broken indexes that reference runs (if any in sqlite_master) ---
    # We don't touch data; just remove references that point to corrupt pages.
    idx_names = [r[0] for r in cur.execute(
        "SELECT name FROM sqlite_master WHERE type='index' AND tbl_name='runs'"
    ).fetchall()]
    for name in idx_names:
        try:
            cur.execute(f'DROP INDEX IF EXISTS "{name}";')
        except sqlite3.DatabaseError:
            # if the index btree itself is corrupted, DROP may fail harmlessly
            pass
    conn.commit()

    # --- 7) Minimal sanity check ---
    try:
        max_id = cur.execute('SELECT max(run_id) FROM runs').fetchone()
        print("max(run_id) in rebuilt runs:", max_id)
    except sqlite3.DatabaseError as e:
        print("Warning: SELECT max(run_id) failed even after rebuild:", e)

    # Optional: quick check. Full integrity_check may still report other tables;
    # we only care that runs is readable again for QCoDeS to start.
    try:
        chk = cur.execute("PRAGMA quick_check;").fetchone()
        print("PRAGMA quick_check:", chk)
    except sqlite3.DatabaseError as e:
        print("PRAGMA quick_check failed:", e)

    conn.close()
    print("Done. Try QCoDeS again (load_by_id). If it works, you can later drop 'runs_bad' to reclaim space.")

if __name__ == "__main__":
    main()
