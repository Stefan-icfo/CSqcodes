"""Repair a corrupted qcodes .db so plottr-inspectr can open it again.

Usage (qcodes env):  python repair_qcodes_db.py <corrupted.db> [<repaired.db>]

Pure-python rebuild: copies schema + all readable rows into a fresh db
(unreadable rows/pages are skipped, batch- then row-wise), restores the qcodes
schema version, prunes runs that still fail to load, then verifies with
plottr's own loader. The original file is never modified.
"""
import os, sqlite3, sys, time

BATCH = 50_000


def copy_table(new, t):
    try:
        new.execute(f'INSERT INTO main."{t}" SELECT * FROM old."{t}"')
        return
    except sqlite3.DatabaseError as e:
        bulk_err = e
    new.execute(f'DELETE FROM main."{t}"')
    try:
        maxrow = new.execute(f'SELECT max(rowid) FROM old."{t}"').fetchone()[0] or 0
    except sqlite3.DatabaseError:
        maxrow = 0
        try:
            for (r,) in new.execute(f'SELECT rowid FROM old."{t}"'):
                maxrow = r
        except sqlite3.DatabaseError:
            pass
        if maxrow:
            maxrow += BATCH
        else:
            print(f"  {t}: unreadable ({bulk_err}), table lost")
            return
    lost, row_err = 0, None
    for a in range(1, maxrow + 1, BATCH):
        b = min(a + BATCH - 1, maxrow)
        try:
            new.execute(f'INSERT INTO main."{t}" SELECT * FROM old."{t}"'
                        ' WHERE rowid BETWEEN ? AND ?', (a, b))
        except sqlite3.DatabaseError:
            for r in range(a, b + 1):
                try:
                    new.execute(f'INSERT INTO main."{t}" SELECT * FROM old."{t}"'
                                ' WHERE rowid=?', (r,))
                except sqlite3.DatabaseError as e:
                    lost += 1
                    row_err = row_err or e
    print(f"  {t}: {lost} rows lost (bulk: {bulk_err}; row: {row_err})")


def recover(src, dst):
    new = sqlite3.connect(f"file:{dst}", uri=True, isolation_level=None)
    new.execute("PRAGMA journal_mode=MEMORY")
    new.execute("PRAGMA synchronous=OFF")
    new.execute("ATTACH ? AS old", (f"file:{src}?mode=ro",))
    uv = new.execute("PRAGMA old.user_version").fetchone()[0]
    schema = new.execute("SELECT type, name, sql FROM old.sqlite_master"
                         " WHERE sql IS NOT NULL AND name NOT LIKE 'sqlite_%'").fetchall()
    tables = [n for t, n, s in schema if t == "table"]
    for t, n, s in schema:
        if t == "table":
            new.execute(s)
    print(f"copying {len(tables)} tables ...")
    t0 = time.time()
    for i, n in enumerate(tables, 1):
        copy_table(new, n)
        if i % 200 == 0 or i == len(tables):
            print(f"  {i}/{len(tables)} tables, {time.time() - t0:.0f} s")
    for t, n, s in schema:
        if t != "table":
            try:
                new.execute(s)
            except sqlite3.OperationalError as e:
                print(f"  {t} {n}: {e}")
    new.execute(f"PRAGMA user_version={uv}")
    new.execute("PRAGMA journal_mode=WAL")
    new.execute("INSERT INTO experiments SELECT DISTINCT exp_id, 'recovered', 'unknown', 0, 0,"
                " '{}-{}-{}', 0 FROM runs WHERE exp_id NOT IN (SELECT exp_id FROM experiments)")
    orph = [r[0] for r in new.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'results-%'"
        " AND name NOT IN (SELECT result_table_name FROM runs)")]
    if orph:
        print(f"  {len(orph)} orphaned result tables kept (runs-row lost),"
              f" e.g. {', '.join(orph[:5])}")
    new.close()
    print(f"rebuild done in {time.time() - t0:.0f} s")


def find_bad_runs(dst):
    from qcodes.dataset.data_set import DataSet
    from qcodes.dataset.sqlite.database import connect
    qcon = connect(dst)
    runs = qcon.execute("SELECT run_id, result_table_name FROM runs ORDER BY run_id").fetchall()
    bad = []
    for rid, tbl in runs:
        try:
            DataSet(conn=qcon, run_id=rid)
            qcon.execute(f'SELECT * FROM "{tbl}" LIMIT 1')
        except Exception as e:
            bad.append((rid, tbl, repr(e)[:150]))
    qcon.close()
    return len(runs), bad


def prune(dst, bad):
    con = sqlite3.connect(dst)
    for rid, tbl, err in bad:
        name = con.execute("SELECT name FROM runs WHERE run_id=?", (rid,)).fetchone()
        print(f"  pruning run {rid} {name}: {err}")
        con.execute(f'DROP TABLE IF EXISTS "{tbl}"')
        con.execute("DELETE FROM dependencies WHERE dependent IN"
                    " (SELECT layout_id FROM layouts WHERE run_id=?) OR independent IN"
                    " (SELECT layout_id FROM layouts WHERE run_id=?)", (rid, rid))
        con.execute("DELETE FROM layouts WHERE run_id=?", (rid,))
        con.execute("DELETE FROM runs WHERE run_id=?", (rid,))
    con.commit()
    con.close()


def main():
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(src)[0] + "_repaired.db"
    if os.path.exists(dst):
        sys.exit(f"output exists: {dst}")
    recover(src, dst)
    total, bad = find_bad_runs(dst)
    print(f"{total} runs recovered, {len(bad)} unloadable")
    if total == 0:
        sys.exit("no runs recovered -- repaired file is an empty shell, do not use")
    if bad:
        if len(bad) > total / 2:
            sys.exit("more than half the runs are bad -- not pruning, something is off")
        prune(dst, bad)
    from plottr.data.qcodes_dataset import get_runs_from_db_as_dataframe
    df = get_runs_from_db_as_dataframe(dst)
    print(f"plottr check OK: lists {len(df)} runs -> {dst}")


if __name__ == "__main__":
    main()
