import sqlite3
import re

# matches: "250ndemod"  OR  "1.45udemod"  (spaces optional)
INT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)\s*([nu])\s*demod", re.IGNORECASE)

def _table_cols(cur, table):
    cur.execute(f"PRAGMA table_info({table});")
    return [r[1] for r in cur.fetchall()]

def _pick_first(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None

def _parse_t_ns(label: str):
    if not label:
        return None
    m = INT_RE.search(label)
    if not m:
        return None
    val = float(m.group(1))
    unit = m.group(2).lower()
    return val if unit == "n" else val * 1000.0  # u -> µs -> 1000 ns

def build_integration_groups_from_db(db_path, run_first, run_last, tol_ns=1e-6, verbose=True):
    con = sqlite3.connect(db_path)
    cur = con.cursor()

    runs_cols = _table_cols(cur, "runs")
    exps_cols = _table_cols(cur, "experiments") if "experiments" in [r[0] for r in cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table';"
    ).fetchall()] else []

    # choose which columns likely contain name text
    run_name_col = _pick_first(runs_cols, ["name", "result_name", "run_name", "captured_run_name"])
    exp_id_col   = _pick_first(runs_cols, ["experiment_id", "exp_id", "experiment", "exp"])
    exp_name_col = _pick_first(exps_cols, ["name", "exp_name", "experiment_name"]) if exps_cols else None

    if verbose:
        print("runs columns:", runs_cols)
        if exps_cols:
            print("experiments columns:", exps_cols)
        print("Chosen run_name_col =", run_name_col)
        print("Chosen exp_id_col   =", exp_id_col)
        print("Chosen exp_name_col =", exp_name_col)

    # We build a query that returns some text label for each run.
    # Priority: runs.<name-like> , else experiments.<name-like> (if join possible)
    rows = []

    if run_name_col is not None:
        cur.execute(
            f"SELECT run_id, {run_name_col} FROM runs WHERE run_id BETWEEN ? AND ? ORDER BY run_id ASC;",
            (run_first, run_last)
        )
        rows = cur.fetchall()

    # If the run label is empty or missing, try joining to experiments (if possible)
    if (not rows or all((r[1] is None or str(r[1]).strip() == "") for r in rows)) and exp_id_col and exp_name_col:
        cur.execute(
            f"""
            SELECT r.run_id, e.{exp_name_col}
            FROM runs r
            LEFT JOIN experiments e ON r.{exp_id_col} = e.exp_id
            WHERE r.run_id BETWEEN ? AND ?
            ORDER BY r.run_id ASC;
            """,
            (run_first, run_last)
        )
        rows = cur.fetchall()

    con.close()

    # parse + group
    groups = []
    cur_t = None
    cur_start = None
    cur_end = None

    for rid, label in rows:
        label = "" if label is None else str(label)
        t_ns = _parse_t_ns(label)

        if verbose and (rid in (run_first, run_first+1, run_first+2, run_last-2, run_last-1, run_last)):
            print(f"[debug run {rid}] label='{label[:120]}'  ->  t_ns={t_ns}")

        if t_ns is None:
            if cur_t is not None:
                groups.append({
                    "label": (f"{cur_t/1000:.3f} us" if cur_t >= 1000 else f"{cur_t:.1f} ns"),
                    "t_ns": float(cur_t),
                    "run_first": int(cur_start),
                    "run_last": int(cur_end),
                })
                cur_t = cur_start = cur_end = None
            continue

        if cur_t is None:
            cur_t = t_ns
            cur_start = rid
            cur_end = rid
        else:
            same_t = abs(t_ns - cur_t) <= tol_ns
            consecutive = (rid == cur_end + 1)
            if same_t and consecutive:
                cur_end = rid
            else:
                groups.append({
                    "label": (f"{cur_t/1000:.3f} us" if cur_t >= 1000 else f"{cur_t:.1f} ns"),
                    "t_ns": float(cur_t),
                    "run_first": int(cur_start),
                    "run_last": int(cur_end),
                })
                cur_t = t_ns
                cur_start = rid
                cur_end = rid

    if cur_t is not None:
        groups.append({
            "label": (f"{cur_t/1000:.3f} us" if cur_t >= 1000 else f"{cur_t:.1f} ns"),
            "t_ns": float(cur_t),
            "run_first": int(cur_start),
            "run_last": int(cur_end),
        })

    return groups
