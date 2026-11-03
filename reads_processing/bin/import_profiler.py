#!/usr/bin/env python3
"""Opt-in, low-impact JSONL profiling of sample database imports.

No SQL mutations or transaction changes are made by this module. Timers use
perf_counter; SQLite diagnostics are read-only. Each record is one JSON object.
"""
import functools
import json
import logging
import os
import sqlite3
import time

_CONTEXT = None
_PLANS_SEEN = set()


def configure(output, sample, tool):
    global _CONTEXT
    _CONTEXT = {"output": output, "sample": sample, "tool": tool}
    _PLANS_SEEN.clear()
    emit("session", status="started", python_sqlite=sqlite3.sqlite_version)


def emit(stage, **data):
    if _CONTEXT is None:
        return
    record = {"sample": _CONTEXT["sample"], "tool": _CONTEXT["tool"],
              "stage": stage, **data}
    try:
        with open(_CONTEXT["output"], "a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, default=str) + "\n")
    except OSError:
        logging.exception("Failed to write profiling data")
    logging.info("IMPORT_PROFILE %s", json.dumps(record, sort_keys=True, default=str))


def _rows(obj):
    if obj is None:
        return 0
    if hasattr(obj, "shape"):
        return int(obj.shape[0])
    if isinstance(obj, tuple):
        return [_rows(item) for item in obj]
    if isinstance(obj, (int, float)):
        return obj
    return None


def profile_functions(namespace, names):
    """Wrap existing functions transparently, retaining their return values."""
    if _CONTEXT is None:
        return
    for name in names:
        original = namespace[name]
        if getattr(original, "_import_profile_wrapped", False):
            continue

        @functools.wraps(original)
        def wrapped(*args, __original=original, __name=name, **kwargs):
            start = time.perf_counter()
            input_rows = _rows(args[0]) if args else None
            status = "success"
            try:
                result = __original(*args, **kwargs)
                return result
            except BaseException:
                status = "failed"
                raise
            finally:
                emit(__name, elapsed_s=round(time.perf_counter()-start, 6),
                     input_rows=input_rows,
                     output_rows=_rows(locals().get("result")),
                     status=status)
        wrapped._import_profile_wrapped = True
        namespace[name] = wrapped


def database_snapshot(conn, db_path, stage):
    if _CONTEXT is None:
        return
    start = time.perf_counter()
    info = {}
    for pragma in ("page_count", "page_size", "freelist_count", "journal_mode",
                   "synchronous", "cache_size", "temp_store"):
        info[pragma] = conn.execute("PRAGMA " + pragma).fetchone()[0]
    try:
        info["file_bytes"] = os.path.getsize(db_path)
    except OSError:
        info["file_bytes"] = None
    info["indexes"] = {}
    for table in ("splicing_events", "junction_sample_stats", "junction_summary_stats"):
        info["indexes"][table] = [row[1] for row in
            conn.execute("PRAGMA index_list(" + table + ")").fetchall()]
    # sqlite_stat1 has approximate cardinalities if ANALYZE has been run.
    has_stats = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE name='sqlite_stat1'"
    ).fetchone()
    if has_stats:
        info["sqlite_stat1"] = [
            {"table": tab, "index": index, "stat": stat}
            for tab, index, stat in conn.execute(
                "SELECT tbl, idx, stat FROM sqlite_stat1"
            ).fetchall()
        ]
    emit(stage, elapsed_s=round(time.perf_counter()-start, 6), **info)


def explain_once(conn, name, sql, bindings=()):
    """Capture EXPLAIN QUERY PLAN once per process/query type, read-only."""
    if _CONTEXT is None or name in _PLANS_SEEN:
        return
    _PLANS_SEEN.add(name)
    start = time.perf_counter()
    try:
        plan = conn.execute("EXPLAIN QUERY PLAN " + sql, bindings).fetchall()
        emit("query_plan", name=name, plan=[row[3] for row in plan],
             elapsed_s=round(time.perf_counter()-start, 6))
    except sqlite3.Error as exc:
        emit("query_plan", name=name, error=str(exc))


def finish():
    emit("session", status="finished")
