#!/usr/bin/env bash
# Shared helpers for the dev/refdb build scripts. Sourced, not executed.
#
# See dev/design/18_REFDB_REBUILD.md for the spec this implements.
set -Eeuo pipefail
export LC_ALL=C

# Defense in depth: a script that sources this and forgets its own trap still
# gets a message instead of a silent bare exit.
trap 'echo "lib.sh: unexpected error at line $LINENO" >&2' ERR

# rp_log MESSAGE...
# Timestamped status line, scrubbed on stderr and (scrubbed) to RP_LOG_FILE if
# set. Scrubbed here as well as by rp_scrub_output so a script that sources
# this and never installs the stream filter still cannot log a key.
rp_log() {
  local msg
  msg="[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
  rp_scrub "$msg" >&2
  if [ -n "${RP_LOG_FILE:-}" ]; then
    rp_scrub "$msg" >> "$RP_LOG_FILE"
  fi
}

# rp_die MESSAGE...
# Log and exit 1. Every abort path in the orchestrator funnels through here so
# every abort reason is scrubbed and timestamped the same way.
rp_die() {
  rp_log "ERROR: $*"
  exit 1
}

# rp_need TOOL
# Abort naming the missing tool, since a bare "command not found" later, deep
# inside a stage, does not say which of a dozen dependencies is absent.
rp_need() {
  local tool="$1"
  command -v "$tool" >/dev/null 2>&1 \
    || rp_die "required tool not found on PATH: $tool (install it, e.g. via conda/bioconda, and retry)"
}

# rp_scrub [TEXT...]
# Redact the API key from stdin, or from TEXT if given. api_key= and
# NCBI_API_KEY= are the same secret in two spellings, so both are matched,
# either case, with or without the quote an env dump adds. The historical
# harvest.log had 35 plaintext API keys in it.
RP_SCRUB_RE='s/([Aa][Pp][Ii][_-]?[Kk][Ee][Yy][=:] *"?)[A-Za-z0-9]+/\1REDACTED/g'
rp_scrub() {
  if [ "$#" -gt 0 ]; then
    printf '%s\n' "$*" | sed -E "$RP_SCRUB_RE"
  else
    sed -E "$RP_SCRUB_RE"
  fi
}

# rp_scrub_output
# Route this shell's own stdout and stderr through rp_scrub, so no command can
# leak the key to a terminal or to `2>&1 | tee build.log`, including ones this
# script does not print itself (set -x traces, env dumps, EDirect URLs).
# Idempotent.
rp_scrub_output() {
  [ -z "${RP_SCRUB_PIDS:-}" ] || return 0
  local out err
  exec > >(rp_scrub)
  out="$!"
  exec 2> >(rp_scrub >&2)
  err="$!"
  RP_SCRUB_PIDS="$out $err"
}

# rp_scrub_flush
# Close the filtered streams and wait for the filters, so the last lines are
# not still in flight when the script exits and a reader sees a truncated log.
rp_scrub_flush() {
  [ -n "${RP_SCRUB_PIDS:-}" ] || return 0
  local pids="$RP_SCRUB_PIDS"
  RP_SCRUB_PIDS=""
  exec >&- 2>&-
  # shellcheck disable=SC2086
  wait $pids 2>/dev/null || true
}

# rp_lock LOCKFILE
# PID lockfile; refuses a second concurrent run. Never pgrep, never match on
# a command line: the original chain waiter deadlocked because `pgrep -f`
# matched its own argv. Caller is responsible for removing LOCKFILE on exit.
rp_lock() {
  local lockfile="$1" held_pid claim attempt
  for attempt in 1 2 3; do
    # noclobber makes "create only if absent" one atomic step. Testing for the
    # file and then writing it is two, and two concurrent runs both pass the
    # test and both write, which is the interleaved build tree the lock exists
    # to prevent.
    if (set -C; printf '%s\n' "$$" > "$lockfile") 2>/dev/null; then
      return 0
    fi
    held_pid="$(cat "$lockfile" 2>/dev/null || true)"
    if [ -n "$held_pid" ] && kill -0 "$held_pid" 2>/dev/null; then
      rp_die "another build is already running: pid $held_pid holds $lockfile"
    fi
    # Stale (recorded pid gone) or empty (holder died mid-write): reclaim it,
    # rather than deadlock forever on a lock nobody holds. Renaming it away
    # first is atomic too, so a racer reclaiming at the same instant can see it
    # took a live lock, put it back, and refuse instead of both proceeding.
    claim="$lockfile.stale.$$"
    if mv "$lockfile" "$claim" 2>/dev/null; then
      held_pid="$(cat "$claim" 2>/dev/null || true)"
      if [ -n "$held_pid" ] && kill -0 "$held_pid" 2>/dev/null; then
        mv "$claim" "$lockfile"
        rp_die "another build is already running: pid $held_pid holds $lockfile"
      fi
      rm -f "$claim"
      rp_log "stale lock $lockfile (pid ${held_pid:-unknown} not running); removed"
    fi
  done
  rp_die "could not acquire $lockfile after 3 attempts"
}

# rp_marker_ok STAGE INPUT_HASH
# True only when RP_STATE_FILE has a line for STAGE recording the same hash.
# RP_STATE_FILE is one small text file (stage, hash, timestamp per line), so
# `status` never rescans the build tree.
rp_marker_ok() {
  local stage="$1" hash="$2" line tab
  tab="$(printf '\t')"
  [ -n "${RP_STATE_FILE:-}" ] && [ -f "$RP_STATE_FILE" ] || return 1
  line="$(grep "^${stage}${tab}" "$RP_STATE_FILE" || true)"
  [ -n "$line" ] || return 1
  [ "$(printf '%s' "$line" | cut -f2)" = "$hash" ]
}

# rp_marker_present STAGE
# True when RP_STATE_FILE records STAGE at all, whatever hash it holds. Distinct
# from rp_marker_ok, which conflates "no marker" with "marker disagrees": a stage
# that must refuse on a DISAGREEING marker still has to run on a tree that has
# none, so the two cases cannot share an answer.
rp_marker_present() {
  local stage="$1" tab
  tab="$(printf '\t')"
  [ -n "${RP_STATE_FILE:-}" ] && [ -f "$RP_STATE_FILE" ] || return 1
  grep -q "^${stage}${tab}" "$RP_STATE_FILE"
}

# rp_marker_set STAGE INPUT_HASH
# Record STAGE as complete for INPUT_HASH, replacing any prior line for STAGE.
rp_marker_set() {
  local stage="$1" hash="$2" tmp tab
  tab="$(printf '\t')"
  [ -n "${RP_STATE_FILE:-}" ] || rp_die "rp_marker_set: RP_STATE_FILE is not set"
  mkdir -p "$(dirname "$RP_STATE_FILE")"
  tmp="$(mktemp "${RP_STATE_FILE}.XXXXXX")"
  if [ -f "$RP_STATE_FILE" ]; then
    grep -v "^${stage}${tab}" "$RP_STATE_FILE" > "$tmp" || true
  else
    : > "$tmp"
  fi
  printf '%s\t%s\t%s\n' "$stage" "$hash" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$tmp"
  mv "$tmp" "$RP_STATE_FILE"
}

# rp_marker_clear STAGE
# Drop only STAGE's line. Deleting the whole state file instead would re-run
# stages whose outputs are still on disk, including the 2 h 20 m fetch.
rp_marker_clear() {
  local stage="$1" tmp tab
  tab="$(printf '\t')"
  [ -n "${RP_STATE_FILE:-}" ] && [ -f "$RP_STATE_FILE" ] || return 0
  tmp="$(mktemp "${RP_STATE_FILE}.XXXXXX")"
  grep -v "^${stage}${tab}" "$RP_STATE_FILE" > "$tmp" || true
  mv "$tmp" "$RP_STATE_FILE"
}
