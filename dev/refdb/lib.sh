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
# Timestamped status line to stderr, and (scrubbed) to RP_LOG_FILE if set.
rp_log() {
  local msg
  msg="[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
  echo "$msg" >&2
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
# Redact api_key=... from stdin, or from TEXT if given. Applied to everything
# written to a log; the historical harvest.log had 35 plaintext API keys in it.
rp_scrub() {
  if [ "$#" -gt 0 ]; then
    printf '%s\n' "$*" | sed -E 's/api_key=[A-Za-z0-9]*/api_key=REDACTED/g'
  else
    sed -E 's/api_key=[A-Za-z0-9]*/api_key=REDACTED/g'
  fi
}

# rp_lock LOCKFILE
# PID lockfile; refuses a second concurrent run. Never pgrep, never match on
# a command line: the original chain waiter deadlocked because `pgrep -f`
# matched its own argv. Caller is responsible for removing LOCKFILE on exit.
rp_lock() {
  local lockfile="$1" held_pid
  if [ -e "$lockfile" ]; then
    held_pid="$(cat "$lockfile" 2>/dev/null || true)"
    if [ -n "$held_pid" ] && kill -0 "$held_pid" 2>/dev/null; then
      rp_die "another build is already running: pid $held_pid holds $lockfile"
    fi
    rp_log "stale lock $lockfile (pid ${held_pid:-unknown} not running); removing"
  fi
  printf '%s\n' "$$" > "$lockfile"
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
