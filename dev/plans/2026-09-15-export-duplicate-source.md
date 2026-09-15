# Export: duplicate source info across contigs of one sample

Status: parked. Ideas only, nothing built.

## Problem

NCBI rejects a submission when two sequences share the same source
modifier combination. A sample with several exported contigs (path 1,
scaffolds 1 and 2) ships the same Organism and specimen-voucher on each:

    Sequence_ID             Organism               specimen-voucher
    22030FL-06-01-12_p1_s1  Monomitopus agassizii  USNM:FISH:407074
    22030FL-06-01-12_p1_s2  Monomitopus agassizii  USNM:FISH:407074

## Option 1: distinguishing modifier in the header template (preferred)

Cheapest. The FASTA header is a glue template (export_opts.R) and the header
data already carries `seqid`, `path`, and `scaffold` per unit
(export_header_data). Add a modifier that differs per contig, e.g.

    [isolate={seqid}]        or        [clone={ID}_p{path}_s{scaffold}]

Choices:
- Change DEFAULT_FASTA_HEADER so new projects get it. Pre-release, so no
  migration for old projects (see memory: no schema migrations yet).
- Or only add it when a sample exports more than one unit. Needs a per-unit
  conditional in the template renderer; more code, cleaner headers for the
  common one-contig case.
- Which qualifier: `isolate` definitely counts toward NCBI uniqueness;
  `note` may not. `clone` is semantically closer for multiple copies from one
  voucher. Confirm with a test submission before settling.

Also worth an export check (header_status_ui style): warn when two units in
one export share identical modifier sets, so the user sees it before NCBI does.

## Option 2: auto-split into export groups

Assign contigs of the same sample to different export groups so each
submission has one contig per sample. Bigger: group assignment logic, N
submissions per N-contig sample, the user has to track them all. Only worth it
if NCBI refuses option 1.

## Open question

Does NCBI accept the same voucher with different isolate/clone values for
rRNA contigs? If not, option 2 is the fallback.
