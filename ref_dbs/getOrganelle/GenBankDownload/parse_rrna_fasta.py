#!/usr/bin/env python3
"""
Parse downloaded nuclear ribosomal RNA (nrDNA) GenBank FASTA records into
GetOrganelle seed + label databases.

Records are classified from their description line (GenBank title), since most
nuclear rRNA submissions are feature-less single- or multi-element sequences:

  "... 18S ribosomal RNA gene, partial sequence"                -> 18S
  "... 28S ribosomal RNA gene ..."                              -> 28S
  "... 5.8S ribosomal RNA gene ..."                             -> 5.8S
  "... 18S ... internal transcribed spacer 1, 5.8S ... 28S ..." -> multi (cistron)

Outputs:
  seeds  : read-recruiting sequences (18S / 28S / multi-element cistron),
           headers ">ACCESSION_Genus_species" (mirrors MitoPilot seed format)
  labels : per-gene sequences for contig classification, headers ">18S" etc.

Adapted from GenBankMitoReview/filter_gb.py.
"""

import argparse
import hashlib
import re
from Bio import SeqIO

# Per-element minimum lengths (bp) to keep noise out.
SEED_MIN = {"18S": 1000, "28S": 600, "multi": 800}
LABEL_MIN = {"18S": 1000, "28S": 600, "5.8S": 100}
LABEL_MAX = {"5.8S": 400}  # 5.8S is ~160 bp; reject anything oversized


def classify(desc):
    """Return (set_of_rRNA_elements, has_ITS) from a description line."""
    d = desc.lower()
    elements = set()
    if re.search(r"\b18s\b", d) or "small subunit ribosomal" in d:
        elements.add("18S")
    if re.search(r"\b(28s|26s|25s)\b", d) or "large subunit ribosomal" in d:
        elements.add("28S")
    if re.search(r"\b5\.8s\b", d):
        elements.add("5.8S")
    has_its = ("internal transcribed spacer" in d) or bool(re.search(r"\bits[12]?\b", d))
    return elements, has_its


def genus_species(desc):
    """First two alphabetic tokens after the accession -> 'Genus_species'."""
    toks = desc.split()
    name = [t for t in toks[1:] if re.match(r"^[A-Za-z][A-Za-z.-]*$", t)][:2]
    if not name:
        return "unknown"
    return "_".join(re.sub(r"[^A-Za-z0-9]", "", t) for t in name)


def accession(rec_id):
    return rec_id.split()[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="input FASTA (from download_rrna.sh)")
    ap.add_argument("--seeds", required=True, help="output seeds FASTA")
    ap.add_argument("--labels", required=True, help="output labels FASTA")
    ap.add_argument("--max-seeds", type=int, default=3000,
                    help="cap total seeds, keeping longest (default 3000)")
    ap.add_argument("--max-labels-per-gene", type=int, default=3000,
                    help="cap labels per gene, keeping longest (default 3000)")
    args = ap.parse_args()

    seed_recs = []            # list of (length, header, seq)
    label_recs = {"18S": [], "28S": [], "5.8S": []}
    seen_seed = set()
    seen_label = {"18S": set(), "28S": set(), "5.8S": set()}

    n_total = n_skip = 0
    for rec in SeqIO.parse(args.input, "fasta"):
        n_total += 1
        seq = str(rec.seq).upper().replace("\n", "")
        if not seq:
            continue
        elements, has_its = classify(rec.description)
        if not elements and not has_its:
            n_skip += 1
            continue

        acc = accession(rec.id)
        sp = genus_species(rec.description)
        L = len(seq)
        h = hashlib.md5(seq.encode()).hexdigest()

        multi = has_its or len(elements) > 1

        if multi:
            # whole cistron / multi-element -> seed only
            if L >= SEED_MIN["multi"] and h not in seen_seed:
                seen_seed.add(h)
                seed_recs.append((L, f"{acc}_{sp}", seq))
        else:
            gene = next(iter(elements))
            # single-gene -> label (all three) + seed (18S/28S)
            lo = LABEL_MIN.get(gene, 0)
            hi = LABEL_MAX.get(gene, 10**9)
            if lo <= L <= hi and h not in seen_label[gene]:
                seen_label[gene].add(h)
                label_recs[gene].append((L, gene, seq))
            if gene in ("18S", "28S") and L >= SEED_MIN[gene] and h not in seen_seed:
                seen_seed.add(h)
                seed_recs.append((L, f"{acc}_{sp}", seq))

    # cap, preferring longer (more informative) sequences
    seed_recs.sort(key=lambda x: -x[0])
    seed_recs = seed_recs[: args.max_seeds]
    for g in label_recs:
        label_recs[g].sort(key=lambda x: -x[0])
        label_recs[g] = label_recs[g][: args.max_labels_per_gene]

    with open(args.seeds, "w") as fh:
        for _, header, seq in seed_recs:
            fh.write(f">{header}\n{seq}\n")
    with open(args.labels, "w") as fh:
        for g in ("18S", "5.8S", "28S"):
            for _, header, seq in label_recs[g]:
                fh.write(f">{header}\n{seq}\n")

    n_lab = {g: len(v) for g, v in label_recs.items()}
    print(f"Parsed {n_total} record(s); skipped {n_skip} non-rRNA.")
    print(f"Seeds written: {len(seed_recs)}")
    print(f"Labels written: 18S={n_lab['18S']}, 5.8S={n_lab['5.8S']}, 28S={n_lab['28S']}")


if __name__ == "__main__":
    main()
