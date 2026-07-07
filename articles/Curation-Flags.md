# Understanding the automatic quality checks

After it annotates a sample, RiboPilot runs a set of automatic quality
checks and raises a **flag** for anything that looks unusual. Flags are
there to guide your attention: they point out samples worth a closer
look in the app’s Annotate tab. They do not change your data, and a flag
does not always mean something is wrong. It just means “please check
this”.

You will see flags in two places in the app: a small count next to each
sample in the table, and a detailed **Curation flags** panel when you
open a sample. Each flag links to the exact feature it concerns, so you
can jump straight to it.

## What the flags mean

| Flag | Applies to | What it means |
|----|----|----|
| `missing_18S`, `missing_5.8S`, `missing_28S` | sample | One of the three ribosomal genes was not found. |
| `missing_ITS1`, `missing_ITS2` | sample | A spacer was not recovered between its neighboring genes. |
| `out_of_order` | sample | The five features are not in the expected order (18S, ITS1, 5.8S, ITS2, 28S). |
| `overlap` | feature | Two features overlap when they should sit side by side. |
| `wrong_strand` | feature | A feature was found on the opposite strand from what is expected. |
| `length` | feature | A feature is shorter or longer than typical for your group. |
| `low_identity` | feature | A conserved gene is less similar to reference sequences than expected (could signal an assembly problem or a genuinely unusual sample). |
| `contaminant` | sample | The closest published match belongs to an unrelated group, which can indicate contamination. |

The last two checks (`low_identity` and `contaminant`) rely on looking
up a reference match on NCBI during assembly. If that lookup was
unavailable at the time, those two checks are simply skipped and the
others still apply.

## How the checks are tuned per group

Some checks, especially the expected lengths, depend on the group of
organisms. RiboPilot ships two settings:

- **Generic animal** (the default): deliberately wide length ranges that
  suit a broad range of animals.
- **Fish**: tighter, fish-specific ranges for the ribosomal genes, with
  the spacers kept wide because they vary even between close relatives.

For reference, the generic animal length ranges are:

| Feature | Typical length (bases)                                           |
|---------|------------------------------------------------------------------|
| 18S     | 1600 - 2100                                                      |
| 5.8S    | 140 - 180                                                        |
| 28S     | 2000 - 6000 (the recovered portion varies with sequencing depth) |
| ITS1    | 120 - 1600                                                       |
| ITS2    | 120 - 1600                                                       |

A `length` flag simply means a feature fell outside these ranges. That
is common and often fine, especially for the fast-evolving spacers; it
is a prompt to look, not a verdict.

## Adding tuning for a new group

If you work on a group that would benefit from its own length ranges,
that can be added. This is a small change to RiboPilot’s settings rather
than new analysis code. If you are comfortable with R you can add it
yourself (see `rrna_params()` in the source), but you are also very
welcome to [open an
issue](https://github.com/dmacguigan/RiboPilot/issues) describing your
group, and we can help set it up. In the meantime, the generic animal
setting works for any group; it is just more lenient.
