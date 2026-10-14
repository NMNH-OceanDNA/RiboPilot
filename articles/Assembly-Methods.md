# Choosing an assembly method

RiboPilot can build the ribosomal repeat from your reads in three ways.
[GetOrganelle](https://github.com/Kinggerm/GetOrganelle) is the default.
The other two use [SPAdes](https://ablab.github.io/spades/) and are
there for samples GetOrganelle cannot assemble: **baited** assembles
only the reads that match the reference database, and **de-novo**
assembles every read.

The method is chosen per sample, under **Assembler** in the **Assembly
Options** window of the Assemble step. Everything after assembly
(quality checks, annotation, export) works the same whichever you pick.

If you already have an assembly for a sample, you can skip assembly
altogether. See [Already have an
assembly?](https://nmnh-oceandna.github.io/RiboPilot/articles/Your-Own-Project.html#already-have-an-assembly).

## At a glance

|  | GetOrganelle (default) | SPAdes, baited | SPAdes, de-novo |
|----|----|----|----|
| **How it works** | Recruits ribosomal reads with the seed database, then extends and assembles them over several rounds | Maps reads to the seed database, then assembles only the pairs that match | Assembles all reads, then keeps the contigs that match the label database |
| **Uses the reference database for** | Recruiting reads and picking ribosomal contigs | Recruiting reads and picking ribosomal contigs | Picking ribosomal contigs only |
| **Speed and memory** | Moderate | Fastest, lowest memory | Slowest, highest memory |
| **Best for** | Most samples | Samples where GetOrganelle fails or returns only short contigs | Samples where both of the others fail, or a cross-check on a doubtful result |
| **Watch out for** | Can fail to extend when coverage is low or uneven | Only recovers what resembles the reference, so contigs can stop short in the variable spacer regions | Runtime and memory; more contigs to review |

## How long each one takes

One sample from the test project, run three times on the same desktop
with each assembler, gives a sense of scale. The sample (SRR14163830, a
pufferfish) has 1.4 million read pairs of 2 x 150 bp. Each run had 6
CPUs and up to 24 GB of memory. The times and memory are for the
assembly step alone; the rest of the workflow (read cleanup, quality
checks, annotation) added about two minutes in every case.

|  | GetOrganelle | SPAdes, baited | SPAdes, de-novo |
|----|----|----|----|
| **Assembly time** | 1 min 20 s | 19 s | 38 min 30 s |
| **Peak memory** | 2.9 GB | 0.3 GB | 4.1 GB |
| **Contigs** | 1 (8,960 bp) | 3 (5,031, 2,223, and 564 bp) | 11 (longest 6,159 bp) |

GetOrganelle gave the whole repeat in one piece. Baited SPAdes was the
fastest by far, but split the same repeat into three contigs. De-novo
SPAdes took about thirty times longer than GetOrganelle, and about a
hundred times longer than baited, and left more fragments to review.
Larger read sets scale these numbers up, and de-novo scales up the most.

**Where to start.** Use GetOrganelle for the whole project. Switch a
sample to baited SPAdes when it fails. Keep de-novo SPAdes for the few
samples that neither of the others can assemble.

## GetOrganelle

GetOrganelle starts from a **seed database** of known ribosomal
sequences, recruits the reads that match, and extends outward from them
round after round. A **label database** then picks out which parts of
the result are ribosomal. Both databases come from the [Reference
step](https://nmnh-oceandna.github.io/RiboPilot/articles/Reference-Databases.md),
so a database built for the right group matters more here than any
option.

The **GetOrganelle options** box carries its command-line flags. The
[GetOrganelle wiki](https://github.com/Kinggerm/GetOrganelle/wiki)
explains them.

## SPAdes, baited

Baited assembly is the first thing to try on a sample that failed. Reads
are mapped to the seed database with [Bowtie
2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), only the
read pairs that match are kept, and SPAdes assembles that small pool.
Contigs that match the label database are kept; the rest are discarded.

Because the read pool is small, this is quick and light on memory. The
cost is that a read must resemble the reference to be kept. The genes
(18S, 5.8S, 28S) are conserved and recruit well. The spacers between
them vary much more, so on a sample distant from everything in the
database a contig can end early or break at a spacer.

Two boxes control it:

- **bowtie2 options.** How reads are matched to the seed database. The
  default is `--local`.
- **SPAdes options.** The default is `--isolate`. If SPAdes rejects the
  reads with these flags, the run is retried once with no extra flags.

## SPAdes, de-novo

De-novo assembly skips the baiting step and gives SPAdes every read. The
reference database is used only afterwards, to keep the contigs that
match the label database. Nothing is filtered up front, so this can
recover a repeat whose reads the seed database fails to recruit.

It is also by far the most expensive method, since it assembles all of
the sequencing reads. Raise **Memory (GB)** in the Assembly Options
window if a sample fails for lack of memory. Only the **SPAdes options**
box applies.

## Mixing methods in one project

Since assembly options are saved as named sets, a project can use
GetOrganelle for most samples and SPAdes for the difficult ones. Open a
sample’s **Assembly Options** cell, type a new **Parameter set name**,
tick **Edit**, choose the assembler, and press **Save**. The [Assemble
walkthrough](https://nmnh-oceandna.github.io/RiboPilot/articles/Test-Project-Assemble.html#rescue-a-sample-that-failed)
shows this on a failed test sample.

A GenBank submission records one assembly method. If an export group
mixes assemblers, the [Export
step](https://nmnh-oceandna.github.io/RiboPilot/articles/Test-Project-Export.md)
offers to split it.
