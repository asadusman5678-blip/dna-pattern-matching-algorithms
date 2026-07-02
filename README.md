# DNA Pattern Matching Algorithms

Implementations of classic exact and approximate string-matching algorithms applied to real genomic and sequencing data — comparing brute-force search against the optimized Boyer-Moore algorithm, and analyzing read quality from a FASTQ file.

## Overview

Finding a short pattern (e.g. a primer, motif, or read) inside a much longer reference sequence is a fundamental operation in bioinformatics. This repo walks through:

1. **Naive matching** — straightforward brute-force search, including reverse-complement-aware search and approximate (mismatch-tolerant) search.
2. **Boyer-Moore matching** — a preprocessing-based algorithm that skips alignments that can't possibly match, verified against naive matching on a real excerpt of human chromosome 1.
3. **FASTQ quality analysis** — using Phred quality scores to detect a systematically low-quality sequencing cycle.

## Notebooks

### `01_naive_matching_and_read_quality.ipynb`

| Function | What it does |
|---|---|
| `naive(p, t)` | Brute-force search for exact occurrences of pattern `p` in text `t`. Checks every possible alignment, character by character. |
| `reverseComplement(s)` | Returns the reverse complement of a DNA string (A↔T, C↔G) — needed because a pattern can match either strand of DNA. |
| `naive_with_rc(p, t)` | Runs `naive` on both the pattern and its reverse complement, so matches on either strand are found. |
| `naive_2mm(p, t)` | Like `naive`, but allows up to 2 character mismatches per alignment — approximate matching. |
| `readGenome(filename)` | Reads a FASTA file into a single sequence string, skipping the header line. |
| `readFastq(filename)` | Parses a FASTQ file into two lists: read sequences and their per-base quality strings. |
| `phred33ToQ(qual)` | Converts a Phred+33 quality string into a list of numeric quality scores. |

**What the notebook demonstrates:**
- Searching the lambda virus genome (`lambda_virus.fa`) for the pattern `AGTCGA`, including its reverse complement.
- Approximate matching of `AGGAGGTT` allowing up to 2 mismatches — 215 matches found, versus far fewer for exact matching.
- Analyzing `ERR037900_1.first1000.fastq` to find a sequencing cycle (base position across reads) with a systematically low average quality score — flags miscalibrated or unreliable sequencer cycles.

### `02_boyer_moore_matching.ipynb`

| Function | What it does |
|---|---|
| `naive_with_counts(p, t)` | Same as `naive`, but also tracks the number of alignments tried and character comparisons made. |
| `boyer_moore_with_counts(p, p_bm, t)` | Boyer-Moore search using the bad character rule and good suffix rule (via the `BoyerMoore` class in `bm_preproc.py`) to skip ahead on mismatches, tracking alignments and comparisons for direct comparison against naive. |
| `read_fasta(filename)` | Reads a FASTA file into a single sequence string. |

**What the notebook demonstrates:**

Searching a 48bp pattern against an excerpt of human chromosome 1 (GRCh38), comparing efficiency:

| Algorithm    | Alignments | Character Comparisons |
|--------------|-----------:|-----------------------:|
| Naive        | 799,954    | 984,143                 |
| Boyer-Moore  | 127,974    | 165,191                 |

Boyer-Moore needs **~6x fewer alignments and comparisons** to find the same match — the preprocessing overhead pays off fast on real genomic-scale data.

### `03_edit_distance_and_overlap_graphs.ipynb`

| Function | What it does |
|---|---|
| `global_edit_distance(x, y)` | Computes the standard (Levenshtein) edit distance between two full sequences using dynamic programming — the minimum number of insertions, deletions, and substitutions to turn one into the other. |
| `approximate_substring_match(x, y)` | A modified edit-distance DP that allows the match to start anywhere in `y` without penalty, finding the best approximate occurrence of pattern `x` inside a longer text `y`. |
| `overlap(a, b, min_length)` | Returns the length of the longest suffix of `a` that matches a prefix of `b`, if at least `min_length` — the basic building block for detecting overlaps between reads. |
| `read_fastq(filename)` | Parses a FASTQ file into a list of read sequences. |
| `read_fasta(filename)` | Reads a FASTA file into a single sequence string. |

**What the notebook demonstrates:**
- **Global edit distance**: aligning two full sequences end-to-end (e.g. distance of 11 between two 8bp/16bp test sequences).
- **Approximate substring matching**: finding where a pattern best fits inside a much longer text, allowing mismatches — including scanning a full 800,000-base excerpt of chromosome 1 for a 16bp pattern, finding a best match with only 2 edits.
- **Overlap graph construction for genome assembly**: using a k-mer index (k=30) to efficiently find which sequencing reads overlap with each other by at least `k` bases — the first step in assembling reads into a full genome. On `ERR266411_1.for_asm.fastq`, this finds 904,746 overlap edges across 7,161 reads with at least one outgoing edge. The k-mer index avoids the O(n²) cost of comparing every read against every other read directly.

## Repository structure

```
dna-pattern-matching-algorithms/
├── 01_naive_matching_and_read_quality.ipynb
├── 02_boyer_moore_matching.ipynb
├── 03_edit_distance_and_overlap_graphs.ipynb
├── bm_preproc.py              # Boyer-Moore preprocessing (bad char & good suffix tables)
├── kmer_index.py               # k-mer indexing helper for approximate matching
├── lambda_virus.fa             # ~49KB reference genome used in notebook 1
├── ERR037900_1.first1000.fastq # first 1000 reads, used for quality analysis
├── chr1.GRCh38.excerpt.fasta   # excerpt of human chr1, used in notebooks 2 & 3
├── ERR266411_1.for_asm.fastq   # sequencing reads used for overlap graph in notebook 3
├── requirements.txt
├── .gitignore
└── README.md
```

## Setup & usage

```bash
git clone https://github.com/asadusman/dna-pattern-matching-algorithms.git
cd dna-pattern-matching-algorithms
pip install -r requirements.txt
jupyter notebook
```

Open either notebook and run all cells top to bottom. Both `bm_preproc.py` and `kmer_index.py` must stay in the same directory as the notebooks — they're imported directly, not installed as packages. All data files are small excerpts, not full genomes/datasets.

## Background

These implementations were developed while studying classic and modern algorithms for genomic sequence alignment — the same underlying ideas (exact matching, indexing, approximate matching) used in production tools like BWA and Bowtie for real-world read alignment.

