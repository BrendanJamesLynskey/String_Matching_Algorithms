# String Matching Algorithms

**Computer Science Fundamentals Series**

KMP · Rabin-Karp · Boyer-Moore · Suffix arrays · Aho-Corasick · Pattern matching

*Mid-level software engineer track -- 20 slides*

---

## Table of Contents

1. [The String Matching Problem](#slide-02--the-string-matching-problem)
2. [Naive String Matching](#slide-03--naive-string-matching)
3. [KMP -- Intuition & Failure Function](#slide-04--kmp--intuition--failure-function)
4. [KMP -- Algorithm Walkthrough](#slide-05--kmp--algorithm-walkthrough)
5. [Rabin-Karp -- Rolling Hash](#slide-06--rabin-karp--rolling-hash)
6. [Rabin-Karp -- Spurious Hits & Multi-Pattern](#slide-07--rabin-karp--spurious-hits--multi-pattern)
7. [Boyer-Moore -- Bad Character Rule](#slide-08--boyer-moore--bad-character-rule)
8. [Boyer-Moore -- Good Suffix Rule](#slide-09--boyer-moore--good-suffix-rule)
9. [Z-Algorithm](#slide-10--z-algorithm)
10. [Aho-Corasick -- Multi-Pattern Matching](#slide-11--aho-corasick--multi-pattern-matching)
11. [Suffix Arrays](#slide-12--suffix-arrays)
12. [Suffix Trees](#slide-13--suffix-trees)
13. [Regular Expression Matching](#slide-14--regular-expression-matching)
14. [Approximate String Matching -- Edit Distance](#slide-15--approximate-string-matching--edit-distance)
15. [Longest Common Substring](#slide-16--longest-common-substring)
16. [Application -- Text Editors & Search](#slide-17--application--text-editors--search)
17. [Application -- DNA Sequencing](#slide-18--application--dna-sequencing)
18. [Application -- Plagiarism Detection](#slide-19--application--plagiarism-detection)
19. [Complexity Comparison & Summary](#slide-20--complexity-comparison--summary)

---

## Slide 02 -- The String Matching Problem

### Definition

Given a text `T` of length `n` and a pattern `P` of length `m`, find all occurrences of `P` in `T`.

- **Input:** text `T[0..n-1]`, pattern `P[0..m-1]`
- **Output:** all positions `i` where `T[i..i+m-1] == P`

### Why it matters

String matching is one of the most fundamental problems in computer science. Every time you press Ctrl+F, run `grep`, query a database with `LIKE`, or align a DNA sequence, a string matching algorithm is at work.

### Terminology

| Term | Meaning |
|------|---------|
| **Text** | The string being searched (haystack) |
| **Pattern** | The string being searched for (needle) |
| **Alphabet** | The set of valid characters, size `|Σ|` |
| **Shift** | An alignment position of pattern against text |
| **Valid shift** | A shift where every character matches |

> The naive approach checks every shift. Efficient algorithms skip provably invalid shifts.

---

## Slide 03 -- Naive String Matching

### The brute-force approach

For each position `i` in the text, compare `P[0..m-1]` against `T[i..i+m-1]` character by character.

```
NAIVE-STRING-MATCH(T, P):
  n = len(T), m = len(P)
  for i = 0 to n - m:
    if T[i..i+m-1] == P[0..m-1]:
      report match at i
```

### Complexity

- **Worst case:** `O(n * m)` -- e.g. `T = "aaaaaa"`, `P = "aaab"`
- **Best case:** `O(n)` -- first character of pattern rarely matches
- **Average case:** `O(n)` for random text over a large alphabet
- **Space:** `O(1)` -- no preprocessing

### When it is good enough

- Short patterns (`m < 10`) on moderate text
- Large alphabets (English text) where mismatches occur early
- One-off searches where preprocessing overhead is not justified

> The naive algorithm discards all information gained from a partial match. Every advanced algorithm exploits that wasted knowledge.

---

## Slide 04 -- KMP -- Intuition & Failure Function

### Core insight

When a mismatch occurs at `P[j]` after matching `P[0..j-1]`, the naive algorithm restarts from scratch. KMP asks: what is the longest proper prefix of `P[0..j-1]` that is also a suffix? We can resume matching from there.

### The failure function (prefix table)

`π[j]` = length of the longest proper prefix of `P[0..j]` that is also a suffix.

```
Pattern: A B A B A C
Index:   0 1 2 3 4 5
π:       0 0 1 2 3 0
```

### Building the prefix table -- O(m)

```
COMPUTE-PREFIX(P):
  m = len(P)
  π[0] = 0
  k = 0
  for q = 1 to m - 1:
    while k > 0 and P[k] != P[q]:
      k = π[k - 1]
    if P[k] == P[q]:
      k = k + 1
    π[q] = k
  return π
```

> The prefix table is the heart of KMP. It encodes the internal structure of the pattern so that no text character is ever re-examined.

---

## Slide 05 -- KMP -- Algorithm Walkthrough

### The algorithm

```
KMP-SEARCH(T, P):
  n = len(T), m = len(P)
  π = COMPUTE-PREFIX(P)
  q = 0                       // characters matched so far
  for i = 0 to n - 1:
    while q > 0 and P[q] != T[i]:
      q = π[q - 1]            // fall back in pattern
    if P[q] == T[i]:
      q = q + 1
    if q == m:
      report match at i - m + 1
      q = π[q - 1]            // look for next match
```

### Complexity

| Metric | Value |
|--------|-------|
| **Preprocessing** | `O(m)` -- build prefix table |
| **Matching** | `O(n)` -- each text character examined at most twice |
| **Total** | `O(n + m)` |
| **Space** | `O(m)` -- prefix table |

### Trace example

```
Text:    A B A B A B A C A B
Pattern: A B A B A C
π:       0 0 1 2 3 0

i=4: matched ABABA, mismatch at P[5]='C' vs T[5]='B'
     q = π[4] = 3, resume comparing P[3] with T[5]
     → skipped re-examining T[0..2]
```

> KMP guarantees exactly `O(n)` comparisons in the worst case -- no algorithm can do fewer for single-pattern exact matching.

---

## Slide 06 -- Rabin-Karp -- Rolling Hash

### Idea

Instead of comparing characters, compute a hash of each m-length window of text and compare it to the hash of the pattern. A *rolling hash* lets us update the hash in `O(1)` per shift.

### Rolling hash (polynomial)

```
hash(S[i..i+m-1]) = Σ S[i+j] * d^(m-1-j)  mod q

where d = |Σ| (alphabet size), q = large prime
```

Sliding the window:

```
hash(S[i+1..i+m]) = (hash(S[i..i+m-1]) - S[i] * d^(m-1)) * d + S[i+m]  mod q
```

### Algorithm sketch

```
RABIN-KARP(T, P):
  hp = hash(P[0..m-1])
  ht = hash(T[0..m-1])
  for i = 0 to n - m:
    if ht == hp:
      if T[i..i+m-1] == P:    // verify to eliminate false positive
        report match at i
    ht = roll(ht, T[i], T[i+m])
```

> The power of Rabin-Karp is its simplicity and its natural extension to multi-pattern search -- hash each pattern, store in a set, check every window against the set.

---

## Slide 07 -- Rabin-Karp -- Spurious Hits & Multi-Pattern

### Spurious hits

Two different strings can have the same hash (collision). Each collision triggers a full `O(m)` character comparison. With a good prime `q`, the expected number of spurious hits is `O(n/q)`.

### Complexity

| Case | Time |
|------|------|
| **Best / average** | `O(n + m)` |
| **Worst case** | `O(n * m)` -- many collisions (pathological hash) |
| **Space** | `O(1)` beyond input (single pattern) |

### Multi-pattern variant

To search for `k` patterns simultaneously:

1. Hash all `k` patterns and store in a hash set
2. Slide window over text; check each window hash against the set
3. Expected time: `O(n + km)` -- far better than running `k` separate searches

### Practical tips

- Use double hashing (two independent hash functions) to reduce false positives to near zero
- Rabin-Karp is the basis of many plagiarism detectors (compare document fingerprints)
- Works well when pattern length is fixed and known

> Rabin-Karp trades deterministic guarantees for practical speed and multi-pattern flexibility.

---

## Slide 08 -- Boyer-Moore -- Bad Character Rule

### Why Boyer-Moore?

Boyer-Moore is often the fastest single-pattern algorithm in practice. It scans the pattern *right to left* and uses two heuristics to skip large portions of text.

### Bad character rule

When a mismatch occurs at text character `c` and pattern position `j`:
- If `c` does not appear in the pattern, shift the pattern past the mismatch entirely
- If `c` appears at position `k < j` in the pattern, align that occurrence with the text

```
Text:    ...X E R C I S E...
Pattern:    E X A M P L E
            ←───────── comparing right to left
Mismatch: T[3]='R' vs P[3]='M'
'R' not in pattern → shift pattern 4 positions right
```

### Bad character table

Precompute for each character in the alphabet: its rightmost position in the pattern (or -1 if absent).

```
BUILD-BAD-CHAR(P):
  for each c in Σ: bc[c] = -1
  for j = 0 to m - 1:
    bc[P[j]] = j
  return bc
```

> The bad character rule alone gives sub-linear average performance on large alphabets -- many characters in the text simply do not appear in the pattern.

---

## Slide 09 -- Boyer-Moore -- Good Suffix Rule

### Good suffix rule

When characters `P[j+1..m-1]` have matched but `P[j]` mismatches:
- Find the rightmost occurrence of the matched suffix elsewhere in the pattern
- If no full occurrence exists, find the longest prefix of the pattern that matches a suffix of the matched portion
- Shift accordingly

### Combined shift

Boyer-Moore uses the *maximum* of the bad character and good suffix shifts at each mismatch.

```
shift = max(bad_char_shift, good_suffix_shift)
```

### Complexity

| Metric | Value |
|--------|-------|
| **Preprocessing** | `O(m + |Σ|)` |
| **Best case** | `O(n/m)` -- sub-linear! |
| **Worst case** | `O(n * m)` -- but `O(n)` with Galil's improvement |
| **Average** | `O(n/m)` for large alphabets |

### Boyer-Moore in practice

- Default algorithm in GNU `grep`
- Faster than KMP for long patterns on natural-language text
- The `O(n/m)` best case means longer patterns can actually be *faster* to find

> Boyer-Moore is the rare algorithm where increasing the pattern length improves performance.

---

## Slide 10 -- Z-Algorithm

### Z-array definition

For a string `S`, `Z[i]` = length of the longest substring starting at `S[i]` that matches a prefix of `S`.

```
S:      a a b x a a b
Z:      - 1 0 0 3 1 0
```

### Algorithm -- O(n)

Maintain a Z-box `[l, r]` representing the rightmost interval matching a prefix. For each position `i`:
- If `i > r`, compute `Z[i]` naively and update `[l, r]`
- If `i <= r`, use previously computed values to initialise `Z[i]`, then extend if needed

### String matching with Z-algorithm

Concatenate `P$T` (where `$` is not in the alphabet) and compute the Z-array. Any position `i` where `Z[i] == m` is a match.

```
P = "aab", T = "aabxaab"
S = "aab$aabxaab"
Z = [- 1 0 0 3 1 0 0 3 1 0]
         matches at i=4 and i=8 → text positions 0 and 4
```

### Complexity

| Metric | Value |
|--------|-------|
| **Time** | `O(n + m)` |
| **Space** | `O(n + m)` |

> The Z-algorithm is conceptually simpler than KMP and equally powerful. It is also the foundation for many suffix-based constructions.

---

## Slide 11 -- Aho-Corasick -- Multi-Pattern Matching

### Problem

Given `k` patterns and a text of length `n`, find all occurrences of all patterns in one pass.

### Construction

1. **Build a trie** from all `k` patterns
2. **Add failure links** (analogous to KMP's failure function) -- on mismatch, follow the failure link to the longest proper suffix that is also a prefix of some pattern
3. **Add output links** -- chain patterns that end at the same trie node via suffixes

### Search -- O(n + m + z)

Feed the text through the automaton character by character. At each state, follow output links to report all matching patterns. `z` = total number of matches reported.

```
Patterns: {he, she, his, hers}
Text:     "ushers"

Matches found: she (at 1), he (at 2), hers (at 2)
```

### Complexity

| Metric | Value |
|--------|-------|
| **Build** | `O(m · |Σ|)` where `m` = total pattern length |
| **Search** | `O(n + z)` -- linear in text + output |
| **Space** | `O(m · |Σ|)` for the automaton |

> Aho-Corasick is the multi-pattern generalisation of KMP. It powers network intrusion detection (Snort), virus scanners, and `fgrep`.

---

## Slide 12 -- Suffix Arrays

### Definition

A suffix array `SA` for a string `S` of length `n` is a sorted array of all suffix indices `[0, 1, ..., n-1]` ordered by the lexicographic order of the suffixes they represent.

```
S = "banana$"
Suffixes sorted:     SA:
$                    6
a$                   5
ana$                 3
anana$               1
banana$              0
na$                  4
nana$                2
```

### Pattern search via binary search

To find pattern `P` in `S`, binary search the suffix array. Each comparison is `O(m)`, so search is `O(m log n)`.

With an **LCP array** (longest common prefix between adjacent suffixes), this improves to `O(m + log n)`.

### Construction

| Algorithm | Time | Space |
|-----------|------|-------|
| **Naive sort** | `O(n² log n)` | `O(n)` |
| **Prefix doubling** (Karp-Miller-Rosenberg) | `O(n log n)` | `O(n)` |
| **SA-IS / DC3** | `O(n)` | `O(n)` |

> Suffix arrays use 4-8x less memory than suffix trees while supporting most of the same queries. They are the preferred structure in modern bioinformatics tools.

---

## Slide 13 -- Suffix Trees

### Definition

A suffix tree for string `S` is a compressed trie of all suffixes of `S$`. Every internal node (except root) has at least two children. Edge labels are substrings of `S`.

### Key properties

- Contains `n` leaves (one per suffix)
- At most `n - 1` internal nodes
- Can be built in `O(n)` time (Ukkonen's algorithm)
- Every substring of `S` corresponds to a path from the root

### Operations

| Operation | Time |
|-----------|------|
| Pattern search | `O(m)` |
| Longest repeated substring | `O(n)` -- deepest internal node |
| Longest common substring (two strings) | `O(n + m)` -- generalised suffix tree |
| Number of distinct substrings | `O(n)` -- sum of edge lengths |
| Shortest unique substring | `O(n)` |

### Trade-offs

- **Advantage:** `O(m)` pattern matching (no log factor), rich query support
- **Disadvantage:** 10-20x memory overhead per character; complex to implement
- **In practice:** suffix arrays + LCP arrays have largely replaced suffix trees due to cache efficiency and lower memory

> Suffix trees are the Swiss army knife of string algorithms -- nearly every string problem has an `O(n)` suffix tree solution.

---

## Slide 14 -- Regular Expression Matching

### From regex to automaton

1. Parse regex into an **NFA** (Thompson's construction) -- `O(m)` states for a pattern of length `m`
2. Simulate NFA on input or convert to **DFA** for faster matching

### NFA simulation -- Thompson's approach

Maintain a *set* of active states. For each input character, compute the next set of states via transitions. Match if any final state is in the set.

- **Time:** `O(n * m)` -- `n` characters, up to `m` states per step
- **Space:** `O(m)`
- **Guarantee:** always linear in `n * m`, no pathological backtracking

### DFA approach

Convert NFA to DFA (subset construction). Each input character requires exactly one state transition.

- **Matching time:** `O(n)` -- one transition per character
- **Construction:** up to `O(2^m)` states in worst case (exponential blowup)
- **Practical:** lazy DFA construction builds states on demand, caching frequently used ones

### Backtracking engines (PCRE, Python `re`)

Most real-world regex engines use backtracking, which is `O(2^n)` worst case. Patterns like `(a+)+b` cause catastrophic backtracking.

> Use RE2, Rust's `regex`, or Go's `regexp` for guaranteed linear-time matching. Avoid backtracking engines on untrusted input.

---

## Slide 15 -- Approximate String Matching -- Edit Distance

### Edit distance (Levenshtein distance)

The minimum number of single-character operations (insert, delete, substitute) to transform string `A` into string `B`.

```
edit_distance("kitten", "sitting") = 3
  kitten → sitten (substitute k→s)
  sitten → sittin (substitute e→i)
  sittin → sitting (insert g)
```

### Dynamic programming -- O(nm)

```
DP[i][j] = min(
  DP[i-1][j] + 1,        // delete from A
  DP[i][j-1] + 1,        // insert into A
  DP[i-1][j-1] + cost    // substitute (cost=0 if A[i]==B[j])
)
```

### Variants

| Variant | Allowed operations |
|---------|--------------------|
| **Levenshtein** | Insert, delete, substitute |
| **Damerau-Levenshtein** | + transposition of adjacent characters |
| **Hamming distance** | Substitution only (equal-length strings) |
| **Longest common subsequence** | Insert and delete only (no substitution) |

### Applications

- Spell checkers and autocorrect
- Fuzzy search (find strings within edit distance `k`)
- DNA/protein sequence alignment (with weighted operations)

> Approximate matching with `O(nm)` DP is the foundation of bioinformatics alignment algorithms like Smith-Waterman and Needleman-Wunsch.

---

## Slide 16 -- Longest Common Substring

### Problem

Given two strings `A` (length `n`) and `B` (length `m`), find the longest string that is a contiguous substring of both.

### DP approach -- O(nm)

```
DP[i][j] = DP[i-1][j-1] + 1   if A[i] == B[j]
         = 0                    otherwise

Answer = max(DP[i][j]) for all i, j
```

### Suffix array approach -- O((n+m) log(n+m))

1. Concatenate `A$B#` (with unique separators)
2. Build suffix array and LCP array
3. Answer = maximum LCP between adjacent suffixes where one belongs to `A` and the other to `B`

### Suffix tree approach -- O(n + m)

Build a generalised suffix tree for both strings. The deepest internal node with descendants from both strings gives the longest common substring.

### Comparison

| Method | Time | Space |
|--------|------|-------|
| **DP** | `O(nm)` | `O(nm)` or `O(m)` with rolling |
| **Suffix array** | `O((n+m) log(n+m))` | `O(n + m)` |
| **Suffix tree** | `O(n + m)` | `O(n + m)` (large constant) |

> Longest common substring is a building block for diff tools, plagiarism detection, and DNA sequence comparison.

---

## Slide 17 -- Application -- Text Editors & Search

### Find & replace

Every text editor uses string matching. The choice of algorithm depends on the use case:

- **Short patterns, interactive:** naive or simplified Boyer-Moore -- low latency, no preprocessing overhead
- **Long patterns, large files:** full Boyer-Moore -- sub-linear scanning
- **Regex search:** Thompson NFA or lazy DFA (as in Vim, VS Code, ripgrep)

### grep and ripgrep

| Tool | Algorithm | Key feature |
|------|-----------|-------------|
| `grep -F` (fgrep) | Aho-Corasick | Multi-pattern literal search |
| `grep -E` | DFA-based regex | POSIX extended regex |
| GNU `grep` | Boyer-Moore + regex | Literal prefix optimisation |
| `ripgrep` | Aho-Corasick + lazy DFA | Parallel directory walk, `.gitignore`-aware |

### Full-text search engines

Elasticsearch, Lucene, and Typesense build **inverted indexes** mapping terms to document IDs. At query time, they intersect posting lists -- a form of multi-pattern matching at scale.

> ripgrep demonstrates modern string matching: literal optimisations (Aho-Corasick for multiple literals, memchr for single bytes) layered with regex when needed.

---

## Slide 18 -- Application -- DNA Sequencing

### The alignment problem

DNA sequencing produces millions of short reads (50-300 base pairs). Each must be aligned to a reference genome (3 billion base pairs for humans). Speed is critical.

### BWT-based aligners

Modern aligners (BWA, Bowtie2) use the **Burrows-Wheeler Transform** and **FM-index** -- compressed representations of suffix arrays that support `O(m)` exact matching with minimal memory.

### Approximate alignment

DNA has mutations, so exact matching is insufficient. Tools combine:
- **Seed-and-extend:** find exact short matches (seeds), then extend with edit-distance DP
- **Smith-Waterman:** local alignment via DP with affine gap penalties
- **BLAST:** heuristic seed-based search with statistical significance scoring

### Scale

| Metric | Value |
|--------|-------|
| Human genome | ~3.2 billion base pairs |
| Typical sequencing run | 100M-1B reads |
| Alignment speed (BWA-MEM2) | ~60M reads/hour |

> String matching algorithms are the computational backbone of genomics. Without BWT/FM-index, whole-genome sequencing would be computationally infeasible.

---

## Slide 19 -- Application -- Plagiarism Detection

### Document fingerprinting

Break documents into overlapping k-grams (k consecutive characters or words), hash each, and compare fingerprint sets between documents.

### Winnowing algorithm

Select a subset of k-gram hashes using a sliding window minimum. Guarantees:
- Any match of length `≥ t` is detected (where `t` = guarantee threshold)
- Fingerprint density is bounded -- not every hash is stored

### MOSS (Measure of Software Similarity)

Stanford's plagiarism detector for code. Uses winnowing with language-aware tokenisation (ignoring variable names, whitespace, comments).

### Similarity metrics

| Metric | Formula |
|--------|---------|
| **Jaccard similarity** | `|A ∩ B| / |A ∪ B|` |
| **Containment** | `|A ∩ B| / |A|` |
| **Cosine similarity** | `A · B / (|A| · |B|)` on TF-IDF vectors |

### Pipeline

1. Tokenise and normalise (lowercase, remove stopwords or comments)
2. Generate k-gram hashes (Rabin-Karp rolling hash)
3. Select fingerprints (winnowing)
4. Compare fingerprint sets across document pairs

> Rabin-Karp's rolling hash is the engine behind most modern plagiarism detection systems, from academic tools like Turnitin to code similarity checkers.

---

## Slide 20 -- Complexity Comparison & Summary

### Algorithm comparison

| Algorithm | Preprocessing | Matching | Space | Best for |
|-----------|--------------|----------|-------|----------|
| **Naive** | `O(0)` | `O(nm)` | `O(1)` | Short patterns, one-off search |
| **KMP** | `O(m)` | `O(n)` | `O(m)` | Streaming text, worst-case guarantee |
| **Rabin-Karp** | `O(m)` | `O(n)` avg | `O(1)` | Multi-pattern, fingerprinting |
| **Boyer-Moore** | `O(m + |Σ|)` | `O(n/m)` avg | `O(m + |Σ|)` | Long patterns, large alphabets |
| **Z-Algorithm** | `O(m)` | `O(n)` | `O(n + m)` | Simple implementation, competitive |
| **Aho-Corasick** | `O(Σm)` | `O(n + z)` | `O(Σm)` | Many patterns simultaneously |
| **Suffix Array** | `O(n)` | `O(m log n)` | `O(n)` | Many queries on same text |
| **Suffix Tree** | `O(n)` | `O(m)` | `O(n)` | Rich substring queries |

### Key takeaways

- No single algorithm wins everywhere -- the choice depends on alphabet size, pattern length, number of patterns, and whether the text is static or streaming
- KMP and Z-algorithm give `O(n + m)` worst-case guarantees for single-pattern matching
- Boyer-Moore is fastest in practice for single long patterns on natural text
- Aho-Corasick is the standard for multi-pattern matching
- Suffix arrays and trees trade preprocessing time for fast repeated queries
- Approximate matching (edit distance) is fundamentally `O(nm)` but can be accelerated with bit-parallelism and filtering
- Real systems combine multiple algorithms -- ripgrep uses Aho-Corasick, memchr, and lazy DFA depending on the query

### Recommended reading

| Source | Description |
|--------|------------|
| **CLRS** | *Introduction to Algorithms* -- chapters on string matching cover KMP, Rabin-Karp, and finite automata |
| **Gusfield** | *Algorithms on Strings, Trees, and Sequences* -- the definitive reference for suffix trees, edit distance, and biological sequence analysis |
| **Sedgewick & Wayne** | *Algorithms* -- practical implementations of KMP, Boyer-Moore, and regex matching |
| **Navarro & Raffinot** | *Flexible Pattern Matching in Strings* -- advanced techniques and bit-parallel algorithms |
| **cp-algorithms.com** | [String Processing](https://cp-algorithms.com/string/) -- excellent implementations of Z, KMP, Aho-Corasick, suffix array |
