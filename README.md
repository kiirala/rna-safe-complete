# Safe and Complete RNA Secondary Structure Prediction

The *secondary structure*, or *folding*, of
[RNA](https://en.wikipedia.org/wiki/RNA) is how the bases of a RNA
sequence pair with each other.

Safe and Complete is a concept related to algorithms that may produce a
large amount of equally good solutions. Here *Safe* refers to
finding parts of solutions that are common to all good solutions,
and *Complete* refers to finding all such parts.

This package implements three different *maximum pairs* predictors for
RNA secondary structure. Each of these attempts to create a folding
with as many as possible pairs, which is not always the energically
most favourable folding, but is a reasonable approximation and fairly
straightforward to implement and extend.

* Nussinov / Zuker algorithm. The Nussinov algorithm finds the maximal
  number of pairs for a sequence, and a single example of such
  folding. The Zuker algorithm is extension of this, and finds multiple
  examples of optimal foldings.
* Wuchty algorithm: this algorithm finds the maximal number of pairs
  and all distinct foldings with this many pairs.
* Safe and Complete algorithm: this is a novel extension of Nussinov
  and Wuchty algorithms. This algorithm finds maximal number of
  pairs. It creates a sample folding and tells which parts of it
  are such that they appear in all optimal foldings. It can also show
  all optimal foldings or create a matrix of how may times each base
  pair appears in them.

There are multiple programs provided in this repository:

* **rnafolding**: Experimental program that runs multiple different
  algorithms and checks their results against each other and against
  some sanity checks
  * Can read a single strand from FASTA file or large number of
    strands from Sprinzl tRNA database file or from
    dot-bracket files downloaded from STRAND database
  * Analyzes either single sequence outputting a large amount of
    information to console, or a number of sequences outputting
    statistics about them
  * Can output the folding details of each analyzed sequence to
    JSON file
* **comparesafety** Fast and memory-efficient program for computing
  safety, built to compare the efficiency of trivial safety algorithm
  and the dynamic programming version
* **trivialsafety** Compute safety metrics from the output of the
  RNAsubopt program of the ViennaRNA package, using a trivial
  algorithm
* **dbtofasta** Convert Sprinzl tRNA database file into set of FASTA
  files
* **fastadump** Dump the internal representation of a sequence read
  from a FASTA file into a JSON file
