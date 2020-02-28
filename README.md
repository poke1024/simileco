# Simileco

A header-only templated C++17-library for fast and versatile sequence alignment computations.

Most libraries for sequence alignments out there are hard-coded to a very small alphabet (usually
one-byte letters, i.e. < 256 letters) and highly optimized for various assumptions.

Simileco is built for large, non-standard alphabets and, by using templates, easy yet efficient
customization at various points in the algorithm (e.g. lazy similarity matrices).

Currently supports:

* Needleman-Wunsch
* Smith-Waterman
* Waterman-Smith-Beyer

You might also be interested in alignment-tools, a suite of tools for analyzing and visualizing
alignment algorithms using Python and Mathematica.

## Usage

Build an aligner for a given maximum sequence size of sequences s and t. This 
will allocate all needed data structures internally once. You will then be able
to use this aligner for any pair of sequence below the given maximum lengths.

```
Aligner aligner(20, 20);
```

Now we define two sequences to align. Note that for simplicity, we use strings here,
but `Aligner` can handle any other type of sequence, such as, for example, vectors of
tuples of floats.

```
const std::string s("CHOCOLATEISTHEANSWER");
const std::string t("LATETHAW");
```

We will use Waterman-Smith-Beyer here. For this, we need to specify a similary matrix
as well as a gap function. The similarity matrix is any class or function that can be
called with two indices and gives a similarity score for the sequence items at these
indices. We use the following function, that gives a score of 1 identical letters, and
-1 for non-matching letters in s and t:

```
auto similarity = [&s, &t] (int i, int j) {
	return s[i] == t[j] ? 1 : -1;
};
```

Now we need a gap function. This function computes the gap cost for a given length.
Note that this should return a positive value - in specifications of the algorithms
this value is usually given as a negative score, such as -1. Again, you can use any
lambda or class you like, for example this one:

```
auto gap_cost = [] (int n) {
	return std::pow(1.25, n);
};
```

Now run the algorithm:

```
aligner.waterman_smith_beyer(
	similarity,
	gap_cost,
	s.length(),
	t.length());
```

Finally print the result:

```
std::cout << aligner.pretty_printed(s, t);
```

which gives:

```
CHOCOLATEISTH--EANSWER
     ||||  ||         
-----LATE--THAW-------
```

Penalizing gaps a bit less and changing our gap function to:

```
auto gap_cost = [] (int n) {
    return std::pow(0.5, n);
};
```

gives:

```
CHOCOLATEISTHEANSWER
     ||||  || |  |  
-----LATE--TH-A--W--
```

More examples can be found in the tests in `test.cpp`.

## Building

simileco needs the <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen library</a>, so you might need to add your custom Eigen headers as system includes. You also need to specify C++17 when compiling.

To build `test.cpp`, you also need the <a href="https://github.com/catchorg/Catch2">Catch framework</a>. Here's an example:

```
 clang++ test.cpp -isystem /usr/local/include/eigen3 -std=c++17
```

## License

MIT.
