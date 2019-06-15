[![GitHub release](https://img.shields.io/github/release/ViennaRNA/Treekin.svg)](https://www.tbi.univie.ac.at/RNA/Treekin/#download)
[![Build Status](https://travis-ci.org/ViennaRNA/Treekin.svg?branch=master)](https://travis-ci.org/ViennaRNA/Treekin)
[![Github All Releases](https://img.shields.io/github/downloads/ViennaRNA/Treekin/total.svg)](https://github.com/ViennaRNA/Treekin/releases)
[![Conda](https://img.shields.io/conda/v/bioconda/treekin.svg)](https://anaconda.org/bioconda/treekin)
[![AUR](https://img.shields.io/aur/version/treekin.svg)](https://aur.archlinux.org/packages/treekin/)

# Treekin

Compute folding dynamics on coarse grained energy landscapes through numeric integration of the underlying Markov process

----

## Description

The program `treekin` computes folding dynamics on coarse grained version of an energy landscape, where all
conformations belonging to the same local minimum have been contracted into a single macro-state. `treekin`
reads in the a list of local minima (macro-states) and effective refolding rates between these minima as
computed by `barriers --rates`. Since the the number of macro-states is small (typically on the order of
100-1000), folding dynamics can be computed by direct numerical integration of the master equation
conformations (i.e. diagonalization of the rate matrix).

For detailed instructions see the man page.

----

### ATTENTION:

use `--rrecover` and `--wrecover` options with caution:

first: run `treekin` with `-w`, i.e. let it diagonalize the input matrix once and write
       the eigenvalues and eigenvectors to the corresponding *.{evals,evecs}.bin files

second: run `treekin` with the `-r` option, i.e let it read those *.{evals,evecs}.bin files
        and do ONLY the iteration (no diagonalization)

NOTE: when using `--rrecover` and `--wrecover` option in combination with the `-a` option,
      `treekin` assumes that the absorbing state between consecutive calls of the program 
      are the _same_ (since the information about the absorbing state influences the 
      transition matrix).
      If you change the `-a` option while calling `treekin` with `--rrecover`, `treekin` will
      produce junk output (since the transition matrix that has been diagonalized earlier 
      is different from the current one, due to the different absorbing state)

----

## References

If you use this software, you may want to cite the follwing publications:

- [M.T. Wolfinger et al.](https://doi.org/10.1088/0305-4470/37/17/005),
  "Efficient computation of RNA folding dynamics",
  J. Phys. A: Math. Gen. 37 4731

----

## License

Please read the copyright notice in the file [COPYING](COPYING)!

If you want to include this software in a commercial product, please contact 
the authors. 
