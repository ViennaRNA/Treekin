[AUTHOR]

Michael T. Wolfinger, Marcel Kucharik, Ivo Hofacker, Christoph Flamm, Andreas Svrcek\-Sailer, Peter Stadler.

[EXAMPLES]

Typically, computation of a reduced dynamics based on the analysis of folding landscapes
requires two steps: Elucidation of the landscape (topology) and \- based on that \- 
calculation of the reduced dynamics.


The first step involves computing the relevant properties of an energy landscape by
\fBbarriers\fR (see \fBbarriers\fR(1) for details). The resulting .bar\-file contains
information on local minima, basins, saddle points as well as thermodynamic properties
of the energy landscape. Additionally, the \-\-rates option in the below example triggers
\fBbarriers\fR to generate another output file (rates.out) containing the transition
rates between all pairs of macrostates (ie. basins of attraction), calculated by summing
over the corresponding microscopic rates. 


.nf
.ft CW
  $ barriers \-\-saddle \-\-bsize \-\-rates < foo.sub > foo.bar
.ft
.fi


In a second step, \fBtreekin\fR is called with options to calculate the macrostate
dynamics on the previously generated landscape by means of applying microscopic
transition rates (option \-m I):


.nf
.ft CW
  $ treekin \-\-p0 2=1 < rates.out
.ft
.fi


In this example, the simulation starts with 100% of the initial population in macrostate
2, i.e. the second lowest minimum in the barrier tree (option \-\-p0 2=1).  The
transition matrix is computed from a set of microscopic rates, read from a rates file
(as computed by \fBbarriers\fR).


Generally, calculation of the macrostate dynamics by means of microscopic rates
(option \-m I) is consiberably more accurate than the simplified Arrhenius\-like
dynamics (option \-m A).


Looking at the default output produced by \fBtreekin\fR, there are two sections:
Overall status information on the computation (marked by hash signs at the beginning
of the line) are printed at the top. Below, the actual data is printed for each time
step in (n+1) space\-separated columns, where n is the number of investigated
(macro)states. The first column lists the current time, whereas all remaining columns
correspond to the population probabilities of individual (macro)states.

[>REFERENCES]

M.T. Wolfinger, W.A. Svrcek\-Seiler, Ch. Flamm, I.L. Hofacker, P.F. Stadler
"Efficient computation of RNA folding dynamics"
J.Phys.A: Math.Gen. 37: 4731\-4741 (2004)

I.L. Hofacker, Ch. Flamm, Ch. Heine, M.T. Wolfinger, G. Scheuermann, P.F. Stadler
"BarMap: RNA folding on dynamic energy landscapes"
RNA: 2010 16: 1308\-1316 (2010)

[SEE ALSO]

barriers(1)
