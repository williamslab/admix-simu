Admixture simulator
===================
Program to simulate admixture between multiple populations using arbitrary
proportions of ancestry from any population in each generation. The model
uses a recombination map and discrete generations to model the location of
switches between haplotypes.

Usage:

    ./simu-mix.pl [in.dat] [in.snp] [out_prefix] [POPULATION FILES]

An example using data for 80 SNPs:

    ./simu-mix.pl example/AA.dat example/yri.snp AA -CEU example/ceu.phgeno -YRI example/yri.phgeno

The above call will create files called AA.phgeno and AA.bp in the current
directory. This call requires that the C++ program 'mixer' is compiled. See
"Compiling" below.

Information on each of the command line arguments as well as output files
resulting from simulation is given below.

------------------------------------------------------

Two modes of simulation
-----------------------

At every recombination breakpoint---including recombinations between haplotypes
from the same population---the script ./simu-mix.pl randomly selects a
different haplotype than the one being copied from to start copying from.
To simulate N haplotypes, at most N+1 haplotypes from each source population
will be needed. However, because ./simu-mix.pl is randomized, some runs will
succeed in simulating N haplotypes even when fewer than N+1 haplotypes are
available for each population.

An alternative script, ./simu-mix-2n.pl, requires 2N haplotypes from each
population in order to simulate N admixed haplotypes. This script assigns
2 haplotypes from each population to be used for simulating each haplotype.
When a recombination occurs between haplotypes in the _same_ population,
copying switches from one of these two haplotypes to the other.

------------------------------------------------------

Simulating only breakpoints
---------------------------

If haplotypes for the simulated samples are not needed, the mixer program
can be used to simulate only breakpoints according to a specified genetic
map and marker set.

Usage:

    ./mixer [in.dat] [in.snp] [out.bp]

See below for information on the format of the breakpoints file.

------------------------------------------------------

Compiling
---------

To compile the 'mixer' C++ program, most Linux/Unix-based users should simply
be able to type

    make

Other setups may require editing of the Makefile or alternate means of
compiling.

------------------------------------------------------

Dat file
--------

The dat file describes the admixture process to be simulated. Example dat files
are available in example/aa.dat and example/continuous.dat. The first line of
the file contains K+2 columns, where K is the number of source populations.
The first field on line one specifies the number of samples to be simulated.
The second field is a label for the admixed population that is produced during
the run. (This population can contribute have non-zero ancestry contribution
after the first generation in which admixture occurs; i.e., beginning on the
second line of the dat file this column can be non-zero). Columns 3 through
K+2 specifies labels for the source populations. These same labels, which are
case-sensitive, should be used in the option names given to the ./simu-mix\*.pl
script.

The second and all subsequent lines have the format:

    Generation#  p_0   p_1  ...  p_K

Here, p_i is the proportion of ancestry for the corresponding population
(beginning with the admixed population at index 0). The sum of all p_i
values on each line must be 1.0.

The first such line describes the first admixture event, where p_0 = 0 is
required. The admixture event occurs at generation 1 and is a pulse admixture
event. After the pulse admixture event, the population mixes with itself up to
the specified generation number. The next line in the file is then processed,
taking the number of generations from the previous line as the starting point
for further admixture events (i.e., each line must have a generation number
greater than the generation number on the previous line).

The following from example/AA.dat:

    40   Admixed   CEU    YRI
    6    0         0.2    0.8

is equivalent to:

    40   Admixed   CEU    YRI
    1    0         0.2    0.8
    2    1         0      0
    3    1         0      0
    4    1         0      0
    5    1         0      0
    6    1         0      0

Both specify a pulse admixture event between CEU and YRI (20% and 80% ancestry
from these populations, respectively) that occurred 6 generations ago with
no further ancestry from these populations after that event, and both will
generate 40 admixed haplotypes.

The file example/continuous.dat will simulate 40 haplotypes in which admixture
first occurred between the CEU and YRI population, with 5% and 95% ancestry
respectively. In generation 2 the admixed individuals contribute 20%
ancestry, CEU contributes 5% ancestry, and YRI contributes 75% ancestry.
These same proportions of ancestry occur in generations 3 through 6. Note
that the final population uses only 20% ancestry from the admixture that
occurred in prior generations and has at least 75% YRI ancestry. Note also
that continuous mixture such as that from this example must have each line
specified since, if the number of generations specified on a line is
more than 1 greater than the previous line, all generations after the first
are modeled as mixture of the population with itself (as in the more lengthy
example file given above compared to example/aa.dat).

Finally, note that a single generation of mixture will not result in distinct
local ancestry tracts on a single haplotype. Recombination in unadmixed
individuals produces chromosomes with a single ancestry composition (although
recombinations from the population to itself will be represented).

------------------------------------------------------

SNP file
--------

Eigenstrat format SNP file that specifies the genetic map to be used for
sampling recombination events. The insert-map.pl script can be used to
add a genetic map to a correctly formatted SNP file.

------------------------------------------------------

Output haplotype file (out_prefix.phgeno)
-----------------------------------------

Simulated haplotypes file in Eigenstrat format. The program will overwrite
any existing file with the name specified on the command line.

------------------------------------------------------

POPULATION FILES
----------------

This is a list of arguments, one for each source population specified in the
dat file. The format is -POP\_LABEL [pop.phgeno]. The option name (POP\_LABEL)
must be identical to the one given in the dat file, and it is case sensitive.

------------------------------------------------------

Output breakpoints file (out_prefix.bp)
---------------------------------------

Each run of the 'mixer' program (invoked by the simu-mix\*.pl scripts)
generates a breakpoint file. This file contains information on ancestry and
breakpoint locations for each sample. The first line contains the population
labels (from the dat file) that were used for admixing. The remainder of the
file contains one line per sample with ancestry and breakpoints specified as
\[pop\_index\]:\[marker\_index\]. Here \[pop\_index\] is a 0 indexed population
number with population 0 being the first population listed on line 1
population 1 the second population, etc.  (Note that the column for the
proportions of admixed individuals in the dat file does not count as a
population.) The \[marker\_index\] is 0 based and is inclusive.

Example bp line:

    0:20 1:30

The above individual has ancestry from population 0 at markers 0 through 20,
inclusive (i.e., 21 markers) and ancestry from population 1 at markers 21
through 30 (10 markers).
