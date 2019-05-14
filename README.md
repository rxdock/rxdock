# rDock: a fast, versatile, and open-source program for docking ligands to proteins and nucleic acids

![rDock Logo](docs/_static/logo.png)

rDock is a fast and versatile **open-source docking program** that can be used
to dock **small molecules** against **proteins** and **nucleic acids**. It is
designed for high-throughput virtual screening (HTVS) campaigns and binding mode
prediction studies.

rDock is mainly written in C++ and accessory scripts and programs are written in
C++, Perl or Python languages.

The full rDock software package requires **less than 50 MB** of hard disk space
and it is compilable (at this moment, **only**) in **all Linux computers**.

Thanks to its design and implementation, it can be installed on a computation
cluster and deployed on an **unlimited number of CPUs**, allowing HTVS campaigns
to be carried out in a **matter of days**.

Besides its main Docking program, the rDock software package also provides a set
of tools and scripts to facilitate **preparation** of the input files and
**post-processing** and **analysis** of results.

## About

### Features

![The above image illustrates the first binding mode solution for ASTEX system
1hwi, with an RMSD of 0.88 Å.](docs/_images/dock1.jpg)

#### Docking preparation

Define cavities using **known binders** or with user-supplied **3D
coordinates**. Allow -OH and -NH2 receptor side chains to rotate. Add explicit
solvent molecules and structural waters. Supply pharmacophoric restraints as a
bias to **guide docking**.

#### Pre-processing of input files

Define common ligand structure for performing **tethered docking** (requires
Open Babel Python bindings). Sort, filter or split ligand files for facilitating
**parallelization**. Find **HTVS protocol** for optimizing calculation time.
Pre-calculate grids to decrease subsequent calculation times.

#### Post-processing and analysis of results

Summarize results in a tabular format. Sort, filter, merge or split results
files. Calculate **RMSD** with a reference structure taking into account
internal symmetries (requires Open Babel Python bindings).

#### Binding mode prediction

Predict how a ligand will bind to a given molecule. The ASTEX non-redundant test
set for proteins and DOCK and rDock test sets for RNA have been used for
validating and comparing rDock with other programs.

#### High-throughput virtual screening

Run for million of compounds in short time by exploiting the capabilities of
computer calculation farms. Ease of **parallelization** in relatively unlimited
CPUs to optimize HTVS running times. The DUD set has been used for validating
rDock and comparing its performance to other reference docking programs.

![In red mesh, definition of the cavity obtained by execution of `rbcavity`
program.](docs/_images/dock2.jpg)

### History

The rDock program was developed from 1998 to 2006 by the software team at
RiboTargets (subsequently [Vernalis (R&D) Ltd](https://www.vernalis.com/)).
In 2006, the software was licensed to the
[University of York](http://www.ysbl.york.ac.uk/) for maintenance and
distribution.

In 2012, Vernalis and the University of York agreed to release the program as
open-source software. This version is developed with support from the
[University of Barcelona](http://www.ub.edu/cbdd/) --
[sourceforge.net/projects/rdock](http://sourceforge.net/projects/rdock/).

### License

rDock is licensed under [GNU LGPL version 3.0](LICENSE.md).

### Contributor concordat

rDock adheres both to No Code of Conduct and Code of Merit.

* [Contributor Code of Conduct](CODE_OF_CONDUCT.md)
* [Code of Merit](CODE_OF_MERIT.md)

### References

If you are using rDock in your research, please cite [rDock2014]. Former
software reference provided for completeness is [RiboDock2004].

* [rDock2014] Ruiz-Carmona, S., Alvarez-Garcia, D., Foloppe, N.,
  Garmendia-Doval, A. B., Juhos S., et al. (2014) rDock: A Fast, Versatile and
  Open Source Program for Docking Ligands to Proteins and Nucleic Acids. PLoS
  Comput Biol 10(4): e1003571.
  [doi:10.1371/journal.pcbi.1003571](https://doi.org/10.1371/journal.pcbi.1003571)

* [RiboDock2004] Morley, S. D. and Afshar, M. (2004) Validation of an empirical
  RNA-ligand scoring function for fast flexible docking using RiboDock®. J
  Comput Aided Mol Des, 18: 189--208.
  [doi:10.1023/B:JCAM.0000035199.48747.1e](https://doi.org/10.1023/B:JCAM.0000035199.48747.1e)
