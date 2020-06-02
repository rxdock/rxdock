Overview
========

|Dock| is a fast and versatile open-source docking program that can be used
against proteins and nucleic acids. It is designed for high-throughput virtual
screening (HTVS) campaigns and binding mode prediction studies.

The rDock program was developed from 1998 to 2006 (formerly known as
RiboDock [RiboDock2004]_) by the software team at RiboTargets (subsequently
Vernalis (R&D) Ltd.). In 2006, the software was licensed to the University of
York for maintenance and distribution under the name rDock. In 2012, Vernalis
and the University of York agreed to release the program as open-source
software. The released version is licensed under GNU LPGL version 3.0 with
support from the University of Barcelona --
`rdock.sourceforge.net <http://rdock.sourceforge.net/>`__. The development of
the open-source code stopped in 2014, so in 2019 RxTx decided to revive it
by forking rDock as RxDock.

The major components of the platform now include fast intermolecular scoring
functions (van der Waals, polar, desolvation) validated against protein and RNA
targets, a Genetic Algorithm-based stochastic search engine, a wide variety of
external structure-based drug discovery (SBDD) derived restraint terms (tethered
template, pharmacophore, noe distance restraints), and novel Genetic
Programming-based post-docking filtering. A variety of scripts are provided to
perform automated validation experiments and to launch virtual screening
campaigns.

This introductory guide is aimed at new users of |Dock|. It describes the minimal
set of steps required to build |Dock| from the source code distribution, and to
run one of the automated validation experiments provided in the test suite
distribution. The instructions assume that you are comfortable with simple Linux
command line administration tasks, and with building Linux application from
makefiles. Once you are familiar with these steps you should proceed to the User
and Reference Guide for more detailed documentation on the usage of |Dock|.
