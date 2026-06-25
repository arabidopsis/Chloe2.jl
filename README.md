# Chloe2.jl
Software for annotating plastid genomes. Follows on from Chloe.jl (https://github.com/ian-small/Chloe.jl), but is more sensitive and accurate (and slower)

## Dependencies
Chloe2 relies on two amazing software packages from Sean Eddy's lab, [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). Many programs from these packages must be accessible via your $PATH for Chloe2 to function. Follow the installation instructions on the [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/) websites to install these packages.
Chloe2 is written in [Julia](https://julialang.org); follow the [download and installation instructions](https://julialang.org/downloads/).

## Installation
Clone this github repo:

`git clone https://github.com/ian-small/Chloe2.jl.git`

Tell Julia to treat the repo as a Julia package:

`julia`

`]`

`dev '~/github/Chloe2.jl'`

`instantiate`

(replace ~/github with whatever path you've cloned the Chloe2.jl repo to)

This will install all the Julia packages that Chloe2 needs, and the models it needs to annotate plastid genomes.

You can now quit the Julia REPL with Ctrl-D

## Running Chloe2
`julia --project=~/github/Chloe2.jl ~/github/Chloe2.jl/src/command.jl --help`                                             
Usage: Chloe2.jl/src/command.jl [options] <FASTA file or directory>

Optional Arguments:

--gff file/dir for .gff output. Required if you want annotations.

--loglevel level of information provided during operation; one of (info,warn,error,debug). Would recommend 'info' if running single-threaded, 'error' if multi-threaded

--edits file/dir for gff input containing edit site information

--pseudo reports incomplete or otherwise problematic features as pseudogenes

--max uses the --max setting for searches with HMMER and Infernal. Not recommended for general use.

Note: Use consistant inputs/outputs. If you wish to annotate a directory of fasta files, ensure that the output options are also directories.

## Examples

To annotate a single genome in .gff format:

`julia --project=~/github/Chloe2.jl ~/github/Chloe2.jl/src/command.jl --gff my_genome.gff my_genome.fasta`

To annotate many fasta files concurrently in the directory 'my_genomes' and save the generated .gff files in the same directory:

`julia --project=~/github/Chloe2.jl -t 16 -- ~/github/Chloe2.jl/src/command.jl --gff my_genomes my_genomes`

In this case `-t 16` indicates julia should use 16 threads concurrently. (Note the `--` to separate julia options from Chloe2's)

Running multiple threads concurrently will greatly increase speed until limited by number of cores, RAM or I/O. Each thread will consume ~1.2 GB of RAM.

## Install Chloe2 as a package

You can also invoke Chloe2 if you have it installed in your project (e.g. with say
`] add https://github.com/ian-small/Chloe2.jl.git`).

`julia --project=. -e 'using Chloe2; main()' -- --gff my_genomes my_genomes`

