# Chloe2.jl

Software for annotating plastid genomes. Follows on from Chloe.jl (https://github.com/ian-small/Chloe.jl), but is more sensitive and accurate (and slower)

## Dependencies

Chloe2 relies on two amazing software packages from Sean Eddy's lab, [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). Many programs from these packages must be accessible via your $PATH for Chloe2 to function. Follow the installation instructions on the [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/) websites to install these packages.
Chloe2 is written in [Julia](https://julialang.org); follow the [download and installation instructions](https://julialang.org/downloads/).

## Installation

For general use do the following

Tell Julia to treat the repo as a Julia package:

`julia`

type the character `]` which should bring up a prompt like: `(@v1.12) pkg>` 

Now type:

`add https://github.com/ian-small/Chloe2.jl`

and if you have a julia version > 1.11 you can also create a commandline script

`app add https://github.com/ian-small/Chloe2.jl`

(This requires you to have `~/.julia/bin` in your `PATH`)

You can now quit the Julia REPL with Ctrl-D

*Note*: Currently `add` and `app add` seem to be separate; you can do one or the other or both as you
wish. But -- say -- just adding the `app` *doesn't* give you access to the package e.g: `julia -m Chloe2` will fail.

## Running Chloe2

To annotate a single genome in .gff format

`julia -m Chloe2 --loglevel info --gff my_genome.gff my_genome.fasta`

*OR* if you have added the app (*and* have `~/.julia/bin` in your `PATH`) simply type:

`chloe2 --loglevel info --gff my_genome.gff my_genome.fasta`

To annotate many fasta files concurrently in the directory 'my_genomes' and save the generated .gff files in the same directory (note: Use consistant inputs/outputs. If you wish to annotate a directory of fasta files, ensure that the output options are also directories):

`julia -t 16 -m Chloe2 --loglevel error --gff my_genomes my_genomes`

In this case `-t 16` indicates julia should use 16 threads concurrently.

Running multiple threads concurrently will greatly increase speed until limited by number of cores, RAM or I/O. Each thread will consume ~1.2 GB of RAM.

For the command line use the bare `--` option to separate julia options from chloe2 options:  `chloe2 -t 16 -- --loglevel error --gff my_genomes my_genomes`

For more options, see

`julia -m Chloe2 --help`                                             

--edits file/dir for gff input containing edit site information

--pseudo reports incomplete or otherwise problematic features as pseudogenes

--max uses the --max setting for searches with HMMER and Infernal. Not recommended for general use.
