# bio-pileup_iterator

[Pileup format](http://samtools.sourceforge.net/pileup.shtml) files are a representation of an alignment/mapping of reads to a reference. This biogem builds on the [bio-samtools biogem](https://github.com/helios/bioruby-samtools) to create enable developers to iterate through columns of a pileup format file, and interrogate possible polymorphisms for e.g. SNP detection. Say we have the pile lines like so

    contig00001	199	A	4	.$...$	>a^>
    contig00001	200	T	2	..+1A	aR

i.e.

    line = "contig00001\t199\tA\t4\t.$...$\t>a^>\ncontig00001\t200\tT\t2\t..+1A\taR"

Then

    piles = Bio::DB::PileupIterator.new(line).to_a
    piles[0].reads #=> An array of 4 pileup reads (Bio::DB::PileupIterator::PileupRead objects)

The first reads ends at the first position

    piles[0].reads[0].sequence #=> 'A'

The second read covers both positions:

    piles[0].reads[1].sequence #=> 'AT'

Note that when you don't use ```to_a```, instead using ```Bio::DB::PileupIterator#each```, there is no "lookahead" (yet), so it doesn't find the T before it has iterated over it:

    Bio::DB::PileupIterator.new(line).each{|pile| puts pile.reads[1].sequence if pile.pos==199} #=> "A"

Directions

    piles[0].reads[1].direction #=> '+'

Insertions

    piles[1].reads[1].insertions #=> {200=>"A"}

Apologies in advance for any missing features (e.g. currently it does handle deletions) and slowness (it wasn't really written with speed in mind).

## Installation

        gem install bio-pileup_iterator

## Developers

To use the library 

        require 'bio-pileup_iterator'

The API doc is online. For more code examples see also the test files in
the source tree.
        
## Project home page

Information on the source tree, documentation, issues and how to contribute, see

  http://github.com/wwood/bioruby-pileup_iterator

## Cite

  If you use this software, please cite http://dx.doi.org/10.1093/bioinformatics/btq475

## Biogems.info

This Biogem is published at http://biogems.info/index.html#bio-pileup_iterator

## Copyright

Copyright (c) 2012 Ben J. Woodcroft. See LICENSE.txt for further details.

