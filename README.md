# uorfs
Software to calculate the effect of a variant in a 5-prime UTR on the upstream open reading frames.

See https://www.biorxiv.org/content/10.1101/543504v2

This software will take a description of a variant, and calculate the set of uORFs in the reference and alternate alleles, then calculate the change in consequence. The software works on SNVs and InDels.

Matthew Wakeling

## Running Java
This code requires the htsjdk library. The easiest way to get everything required is to download the GATK Jar. All operations require this GATK jar and the Uorf.java file to be in the Java classpath, in order for Java to find them.

This can be done in two ways. The first option is to set the CLASSPATH environment variable:
```
export CLASSPATH=/path/to/GenomeAnalysisTK.jar:/path/to/software
```
The ":" character separates the two parts of this path, to specify that code can be found in the two places. The second option is to specify the "-cp" option every time you run java, like this:
```
java -cp /path/to/GenomeAnalysisTK.jar:/path/to/software blah blah blah
```
For all subsequent code fragments, where "java" or "javac" is specified, it is assumed that the classpath is correctly configured as specified above, either by adding the "-cp" option or using the CLASSPATH environment variable.

If you have a large server, it is also sensible to add the -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 options to java, to prevent it creating too many garbage collection threads. This is a small performance enhancement, and is done like this:
```
java -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 blah blah blah
```

## Compiling
Compiling the code is then done by:
```
javac *.java
```
in the directory with Uorf.java.

## Usage
The software should ideally be integrated into a larger system. The calculateUorfEffect method performs the calculation and returns a UorfResult object back with the results. However, the software also has a command-line interface. This is:

```
java Uorf <genome.fasta> <chr> <position> <reference_allele> <alternate_allele> <gene_strand> <5'UTR_start> <5'UTR_end>
```

The genome.fasta is the genome reference fasta file, which must match the genome against which the variant was called.
The chr, position, reference_allele and alternate_allele are the description of the variant, and can be taken directly from a VCF file (they are the first four columns of the file). The gene_strand is 1 if the gene is read in the forward direction and -1 for the reverse direction. The position of the 5'UTR is then described with start and end positions (inclusive). If the UTR is spread over multiple exons, then the start and end should be specified for each one. For example:

```
java Uorf human_g1k_v37.fasta 19 633529 G GGCGCCGCCGCCGCCGCCGCC -1 633513 633568
```
describes an insertion in the 5'UTR of the POLRMT gene. The software produces the following output:

```
Effect: out-of-frame_oORF
Start codon strength: Weak
Start codon distance: 65
ORF finish distance: 0
uORFs in reference: [(EXTENDING, distance: 45 to 0, Weak)]
uORFs in alternate: [(FRAMESHIFT, distance: 65 to 0, Weak)]
Reference UTR bases, offset 0: gg cgc ggg cgc ATG1CGC AGG CGC GGG CCG GTG GGG TGG CCT GGA GCG GCG TGC GTA 
Alternate UTR bases, offset 0: g gcg cgg gcg cat gcg cag gcg cgg gcc ggt ggg gtg gcg gcg gcg gcg gcg gcg gcg cct gga gcg gcg tgc gta 
Reference UTR bases, offset 1: ggc gcg ggc gca tgc gca ggc gcg ggc cgg tgg ggt ggc ctg gag cgg cgt gcg ta
Alternate UTR bases, offset 1: gg cgc ggg cgc ATG1CGC AGG CGC GGG CCG GTG GGG TGG CGG CGG CGG CGG CGG CGG CGC CTG GAG CGG CGT GCG TA
Reference UTR bases, offset 2: g gcg cgg gcg cat gcg cag gcg cgg gcc ggt ggg gtg gcc tgg agc ggc gtg cgt a
Alternate UTR bases, offset 2: ggc gcg ggc gca tgc gca ggc gcg ggc cgg tgg ggt ggc ggc ggc ggc ggc ggc ggc gcc tgg agc ggc gtg cgt a
```
This shows that the variant creates an out-of-frame_oORF, which means that the ORF overlaps with the gene start codon in a different frame. The reference allele also has an overlapping ORF, but it is in frame with the gene, so it is classified as an "EXTENDING" ORF. While the variant does not change the ORF start codon, it does change the frame, so the ORF becomes a "FRAMESHIFT" ORF.

The original ORF has a start codon 45 bases from the gene start, but because there is an insertion the alternate allele has the ORF start 65 bases from the gene start.

The last six lines show a visualisation of the ORFs in the UTR. The reference allele and the alternate allele are shown with spaces grouping bases in codons, with the three different frames. For each, the detected ORFs are indicated by the bases being capitalised. The number 1 that takes the place of the space after the ATG at the beginning of the ORF is an indicator of the strength of the start codon, from 1 (weak) to 3 (strong).

Another example is taken from figure 3c in the paper linked above:
```
java Uorf human_g1k_v37.fasta 22 29999922 A AT 1 29999885 29999987
```
which produces:
```
Effect: out-of-frame_oORF
Start codon strength: Strong
Start codon distance: 98
ORF finish distance: 0
uORFs in reference: [(NON_OVERLAPPING, distance: 97 to 61, Strong)]
uORFs in alternate: [(FRAMESHIFT, distance: 98 to 0, Strong), (NON_OVERLAPPING, distance: 67 to 61, Medium)]
Reference UTR bases, offset 0: g cca cca tgg tgg ccc tga ggc ctg tgc agc aac tcc agg ggg gct aaa ggg ctc aga gtg cag gcc gtg ggg cgc gag ggt ccc ggg cct gag ccc cgc gcc 
Alternate UTR bases, offset 0: gc cac cat ggt ggc cct gag gcc tgt gca gca act cca tgg ggg gct aaa ggg ctc aga gtg cag gcc gtg ggg cgc gag ggt ccc ggg cct gag ccc cgc gcc 
Reference UTR bases, offset 1: gc cac cat ggt ggc cct gag gcc tgt gca gca act cca ggg ggg cta aag ggc tca gag tgc agg ccg tgg ggc gcg agg gtc ccg ggc ctg agc ccc gcg cc
Alternate UTR bases, offset 1: gcc acc ATG3GTG GCC CTG AGG CCT GTG CAG CAA CTC CAT GGG GGG CTA AAG GGC TCA GAG TGC AGG CCG TGG GGC GCG AGG GTC CCG GGC CTG AGC CCC GCG CC
Reference UTR bases, offset 2: gcc acc ATG3GTG GCC CTG AGG CCT GTG CAG CAA CTC CAG GGG GGC TAA agg gct cag agt gca ggc cgt ggg gcg cga ggg tcc cgg gcc tga gcc ccg cgc c
Alternate UTR bases, offset 2: g cca cca tgg tgg ccc tga ggc ctg tgc agc aac tcc ATG2GGG GGC TAA agg gct cag agt gca ggc cgt ggg gcg cga ggg tcc cgg gcc tga gcc ccg cgc c
```
(Note that the 5'UTR doesn't start from 22:29999885 in this gene, but using this value cuts off the left-hand side without changing the result.) As described in the paper, this variant causes a frameshift in the existing short non-overlapping uORF, which creates a new frameshift uORF that overlaps with the gene. It also creates a new start codon, which is in frame with the stop codon of the original uORF, creating a new shorter non-overlapping uORF.
