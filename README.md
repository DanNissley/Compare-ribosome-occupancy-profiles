# Compare-ribosome-occupancy-profiles
Code, inputs, and outputs from comparison of ribosome occupancy profiles between yeast domains. This code tests the hypothesis that translation rate profiles are conserved between protein domains that are structurally and evolutionarily related to one another. 

# Contents of repository
A_site_profiles     : a directory containing the (unnormalized) A-site read profiles for the Nissley, Jan, Weinberg, Williams, and Young data sets. pooled.tab is the sum of the other files at every nucleotide position

inpfiles.zip        : compressed directory containing reference genome files for S288C_reference_genome_R64-2-1_20150113, SUPERFAM domain assignments for S288C, and SCOP2.07 for defining structural classes of domains

pooled_k0_30-80.bck : sample output that results from running the commands in "commands.txt"

src                 : contains the python program compare.py that is loaded by align.py and runner.py to perform comparisons; loading it as a module gives you access to various methods for parsing domain information and performing comparisons in different ways

align.py            : python program that runs analysis up to and including AA and DNA alignments between domains but does not run comparisons (see below for reason)

runner.py           : python program that runs complete analysis excluding alignments to streamline processing

commands.txt        : commands that can be run to generate the data in pooled_k0_30-80/

# Running comparisons between ribosome occupancy profiles

## Download the repository and prepare files

Extract inpfiles.zip to provide access to its contents

Open align.py and runner.py and change the line reading "sys.path.append('/data/icarus/nissley/trans_sched/fullPipeline/update_v4/trans_sched/src')
" to instead direct sys to append the absolute path to the src/ directory

## Perform initial MUSCLE alignments

Use the first command in commands.txt to run align.py. This will generate the directory example/ and run analysis up to and including amino acid and DNA alignments between pairs of domains in the same SUPERFAM family. If you want to run many independent comparisons to different selections of unrelated domains this prevents having to run alignments for each iteration, saving time. That is - just generate the alignments once so your independent comparisons can all use them.

## Compare ribosome occupancy profiles

Use one or all of the other 20 commands in commands.txt to call runner.py and run comparisons between ribosome occupancy profiles for related domains. For the pooled data set this can be time consuming (~12 hours) due to the large number of random domains that must be selected and processed. 

Output files will be written to examples/1_pairwise.txt, examples/2_pairwise.txt, ... , examples/20_pairwise.txt

Each chunk of 20 lines in these files is a comparison between one pair of related domains and between the first related domain and 19 unrelated domains. First column is the integer index of the current pair of related domains. Second column is a second integer, 0 through 19, that keeps track of the 20 comparisons run for each pair of related domains. Column 3 notes whether or not this line corresponds to a pair of related domains or not. Column 4 gives the pairs of domains compared in that line. Column 5 is the fraction of positions within the aligned ribosome occupancy profiles with the same classification as fast, medium, or slow translating (see manuscript for details). 

runner.py will also parse the output in each of these pairwise comparison files to generate a summary file, e.g. pooled_k0_30-80.bck/1_slowMidFast.txt. These files give the number of pairs of related domains placed in ranks 1 through 20 (see manuscript Methods). 

Finally, a summary file can be generated that gives the average number of domain pairs in each rank over the 20 trials (e.g. pooled_k0_30-80.bck/slowMidFast_summary.txt). These data are plotting in Figure 2b.

## Other information

Read function headers in src/compare.py for detailed information on inputs

Many functions in compare.py have the "verbose" argument, which is Boolean False by default. Setting it to Boolean True will provide debugging information. Note that it will print a very large amount of data to the screen and, if piping output to file, quickly generate large files. 

The program src/verbose.py is identical to src/compare.py except it has some hardcoded print statements that provide a massive amount of information on profiles being generated and processed. In combination with setting "verbose=True" for various functions this program is most useful for debugging. Note that you will need to create several folders, starting with new_output/, to which profiles will be written. Search for "new_output" in verbose.py or compare.py for details.
