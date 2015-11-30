
##Phylogenetic Trees

The smallest elements in a distance matrix are not necessarily two neighbors within a tree.


##Neighbor-joining
Takes an *additive* matrix. Each time will pick the closest pair *farthest* from other leaves.


Replace minimum element with new leaf `m`, refilling the matrix using the Additive Phylogeny Equation (used to find insert place of parent). Do this until we have `n-1` leaves (all neighbors (rows) replaced with leaves(?)). Construct a subtgree using this forumla.

Add 2 limbs back to each node using limb forumla (Recopmute D*, keeping and using new D)

Additive Matrix: The distances are equal to the sum aalong the path. May not work if there are unexplained mutations.

Additive Phylogeny: Using a distance matrix, construct new nodes by building internal nodes between random leaves. i.e. attach a common internal/center node.




##K-means clustering
0. (don't repeat) Pick k centers as k random points.
1. Assign all points to clusters by whatever is the closest centroid.
2. Set new centers by finding the middle of each set of cluster points.
3. Repeat unntil center positions stop changing.


#Misc. Software development

Behavior/biologist-driven Design 
* Robot Framework
* Test by compairing a series of files

###Testing
Reasons for and types of testing:
1. Regression (system-level testing)
2. Establish correctness
3. Explain functionality
4. Find bugs
5. Assist design

Dos:
* __Use assertions__
* external tests of correctness ("ground truth")
* Unit test algorithms,
* Integration test "coordinating"/glue code.
* Create unit tests for any found bugs.
* Separate IO, API-specific, and coordinating code, like oil & water.


Inforfmation comes from failed tests; too much test code == lower velocity.
Quck-check

Travis only gets one core.
SimSeq to establish correctness
###Other

bio_pieces:
version .x.y.z  x-> don't change, y-> change when user interaction changes or add a scirpt, z -> bug fix or new parameters.
dev is z=0.
minimize maintenance




### Lab Chemistry
PCR: know how much template (real DNA) there is an associated primer file
primer errors happen at the ends of reads (so SNPs near ends of reads are likely errors)


SNP is stop codon in middle of coding region; indel is not a multilpe of 3
is SNP in other sample/references?
The third nt can change in codon and AA won't change.

invalid/dead particles v. errrors
PIs will alter consensus/VCF itself after reviewing in geneious etc.













Bayesian + ML(maximum likelihood) allow for applying evolutionary model to tree construction;
  model: includes info like transition & transversion rates (ie.e. A->G v. A-T, etc.)
  some areas have high/low variability (e.g. antigen v. internal/strucgtural)
  Dengue is a polypeptide protein. so no stop codons until end.

  1s. parameter to model: rate of mutatiom (transitions & transversions). calculated from dataset.
  3rd param: rate of actual NTs (A,C,T,G)
  Gamma distribution to model 3rd nt on changin in codon (wobble posiition)
  model AA idff. from NT subs
  invairant site can bials alignment/tree (e.g. protein demands "T" at position i)

Use likelihood ratio test to pick model.
Don't want to overparameterize (deg. of freedom)
* FigTree

expect ambiguous bases. 

SOP for phylo tree construction: 
  aligned fasta (.aln)
  want all sequences the same length
  diff. datasets have diff. gene lengths (otherwise will move longer sequnces further away in tree).


Dengue variation gets drowned out by the consensus. 







