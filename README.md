
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
* Separate IO, API-specific, and coordinating code, like oil & water.


Inforfmation comes from failed tests; too much test code == lower velocity.
Quck-check

