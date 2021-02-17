<center> <h1>PangenomeVis</h1> </center>

# GOAL: An interactive, JavaScript-based browser for microbial pan-genome graphs 

Variant discovery is largely conducted by comparing sample genomes to a single canonical reference genome under the assumption that the sample of interest differs only slightly from its reference.  However, many pathogens diversify their genomes so quickly that the sample and reference may differ greatly.  When these differences are sufficiently large, existing algorithms fail to adequately align sample sequence data to the reference, and variation underlying virulence factors (e.g. drug resistance, immune evasion) may go undetected.

The Broad Institute is developing a new variant discovery algorithm based on the emerging “graph genome” paradigm.  A graph genome is a data structure (formally a multigraph of vertices, edges, and sample metadata) that contains many reference genomes as well as sequencing data for the sample(s) of interest.  Subsequences from one or more genomes are represented by vertices, while adjacencies between subsequences are represented by edges.  Each edge is assigned a color that is associated with a specific sample, thus allowing for the data structure to represent a large collection of genomes as a single compact and navigable object.  Theoretically, any type of sequence variant could be discovered by simple graph navigation operations.  For example, a trio of samples (mother, father, and child) may be assigned the colors blue, green, and red (respectively).  A de novo variant in the child is simply a bifurcation of the child’s genome from both parents (depicted in the figure to the right).  Variant discovery in this context is essentially walking the parental and child colors independently from where they diverge to where they coalesce.

Navigating these data structures is conceptually easy, but surprisingly difficult in practice as data may have sequencing errors, gaps, tangles, and other obstacles, making visualizing genome graphs an important step for developing variant discovery algorithms. Thus, we’re seeking collaborators to help develop an interactive web tool to display user-selected portions of a multi-sample genome graph (“subgraphs”).

# For users

An example of a site displayed

<img src="misc/example.png" alt="rGFA example" height="350"/>

MORE TO BE FILLED IN...

# For developers

If interested in contributing, please create tickets, or better yet, join the Slack channel [#dsp-pangenome-vis](https://broadinstitute.slack.com/archives/C018JAW5YF9).

## Technology used

We currently use JavaScript and the React + ReactDOM framework, coupled with [D3](https://d3js.org/). For backend (for which there isn't much support yet), we use [FastAPI](https://fastapi.tiangolo.com/).




