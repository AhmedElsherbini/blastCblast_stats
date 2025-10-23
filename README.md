# blastCblast_stats
**Kindly if you find this repo useful for your work, cite & star this repo**

**What is this script?**

[Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) or [Cblaster](https://github.com/gamcil/cblaster) are great tools for finding homologous and co-located homologous respectively. 
Here, you can break down the long, complicated results  into a simple result per species, which will count the number of hits per species.
Then, it will efetch the number of the assembled genomes per species (via Biopython [Enterez](https://biopython.org/docs/1.75/api/Bio.Entrez.html) based on the NCBI assembly database). Dividing both numbers will hint at the spread of your cluster among different species.
Finally, this script will draw the smallest tree (based on ete3 definition) based on the pre-defined NCBI taxonomy among your species using [ete3 toolkit](http://etetoolkit.org/docs/latest/tutorial/index.html). Finally, merging the results of the database as a pie chart with the tree will give you a nice visualization.

**What do you need?**
>>>For NCBI BlastN (DNA) users, you need the description CSV table as explained here,

![alt text](https://github.com/AhmedElsherbini/blastCblast_stats/blob/main/Ncbi_blast.png)

>>>For NCBI BlastP (protein) users, you need the description CSV table as explained here,

![alt text](https://github.com/AhmedElsherbini/blastCblast_stats/blob/main/blastp_update.png)


>>>For Cblaster users,
You shall have the binary file as easily as I get it, like this.

```bash
cblaster search --query_file CP018841.1.faa --binary example_binary.csv -bde "," -bhh -bdc 6 -mi 50 -mc 50 -hs 3000
```
PS: <code>hs</code> is handy if you have a lot of results due to low coverage <code>mc</code>, low identity search <code>mi</code>. suggested by the last author in this [issue](https://github.com/gamcil/cblaster/issues/96).

So, type this command effortlessly.

```bash
 python blastcblast_stats.py -i  example_blastn_binary.csv -m strict -og deinococcus_radiodurans
```
```bash
 python blastcblast_stats.py -i  example_blastp_binary.csv -m strict -og deinococcus_radiodurans
```

```bash
 python blastcblast_stats.py -i  example_cblaster_binary.csv -m strict -og deinococcus_radiodurans
```
"-i /--input_dir"  is your path to the directory for your <**BlastN, P ,or Cblaster binary file**>  

"-og /--outgroup" <**optional**> is an outgroup species that I know it is NOT in my results and phylogenetically far from my results. 
*PS: do not forget to use underscore _ in the name of this species.* 

 **What about dependencies?**

Pandas, Biopython, ete3, argparse

Well, for [ete3](http://etetoolkit.org/download/), I recommend installing it via conda env (even if it takes a lot of time) if pip does not work properly.


**What do you get?**

Currently, there are three files.

1. database_percentage_*your_file*.csv (The main output where you can find for each species the count of this species in the binary file, the number of assembled genomes per NCBI assembly database, and the percentage of (count/assembly)*100. 
2. *your_file*_tree.nwk. If you would like to take this tree to a visualization tool (iTOL,FigTree,..)
3. *your_file*_tree.pdf. This is just a basic tree that links your isolates together, but with a pie chart that shows the results of file number 1.
   ![alt text](https://github.com/AhmedElsherbini/Cblaster_stats/blob/main/example_binary_tree_with_pies-1.png)

I hope this helps.

Thanks
