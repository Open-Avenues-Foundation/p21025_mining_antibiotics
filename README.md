# p21025_mining_antibiotics
Tutorials and resources for p21025 "Genome Mining for New Antibiotics"

# Initial setup

In research we're mostly using UNIX based systems for our daily work which includes Linux and MacOS. Windows added the Windows Subsystem for Linux ([WSL](https://docs.microsoft.com/en-us/windows/wsl/about)) to make it easier for Windows users to use Unix based programs as well. In this project we'll be using Visual Studio Code [VScode](https://code.visualstudio.com) since it's available for all platforms and supports [WSL integration](https://docs.microsoft.com/en-us/windows/wsl/tutorials/wsl-vscode)

From the section below follow your system specific instructions and afterwards proceed to the [general instructions](#general-instructions) section.

# Getting started with the terminal

We'll be using the terminal to download data, using a program called `curl` and interact with Python. Go to [ncbi genome](https://www.ncbi.nlm.nih.gov/genome/?term=Aspergillus+nidulans+FGSC+A4). Download the proteins of that organism using the link below (also reachable  by clicking on the **protein** link). `gunzip` extracts the *compressed file*.

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/205/GCF_000149205.2_ASM14920v2/GCF_000149205.2_ASM14920v2_protein.faa.gz
gunzip GCF_000149205.2_ASM14920v2_protein.faa.gz


```

## Fasta format

`.fasta` or `.faa` file are used to collect (mostly protein) sequences. They start with an *identifier line* marking the beginning with `>`. On the following lines you can see the protein sequence in *one letter code* until the next record starts (`>`).

```
>XP_647833.1 hypothetical protein AN0001.2 [Aspergillus nidulans FGSC A4]
MVITDRLTKGVILEGMSEIDSESVAWALVRVLISKHGIPKAITSDRGSQFTRVSPFYLSYRYNLSLFSPTKEVEQLAEEP
AKSPIQKGEAIIQKVKEALD
>XP_647834.1 hypothetical protein AN0002.2 [Aspergillus nidulans FGSC A4]
MGADILAKSKYICWILKAITSDRGSQFTRVSPFYLSYGYNLSLFTPTEEVEHLAEEPTKSPIQKGEAIVQKVKEALDWAQ
ASMAYSQQNTENQANKHRSPATNYQVGDKVWLSLKNIRTDQPSKKLDWKNTKYEVIGLVGSHAVQLNTPPGIHPVFHVDL
LWLASSDPLPSQKNNDSQPPSIMVNGEEEYMVEKILDKRRRRYGRGHRLEYLVKWSGYAQPTWEAATALEEVQALDEWLD
HLSDSAKACAVLKARFLLPNIQGFIQQRIQQAGDSQSSQTYIIAAIVHDLSRRIFLSKRAPGEYQMPSSSTCRHI
```

# Getting started with python

In Python, there are several data types or objects:

* `"Strings"` represent text enclosed by `"` or `'`.
* `1` are **integers** (whole numbers)
* `variables` are **not enclosed with** `"` or `'` and used to store e.g. **strings** or **integers**.
* `#` is only a comment and will not be executed (see below).

E.g. we would like to assign a dna sequence and store it in a variable called `dna_sequence` we would enter the following into the *terminal*.

```python
dna_sequence = 'GTGTCACTTTCGCTTTGGCAGCAGTGTCTTGCCCGATTGCAGGATGAG'
```

## Loops and lists

We can loop over an object, do some calculation, and append out result to a list to keep for later.

```python
al = []
for i in [1, 2, 3]:
    al.append(i+i)

print(al)
```

We can do the same with a **list comprehension**.

```python
al = [i+i for i in [1, 2, 3]]

print(al)
```

### Important objects

**Lists** represented as **square brackets** `[]` are used to store things in order. You can access a list by *index* or *append* elements. Copy paste line for line into your python session execute:

```python
storage = ["a", "b", "c"]  # store "a", "b", "c" in the list called storage
len(storage)  # Return the "length" of the list
storage[0]  # Return value of storage at index 0, i.e. "a"
storage[1]  # Return value of storage at index 1, i.e "b"
storage[-1]  # Return the last element
storage[1:]  # Return all elements starting from the "1st"
storage.append("d")  # Append "d" to storage
len(storage)
print(storage)  # print all values of storage
```

**Dictionaries** represent as **curly braces** `{}` a storage for *fast retrieval*. They consist of *key*, *value* pairs.

> The protein sequences were truncated for better visibility!

```python
proteins = {"AN0001.2": "MGADILAKSKYICWILKAITSDRGSQFTRVSPFYLSYGYNL",
            "AN0002.2": "MGADILAKSKYICWILKAITSDRGSQFTRVSPFYLSYGYNLSL"}

proteins["AN0001.2"]  # "MGADILAKSKYICWILKAITSDRGSQFTRVSPFYLSYGYNL"

# We can add more values using a different syntax
# We cannot just append (like with the list), we need a key and a value:

proteins["AN0003.2"] = "MGADILAKSNLSL"  # adding another key:value pair
```

In general: Use *lists* if you want to work on all elements of it and use *dictionaries* if you want to look up certain elements quickly by their identifier.

### Defining functions

```python
def add_1(x):
    results = x+1
    return result

print(add_1(1))
```

## Installing packages

*Packages*, such as [*BioPython*](https://biopython.org), can easily be installed in python. We will need to create a *virtual environment* to install specific packages and leave your main python distribution untouched. *Virtual environments* are neccesary to keep dependencies isolated over several projects*. In your terminal, try the following:

```bash
python3.8 -m venv genome_mining  # Create a virtual environment for the project
source /path/to/genome_mining/bin/activate  # Activate the virtual environment
pip install biopython
```

### Imports

`import <package_name>` enables us to use functions from packages. In the snippet below we are using the `os` package to list the contents of the current directory:

```python
import os
print(os.listdir)
```

### BioPython

We'll be using the BioPython package during the project. Extensive documentation can be found [here](https://biopython.org/wiki/Documentation). For the project, the most important functions consider the reading of biological sequence data.

```python
from Bio import SeqIO

records = []

for r in SeqIO.parse("Aspergillus_nidulans_FGC4.gbff",
                     "fasta"):
    records.append(r)

print(records)
```

# Why the hassle of programming in the biotech industry?

A common trend today is that companies require their employees to analyze experimental or genomics data in Python. Furthermore, companies become more integrated and build whole pipelines of python scripts that process the company's data.


# Determining genome quality

Our company set out on an expedition to find producers of new antibiotics and our scouting team just came back from collecting samples. After isolation of fungal strains from environmental samples, the molecular biology team is harvesting genomic DNA for [illumina sequencing](https://www.illumina.com/techniques/sequencing/dna-sequencing.html). Reads automatically trigger the assembly and annotation pipeline and you find a new genome, ready for analysis, in the company database.

* What are some general genome statistics you'd like to collect?
  * Considering our target (antibiotic genes), what's an ideal genome sequence for analysis?
* What are the benefits of illumina sequencing?
* What are the pitfalls of illumina sequencing?
* How do you expect our data to look like?

> We need two components to work: data retrieval and data visualization.

# Data visualization (general)




# Data retrieval

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/205/GCF_000149205.2_ASM14920v2/GCF_000149205.2_ASM14920v2_genomic.gbff.gz

gunzip GCF_000149205.2_ASM14920v2_genomic.gbff.gz
```

```python
from Bio import SeqIO

records = []

for r in SeqIO.parse("GCF_000149205.2_ASM14920v2_genomic.gbff", "genbank"):
    records.append(r)

print(records)
print(len(records))

print(records[0])
print(records[1])
```

* What are we calling with `records[0]` and `records[1]`?


```python
lengths = []
for r in records:
    lengths.append(len(r.seq))

lengths  # What is this?
```

# Exercise 1
<!-- Should take around 20 min -->

> Overall deliverable: Visualize the length of contigs

Let us see how long the contigs are and how many proteins they encode.

Below you see a **list comprehension**. If that does not sound familiar, have a look [here](Getting-Started.md)
* `features` is a list of gene features on the contig
* `f.type` checks whether the `type` *attribute* on the current feature is *CDS*.


```python
r = records[0]  # get a single record as example

# go over features and get the coding sequences
cdss = [f for f in r.features if f.type=="CDS"]

# Information in cds.qualifiers:
# /locus_tag="AN6487.2"
# /protein_id="XP_664091.1"
# /product="hypothetical protein"
```

```python
all_cdss = []
for record in records:
    rdata = {}
    rdata["contig"] = record.id
    rdata["contig_length"] = len(record.seq)
    for feature in record.features:
        if feature.type == "CDS":
            all_cdss.append(rdata | dict(feature.qualifiers))  # python3.9
            # all_cdss.append({**rdata, **dict(feature.qualifiers)})  # python3.8
```

```python
df = pd.DataFrame(all_cdss)
# df.groupby('contig')['product'].agg('count')
df['protein_count'] = df.groupby('contig')['protein_id'].transform('count')

# getting the first item in the dataframe:

df.iloc[0]
sub = df[['contig', 'contig_length', 'protein_count']].drop_duplicates()
sns.violinplot(data=sub,
               y='protein_count')
```

# Exercises

## Visualizing BGC content for two fungi:

* Create a new notebook in the project session_2 on deepnote
* Modify the function we created in the last session to read genbank files to:
  * the features of type protocluster
  * get the product of the protocluster (e.g. PKS, NRPS)
* visualize the counts of different cluster per organism for the files: aspergillus_nidulans.gbk, aspergillus_terreus.gbk


- Python documentation for the packages:
  - http://biopython.org/DIST/docs/tutorial/Tutorial.html
  - https://pandas.pydata.org
  - https://seaborn.pydata.org


## Lovastatin:

Lovastatin is an important cholesterol lowering drug. It’s produced by Aspergillus terreus a fungus in the genus Aspergillus which contains prolific producers of secondary metabolites. The biosynthetic gene clusters (BGC) producing bioactive compounds (such as Lovastatin) were annotated using the web tool antiSMASH (https://antismash-db.secondarymetabolites.org/output/GCF_000149615.1/index.html#). The microbiology team wants to express the BGC to produce Lovastatin. Using the information in the literature, can you point out which BGC produces Lovastatin?

Helpful articles:
* https://onlinelibrary.wiley.com/doi/epdf/10.1002/bip.21428
* https://academic.oup.com/nar/article/49/W1/W29/6274535


## Visualizing genome statistics and biosynthetic gene cluster data

Start a blank deepnote notebook and work on the exercises below. You should reuse code that we used in previous session. Be sure to use `Code` cells when you want to execute code and `Markdown` (You can find markdown by creating a new cell with `+Block` and then clicking on `Markdown`) cells when you want to document your steps.


### Create an overview of genome statistics

Create an overview of genome statistics for Aspergillus and Penicillium species from the [antiSMASH databse](https://antismash-db.secondarymetabolites.org).

*  Visit the [browse section of the antismash database](https://antismash-db.secondarymetabolites.org/browse.html)
* Go through the taxonomy to get the Aspergillus and Penicillium species. If you don't know the taxonomy you can look it up on [NCBI](https://www.ncbi.nlm.nih.gov/taxonomy)
  * Download the corresponding antismash files for *three* organisms of your choice by
    * clicking on the organism name in the antismash overview
    * Now that you are in the overview, click on the download icon and select *Download genbank summary file*
    * Rename the file to something that's easy to use in scripts, e.g. the species name
    * upload the file to deepnote
* Create an overview of genome statistics for each organism (does not need to be combined), you choose how to plot it or create a table similar to *Table 1* in [this publication](https://www.pnas.org/content/115/4/E753)

Conditions:
* It's sufficient to only display Genome size in Mbp, number of proteins, number of contigs and number of scaffolds > 2kbps
  * Go through the [pandas tutorial](https://pandas.pydata.org/docs/getting_started/intro_tutorials/index.html), especially [Dataframe aggregation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.agg.html) and [Series aggregation](https://pandas.pydata.org/docs/reference/api/pandas.Series.agg.html)
  * to calculate the counts you'll need to select the [correct columns/ subset of the dataframe](https://pandas.pydata.org/docs/getting_started/intro_tutorials/03_subset_data.html) use pandas functions like `drop_duplicates()` and `agg`

For this task you'll need to use the functions we created to read the CDS data from *Biopython* `records`.

Extra task:
* How can you show the generated statistics for all organisms in the same table/ plot?

###  Create an overview of secondary metabolite genes

The files you downloaded from the antiSMASH database not only contain CDS annotation, but also antismash cluster annotations.

* Use the functions extracting the antiSMASH *"product"* from the Biopython `record`.
* For each organism, create a table or a barplot using `seaborn` of the *"products"*.

Extra task:

Extra task:
* How can you show the generated statistics for all organisms in the same table/ plot?

### Find the promising drug candidate Malformin

df_tubingensis_sub[
    df_tubingensis_sub['nrps_length'] > 5000
    ]

> Aspergillus tubingensis is a producer of Malformin: how do we find the gene cluster responsible for Malformin production?

Questions:
* What domain is specific to NRPS enzymes?
* How do we find it in antismash data?
* How many NRPSs does Aspergillus tubingensis contain?
* How do we know what NRPS will produce what non-ribosomal peptide?
* What's the number of amino acids for a single NRPS module?

Solution:
> How to find the Malformin gene cluster in Aspergillus tubingensis:
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0140069
# See e.g. https://pubchem.ncbi.nlm.nih.gov/compound/Malformin-A#section=3D-Status Malformin A

Malformin is a pentapeptide, thus we can infer NRPS-size as a criterium:
*  Cys-Cys-Val-Leu-eIle
<!-- * cyclo[DL-cysteinyl-DL-cysteinyl-DL-valyl-DL-leucyl-DL-isoleucyl] (1->2)-disulfide -->
# NRPS proteins have an AMP-binding domain. Try to search for it

# df_tubingensis_sub.size  total number of AMP-binding containing proteins
# solution: https://antismash-db.secondarymetabolites.org/output/GCF_013340325.1/index.html#r5c3

Python code (partial)

```python
tubingensis = list(SeqIO.parse(
    Path.home().joinpath("deepnote/aspergillus_tubingensis.gbk"), "genbank")
    )
df_tubingensis = get_df(tubingensis)
df_tubingensis['org'] = 'A. tubingensis'

df_tubingensis_sub = df_tubingensis[df_tubingensis.sec_met_domain.str.contains('AMP-binding', na=False)][['locus_tag', 'contig_id', 'protein_id', 'translation']]

df_tubingensis_sub['nrps_length'] = df_tubingensis_sub.translation.str.len()  # length should be 4-5 * 1200 AA for each module

df_tubingensis_sub[
    df_tubingensis_sub['nrps_length'] > 5000
    ]
```

The four candidates are:
* XP_035354456.1
* XP_035354521.1
* XP_035357159.1
* XP_035359919.1
