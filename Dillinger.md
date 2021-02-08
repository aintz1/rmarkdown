---
title: "NGS assesment 2"
author: "Aintzane Díaz Mazquiarán"
date: "8th February 2021"
output: 
  html_document:
    toc: yes
    toc_depth: 5
    toc_float:
        collapsed: yes
        smooth_scroll: yes
---
```{bash setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, tidy = TRUE, warning = FALSE, eval=FALSE)
```
# Coursework 2: Bacterial genome comparison - identifying virulence-associated genes
### Background

We have received the following letter:
*Dear Universidad de Navarra students,*
*We are terribly grateful for your advice on antibiotic prescriptions for our E. coli AFPN02.1 infected patients. After a few doses we have controlled the progression of the outbreak in most of the infected - your advice is saving human lives! Unfortunately, we have some extreme cases where the infection has systemically progressed and now is unstoppable with normal prescriptions. We need you to further identify the nature of the E. coli AFPN02.1 virulence-associated genes to treat our patients.*

For this task, we have been asked to help identify **virulence-associated genes** in the *E. coli* AFPN02.1 genome (the strain responsible for a big outbreak in Germany some years ago) that could form targets of antibiotic drugs. Virulence genes aid the spread of infection and are thus good candidates for targeting by drugs. To avoid targeting non-pathogenic bacteria and to limit the potential for evolution of resistance, drugs should target genes that are either unique or differ substantially in sequence to homologues in closely related organisms.

Therefore, our aim in this coursework will be to study virulence-associated genes in the AFPN02.1 *E. coli* strain using the tools proposed on the Genomic Comparison practical, so that we find a potential therapeutic target to avoid the resistance to antibiotics of these bacteria, ensuring a better management of the outbreak.

## Question 1 

In this part of the coursework, we are asked to compile a set of virulence-related genes in the AFPN02.1 strain and other *E. coli* strains and compare them, following the same approach we took in the session 3 practical. 

- ##### First task
The first task to be done will be to build a set of virulence-associated genes from any source or scientific publication. But first of all, let's explain what **virulence factors** are. Virulence factors refer to the properties (i.e., gene products) that enable a microorganism to establish itself on or within a host of a particular species and enhance its potential to cause disease. Among virulence factors we can find bacterial toxins, cell surface proteins that mediate bacterial attachment, cell surface carbohydrates and proteins that protect a bacterium, and hydrolytic enzymes that may contribute to the pathogenicity of the bacterium. The genes that encode for these virulence factors are called **virulence-associated genes**. 

In this coursework, our bacterium of study belongs to an *Escherichia coli* strain. Even if good non-pathogenic *E. coli* strains exist, there are several pathogenic *E. coli* bacteria that cause various diseases in humans, including several types of diarrhea, urinary tract infections, sepsis, and meningitis.

*E. coli* strains that cause human diarrhea of varying severity have been divided into six major categories: enterotoxigenic *E. coli* (**ETEC**), enteroinvasive *E. coli* (**EIEC**), enteropathogenic *E. coli* (**EPEC**), enterohemorrhagic *E. coli* (**EHEC**), enteroaggregative *E. coli* (**EAEC**), and diffusely adherent *E. coli* (**DAEC**). Urinary tract infections (UTIs) are the most common extraintestinal *E.coli* infections and are caused by uropathogenic *E.coli* (**UPEC**). In addition, *E. coli* is the most common Gram-negative bacterium that causes meningitis, particularly during the neonatal period. The pathotype responsible for meningitis and sepsis is called neonatal meningitis-associated *E.coli* (**NMEC**).
 
In our case problem, the bacterium strain to which our sequence belongs is an *E. coli* serotype O104:H4, which seems to be an unusual strain similar to EAEC of serotype O104:H4, with the difference that there is also the presence of a prophage encoding the **Shiga toxin**, a characteristic of EHEC strains. Therefore, in this bacterium there is a **combination of features associating both EAEC and EHEC**, representing a new pathotype. 

In the following table, we can find a set of virulence-related genes in the AFPN02.1 strain and other *E. coli* strains, downloaded from VFDB, a virulence factor database (http://www.mgc.ac.cn/cgi-bin/VFs/compvfs.cgi?Genus=Escherichia):
IMAGEEEEEEEEEN

Following the aim of this first task, we have downloaded the FASTA sequences of several virulence-associated genes from the VFDB database (http://www.mgc.ac.cn/VFs/main.htm) to the *course_materials* file. Let's take a look at the VFDB file:
```{bash}
# Create and change directory
mkdir -p /home/manager/Desktop/course_materials/results_assessment2/virulence_factors
cd /home/manager/Desktop/course_materials/results_assessment2/virulence_factors

# Move VFDB file from share_folder to just created directory
mv  /home/manager/Desktop/sf_share_folder/VFDB_setA_nt.fas .

# Read file obtained from VFDB server containing virulence-associated genes in fasta format
less VFDB_setA_nt.fas
```
IMAGEEEEEEEEEN

Here we can observe how virulence-associated genes are ordered by their corresponding header and FASTA sequence. To better examine which genes and how many are listed in this file, we will run the following commands:
```{bash}
####################################################
# Examine VFDB file
####################################################
# Take a look at the headers
cat VFDB_setA_nt.fas | grep '^>' 

# Count the number of genes listed
cat VFDB_setA_nt.fas | grep '^>' | wc -l
```
Now we know that there are 3688 virulence-associated genes listed in this file. In the header for each one of them, we see there is more information than needed. In fact, we can see that there are virulence factors corresponding to species which are not of interest, as we are only interested in genes which belong to *E. coli* strains. To obtain just the genes belonging to *E. coli*, we will make use of commands including **cat** and **grep**. In other words, we will parse virulence-associated genes in *E. coli*: 
```{bash}
# Output all the lines that contain "Escherichia" string
cat VFDB_setA_nt.fas | grep 'Escherichia' 
# Number of genes belonging to Escherichia coli
cat VFDB_setA_nt.fas | grep 'Escherichia' | wc -l
# There are 309 genes. Create a new file with virulence-associated genes of E. coli appearing the first line of the FASTA sequence
cat VFDB_setA_nt.fas | grep -A1 'Escherichia' > Ecoli_vf.fasta
```
IMAGEEEEEEEEEEEEEN
We can see there are 309 virulence-associated genes that belong to the *E. coli* bacterium. 
Nevertheless, we will put special focus on 6 selected virulence-associated genes, for which I have searched some information: shiga-like toxin 2 subunit B (stx2B), Hemolysin A (hlyA), heat-stable enterotoxin 1 (astA), fimbrial chaperone FaeE (faeE), polysialic acid transport protein (kpsM), and UDP-N-acetylglucosamine 2-epimerase (neuC). The selection could be done using grep -E, but it does not work in this case, perhaps because the header is not separated in different columns. Thus, we will do the selection by using the grep function as many times as necessary.
```{bash}
# Select stb2B fasta sequence
cat Ecoli_vf.fasta | grep -A1 'stx2B' > selected_Ecoli_vf.fasta
# Select Hemolysin A fasta sequence
cat Ecoli_vf.fasta | grep -A1 'Hemolysin A' >> selected_Ecoli_vf.fasta
# Select astA fasta sequence
cat Ecoli_vf.fasta | grep -A1 'astA' >> selected_Ecoli_vf.fasta
# Select faeE fasta sequence
cat Ecoli_vf.fasta | grep -A1 'faeE' >> selected_Ecoli_vf.fasta
# Select kpsM fasta sequence
cat Ecoli_vf.fasta | grep -A1 'kpsM' >> selected_Ecoli_vf.fasta
# Select neuC fasta sequence
cat Ecoli_vf.fasta | grep -A1 'neuC' >> selected_Ecoli_vf.fasta
```
IMAGEEEEEEEEEEEEEN
As we can see, there are 6 different sequences, one for each virulence gene. Furthermore, bear in mind there is no "--" included between sequences. Thus, it will not be necessary to remove it. 
To finish with this first section, we need to clean and modify the fasta header to make it more readable. We will just need the ">" fasta symbol and the abbreviation of the gene name.

```{bash}
# Cut the part where the gene abbreviation appears and create new intermediate file
cut -d ' ' -f2 selected_Ecoli_vf.fasta > intermediate1.fasta 
```
IMAGEEEEEEEN
```{bash}
# Change "(" by the fasta symbol ">" and create another intermediate file
cat intermediate1.fasta | sed 's/(/>/g' > intermediate2.fasta
```
IMAGEEEEEEEN
```{bash}
# Change ")" by nothing and create final file
cat intermediate2.fasta | sed 's/)//g' > final_selected_Ecoli_vf.fasta
```
IMAGEEEEEN
In this image we can see the final list of virulence-associated genes in FASTA format.



- ##### Second task
The second part of this first question will be to create a set of virulence-associated genes from our predicted annotation file. Firstly, we will obtain a consensus sequence for the new strain. In the NGS2 practical we created a fasta-format-consensus sequence using the programs *freebayes* and *vcf2fasta*. This new consensus sequence is called *AFPN02.1_merge_consensus.fa*. We will do a **less** *AFPN02.1_merge_consensus.fa* to have a first impression of the sequence:

```{bash}
# Change directory to fastq folder containing file of interest
cd /home/manager/Desktop/course_materials/results_NGS2
# see fasta format
less AFPN02.1_merge_consensus.fa
```
IMAGEEEEEEEN


We have previously attached some of the virulence-associated genes that are common in *E. coli* strains. Nevertheless, the sequence we are working on can or can not have some of those genes. To know if those genes are present in this sequence, and to know the exact place where they are located, we have to do an annotation of the consensus sequence. For that, we can use the server DFAST, which produces a number of files containing the predicted annotation of genes in our genomic sequence (genbank, gff...). However, we have been given an already annotated sequence, *annotation.zip*, which we will be using during this coursework.

```{bash}
# Go to annotation directory
cd /home/manager/Desktop/course_materials/results_NGS2/annotation
```

Once in the annotation folder, we can observe different annotation formats.

```{bash}
## Explore the file formats present on the annotation folder

head annotation.gbk
head annotation.gff
head cds.fna
head features.tsv
head genome.fna
head protein.faa
head rna.fna
head statistics.txt
```

Each one of these formats stores a specific type of information:
- The **GFF annotation format** (*annotation.gbk*) is designed to store annotation, genomic position, features, metadata, taxonomy, sequence and reference authors.
IMAGEEEEEEEEN
- The **GFF annotation format** (*annotation.gff*) is designed to store genomic position, annotation and metadata.
IMAGEEEEEEEEN
- The **CDS.fna annotation format** (*cds.fna*) is amultifasta format with protein sequence.
IMAGEEEEEEEN
- The **features.tsv annotation format** (*features.tsv*) consists on a table with genomic coordinates and minimal annotation.
IMAGEEEEEEEN


Next step will be the key step which will allow us fulfill our goal. In fact, here we will identify virulence-associated genes and examine whether they are present in our strain and other strains of *E. coli*. For that, we will use the annotation.gff file. Firstly, we will **grep** the virulence-associated genes. 

```{bash}
# Grep lines with "virulence"
cat annotation.gff | grep -E 'virulence' 
# Grep several virulence-associated genes' names and transform to bed with awk
cat annotation.gff | grep -E 'virul|stx|hemoly|astA|faeE|kpsM|neuC' | awk 'BEGIN {FS="\t"} split($9, captured, /[(=);]/) >=10 {print "sequence1" "\t" $4 "\t" $5 "\t" captured[10] "\t" captured[2] "\t" $7}'             
```
IMAGEEEEEEEN

```{bash}
# Create bed file
cat annotation.gff | grep -E 'virul|stx|hemoly|astA|faeE|kpsM|neuC' | awk 'BEGIN {FS="\t"} split($9, captured, /[(=);]/) >=10 {print "sequence1" "\t" $4 "\t" $5 "\t" captured[10] "\t" captured[2] "\t" $7}' > present_in_AFPN02_virulence_genes.bed   
# Move bed file to virulence_factors/ folder
mv present_in_AFPN02_virulence_genes.bed /home/manager/Desktop/course_materials/results_assessment2/virulence_factors/
# Go to virulence_factors directory
cd /home/manager/Desktop/course_materials/results_assessment2/virulence_factors/
```

Once we have the annotation and genomic position, we need to include the fasta sequence. For that, we will make use of bedtools.
```{bash}
# Get fasta sequence with bedtools
bedtools getfasta -name -s -fi /home/manager/Desktop/course_materials/results_NGS2/annotation/genome.fna -bed present_in_AFPN02_virulence_genes.bed -fo present_in_AFPN02_virulence_genes.fasta
# To clean header and modify it to ">" and abbreviation of gene
cat present_in_AFPN02_virulence_genes.fasta | sed 's/\:\:.*//g' > present_in_AFPN02_virulence_genes_acronym.fasta
```
IMAGEEEEEEN

- ##### Third task
In this last task of Question 1, we will use BRIG to visualise which of the virulence genes are present/absent in *E. coli* strains (always following the approach we took in the practical session). First, we will merge the **selected** genes from the VFDB and those queried from the annotation (**present**) in a final file: "**final_comparison_virulence.fasta**".
```{bash}
# Merge virulence files
cat present_in_AFPN02_virulence_genes_acronym.fasta selected_Ecoli_vf.fasta > final_comparison_virulence.fasta
# See final fasta file
less final_comparison_virulence.fasta
```
IMAGEEEEEEEN
```{bash}
# Copy all the new files to wholeGenomeExamples folder
cp present_in_AFPN02_virulence_genes_acronym.fasta selected_Ecoli_vf.fasta final_comparison_virulence.fasta /home/manager/Desktop/course_materials/genomes/wholeGenomeExamples
```
Once having done this merge, it is time to compare virulence-associated genes across different genomes using BRIG. BRIG (BLAST Ring Image Generator) is used to visualize the results of comparing virulence-associated genes in our strain and other *E. coli* strains. The input will be the concatenated list of sequences we created in the orevious sections. 
```{bash}
# Create folder to output BRIG results
mkdir -p /home/manager/Desktop/course_materials/results_assessment2/BRIG_output
# Run BRIG software
cd /home/manager/Desktop/BRIG-0.95-dist 
/usr/bin/java -Xmx1500M -jar BRIG.jar
```
IMAGEEEEEEEN (BRIG1)
Here we have selected the reference sequence. As the aim of this circular BLAST comparison is to look for homology between our reference sequence **final_comaprison_virulence.fasta** obtained in the previous steps and our *E. coli* sequence from the Germany outbreak AFPN02.1, as well as other known *E. coli* bacterial genomes.
The steps to follow are the following:
1. Set the reference sequence as "final_comaprison_virulence.fasta".
2. Set "/home/manager/Desktop/course_materials/genomes/wholeGenomeExamples" as the query sequence folder.
3. Press "add to the pool" to load several items into the pool list.
4. Set the output folder as "/home/manager/Desktop/course_materials/results_assessment2/BRIG_output".
5. Make sure the BLAST option box is blank.
6. Click "Next".
IMAGEEEEEEEEEEN (BRIG2)

Then, we have to configure the rings. For that, we have to set the spacer as 50 and set the legend text for each ring. Next, we have to select the required sequences from the data pool and click "add data" to add. After that, we have to choose  a colour, set upper (90) and lower (70) identity thresholds, and click "add new ring" and repeat the same steps for each ring, which will correspond to another new strain to be compared.

IMAGEEEEEEEEEN (BRIG3)

Finally, we will add one more ring and leave the legend empty, selecting "Add custom features". There, we will select "Multi-FASTA" as imput sequence, change colour to "alternating red-blue", draw features as clockwise-arrow, press "Add" just once, and press "Close".

IMAGEEEEEEEN (BRIG4)

Finally, we will select output format as png, rename the image output, select "Re-do BLAST", and submit.

IMAGEEEEEEN (BRIG5)

The BRIG image obtained after this process is the following:

IMAGEEEEEN !["brig rings" width="200"/>]("file://BRIGcomparison.png")

In this BLAST we have compared the 7th annotation ring (with virulence-associated genes in blue and red) against all the bacterial genomes included in the analysis. Having a hit (where colours of different bacterial species appear) suggests that a homologue of this gene is present on that bacterium (with a similarity of about 70-90%). Thus, returning back to this coursework's aim, if we are asked for drug treatment to kill this new strain, we should recommend drugs that target virulence-associated genes that have a hit on the AFPN02.1 ring, but not in the other ones, in the case that the other *E. coli* strains are not pathogenic and do not provoke any problmatic symptoms). Thus, before making any decisions, we must know a little bit more about the genomes we have chosen to do the comparison.
- The **UTI89** and **CFT073** are uropathogenic *E.coli* (UPEC) strains (https://www.uniprot.org/proteomes/UP000001952, https://pubmed.ncbi.nlm.nih.gov/27214553/). 
- The **O126** is an enteroaggregative *E. coli* strain, a diarrheagenic agent (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3016780/).
- The **HS** was isolated from a healthy human with no disease, and it
has been shown to colonize the human gastrointestinal tract without any apparent clinical symptoms. Thus, it can be considered as a non-pathogenic strain (https://www.genome.jp/kegg-bin/show_organism?org=ecx).
- The **O157** is an enterohemorrhagic strain considered as a major foodborne pathogen causing severe disease in humans worldwide. Healthy cattle are a reservoir of E. coli O157:H7, and bovine food products and fresh produce contaminated with bovine waste are the most common sources for disease outbreaks in the United States. E. coli O157:H7 also survives well in the environment. The abilities to cause human disease, colonize the bovine gastrointestinal tract, and survive in the environment require that E. coli O157:H7 adapt to a wide variety of conditions. Three major virulence factors of E. coli O157:H7 have been identified including Shiga toxins, products of the pathogenicity island called the locus of enterocyte effacement, and products of the F-like plasmid pO157 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3645889/). This strain seems to share similar characteristics of our AFPN02.1 strain.
- The **k12** is used as a model organism strain.

Having said this, we now can conclude that k12 and HS strains should not be killed in this case, whereas the UTI89, CFT073 and O126 strains might be better to target, and without any doubt, the O157 strain should be killed along with our AFPN02.1 German strain, as these are the most pathogenic strains among all. Therefore, looking again at the comaprison BRIG ring, targeting the virulent **Shiga toxins 2A and 2B** encoding genes, which are present in *E. coli* O157 and in our strain, as well as targeting the **sigA** and/or **bfpE** genes, which are mostly present just in our strain, would be the best strategies to stop the outbreak.



## Question 2

- ##### Task 1
The first task in this second question is to select one single gene present in AFPN02.1 and relevant to virulence, and study it in more detail. Firstly, we are asked to find information on the biologial action/mechanism of the gene.
The gene I have selected is stx2A, which encodes for the subunit B of the Shiga-toxin type 2 virulence factor. 

Several serotypes of Escherichia coli produce one or more protein toxins that are closely related to Stx from Shigella dysenteriae serotype 1 (1). Stx from S. dysenteriae was first identified by Kiyoshi Shiga for whom these toxins are named (2). Collectively these E. coli are known as Stx-producing E. coli (STEC).1 E. coli O157:H7 is the STEC responsible for many outbreaks of hemorrhagic colitis or bloody diarrhea in the U.S., Canada, and Japan. In some patients the prodrome of hemor- rhagic colitis may progress to the hemolytic uremic syndrome, which culminates in kidney failure. Two types of Stx may be produced by STEC, Stx1 and/or Stx2. The Stx types all have an AB5 structure, in which a single A-subunit is associated with five B-subunits. The A-subunit embodies the N-glycosidase catalytic activity; it acts by removing a specific adenine base from the 28 S rRNA of the 60 S ribosomal subunit within infected cells. Because this adenine base is on a loop of rRNA that is important for elongation factor binding, the toxin is able to shut down protein synthesis in a targeted cell (https://pubmed.ncbi.nlm.nih.gov/15075327/).

S. dysenteriae and the Stx were identified in the 19th century by Drs. Neisser and Shiga (1) and Conradi (2). Approximately 80 years later the same toxin (now called Stx1 to distinguish it from the toxin produced by S. dysenteriae) was found in a group of E. coli isolates. These bacteria caused bloody diarrhea and a serious sequelae, the hemolytic uremic syndrome (HUS), a condition characterized by thrombocytopenia, hemolytic anemia, and kidney failure (3, 4). 

HUS cases were reported in the literature as early as the 1950s, though the cause was unknown. An understanding of the origins of infection-associated HUS was also complicated by the fact that some cases of the HUS are of genetic origin (see review (112)). However, in 1983 bloody diarrhea and the HUS were linked to certain serogroups of E. coli and O’Brien et al. found that those strains produced a toxin related to the Stx of S. dysenteriae (3, 4, 113). The strong association between the HUS and Stx2a was underscored in the largest outbreak of the HUS in Germany in 2011 in which more than 800 cases of the disease were identified (114). The surprise from German outbreak was that although the implicated strain was an E. coli, unlike typical STEC, the epidemic isolate also encoded virulence factors associated with enteroaggregative E. coli (EAEC). EAEC had only rarely been associated with the HUS previously, and in each of those cases the implicated isolate was found to carry an Stx gene.

##### Structure of Shiga toxins
Shiga toxin family members have an AB5 molecular configuration, as revealed by X-ray crystallography19,20: An enzymatically active monomeric A subunit, StxA (which has a molecular mass of 32 kDa) is non-covalently associated with a pentamer of identical B fragments (each B fragment has a molecular mass of 7.7 kDa) that form the B subunit, StxB, which is responsible for binding to cell surface receptors (Fig. 1a). StxB forms a doughnut-shaped structure with a central pore into which the carboxyl terminus of StxA inserts19. StxA and the StxB fragments are secreted into the bacterial periplasm21, where they assemble non-covalently into the holotoxin, as was initially described for heat-labile enterotoxins from E. coli22 (https://www.nature.com/articles/nrmicro2279).
IMAGEEEEEE ![sadasdasf](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnrmicro2279/MediaObjects/41579_2010_Article_BFnrmicro2279_Fig1_HTML.jpg?as=webp)
a | A cartoon of Shiga holotoxin, consisting of one A subunit (StxA), which is cleaved into fragments A1 and A2, and five B fragments that constitute the homopentameric B subunit (StxB). b | A ribbon diagram of Shiga toxin, highlighting globotriaosylceramide (Gb3)-binding sites on StxB. Gb3 is shown in a ball-and-stick representation. c | An enlargement of StxA at the site of furin cleavage (Arg25-Met252), and showing the disulphide bond (between Cys242 and Cys261) that links the A1 and A2 framents. d | A ribbon diagram of an StxB subunit from the membrane-oriented surface, highlighting the three Gb3-binding sites. Gb3 is shown as a ball-and-stick representation. Note the central pore that is lined by α-helices.


##### Mechanism of action of Shiga toxin (Stx) 
The Stxs (also known as Vero toxins, and previously as Shiga-like toxins) are a group of bacterial AB5 protein toxins of about 70 kDa that inhibit protein synthesis in sensitive eukaryotic cells. Stx is constituted by a pentamer of B subunits bound to a catalytic A subunit. Protein synthesis is blocked by the Stxs through the removal of an adenine residue from the 28S rRNA of the 60S ribosome. This N-glycosidase activity of the toxin resides in the A subunit. The pentamer of identical B subunits mediates toxin binding to the cellular receptor globotriaosylceramide (Gb3) (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4270005/). The B subunits bind to globotriaosylceramide (Gb3) expressed by some eukaryotic cells (1) Stx is internalized by endocytosis (2) Subsequently, Stx undergoes retrograde transport to the trans-Golgi network (TGN) (3) and then to the endoplasmic reticulum (ER) (4) In the ER, Stx encounters its target, the ribosome, inactivating it (4) As a consquence, Stx inhibits protein synthesis, causing cell death by apoptosis (https://pubmed.ncbi.nlm.nih.gov/22919672/).

![saundaiuf ](https://www.ncbi.nlm.nih.gov/pmc/articles/instance/3417539/bin/fcimb-02-00081-g0001.jpg)
Shiga toxin is the prototype of the Shiga toxin family and nearly identical to the E. coli-produced Shiga toxin 1 (Stx1), differing by a single amino acid in the catalytic A subunit of the toxin. STEC can produce Stx1 variants (Stx1 and Stx1c), Stx2 variants (Stx2, Stx2c, Stx2d, Stx2e, Stx2f) or variants of both in a range of combinations (Table 1). However, severe disease has been epidemiologically linked to the presence of Stx2 (Ref. 4). Although Stx1 and Stx2 share a common receptor and possess the same intracellular mechanism of action, they are immunologically distinct and only 56% identical at the amino acid sequence level5. Stx2 variants are 84–99% homologous to Stx2.

- ##### Task 2
Once described the gene selected, the following step will be to study its phylogeny based on BLAST searches against other microorganisms. 

Thus, as BLAST analysis is required, stx2B FASTA sequence will be selected.
```{bash}
# Go to directory of interest
cd course_materials/results_assessment2/virulence_factors
# Grep "stx2B" gene of interest in "present_in_AFPN02_virulence_genes_acronym.fasta" file
cat present_in_AFPN02_virulence_genes_acronym.fasta | grep -A1 'stx2B'
```
IMAGEEEEEEEEN

Here we have our gene of interest and its FASTA sequence. The nucleotide-sequence will be copied to clipboard. Next, we will open blastn web server (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch), and we will paste the FASTA sequence there, which has been extracted from our consensus sequence.
The organisms selected to perform BLAST are: Vibrio, Escherichia coli, Spiribacter, Archaea, Coccidia, Bacilli, and Salmonella. Then, we select “discontiguous megablast” and “show results in a new window”.
IMAGEEEEEEEEN (blast1)

As we can see in the image, in the tree appear strains belonging to *E. coli*. So we will also do a PSI-BLAST with the proteic sequence of stx2B to do iterative analysis and see if we can observe a more complex phylogenetic tree.
The protein sequence has been retrieved from Uniprot  (https://www.uniprot.org/uniprot/A7UQX3) and we have introduced the sequence in BLASTp suite (https://blast.ncbi.nlm.nih.gov/Blast.cgi):
IMAGEEEN
This time, the organisms of interest to be taken into account are Bacteria, Bacillus/Clostridium group, Salmonella enterica, and Bacillus/Staphylococcus group. For each iteration, we will obtain certain information, and we have retrieved the following taxonomic information:
IMAGEEEEN
Finally, this is the phylogenic tree we have captured:
IMAGEEEEEN
As we can observe, stx2B virulence factor is a protein that mostly appears in enterobacteria strains (in green), particularly in *E. coli* strains. There are also stx2B toxins belonging to b-proteobacteria organisms, especifically to *Burkholderia ambifiaria* and *Burkholderia pseudomallei*. Based on this tree, it looks like the Burkholderia was the first type of organism in evolution which adquired the capacity to encode this toxin, and later on, thanks to evolution, enterobacteria like *E. coli* gained the ability to produce stx2B. 











