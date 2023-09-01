# dNdS
Calculation of ratio non sinomic to sinonimic nucleotide substitution (dN/dS ratio)

Procedure contains three steps:

1 - calculation table of probability of substitution between two codons
2 - calculation of substitution between all sequences and referal sequence
3 - calculation of dN/dS ratio of each gene



Step one

Input table: codone_Cap_Low_tab_no-stop
Output table dNdS.tab

Set of spetial instruction:
dN_comb_tab

Run procedure (gawk should be installed):
gawk -f calculation_dN_sites.awk codone_Cap_Low_tab dN_comb_tab > dNdS.tab



Step two

Input tables: dNdS.tab sequences.aln
Output table: sequences.txt

Run procedure
Parameters: name (ancessor genome) name_ref (referal genome)
gawk -f dNdS_simple_cor.awk name=Ng87 name_ref=Ng87 dNdS.tab sequences.aln > sequences_dNdS_ind.tab


Step three

Input: table of substitution of each pair of genes (with referal genome)
Output: dN/dS table

Run of procedure 
gawk -f dNdS_dist_gene_max.awk sequences_dNdS.txt > sequences_dNdS.tab
