#!/bin/bash

#################################################################
#                                                               #
#      Preparation for Circos conf file for Circos plot         #
#  Written by Desirrê Petters-Vandresen & Bruno Janoski Rossi   # 
#                in the 1st Trimester of 2021                   #
#                Last modified = 26/03/2021                     #
#                                                               #
#################################################################

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo 
   echo "                         Circos prep script                             "
   echo
   echo "Written by Desirrê Petters-Vandresen & Bruno Janoski Rossi in March 2021"
   echo "Last modified: 26/03/2021"
   echo
   echo "For updates, check https://github.com/desirrepetters"
   echo
   echo "-------------------------------------------------------------------------"
   echo
   echo "This is a script to prepare and produce Circos plots using the software Circos (http://www.circos.ca/software/)."
   echo
   echo "The final Circos plots will show scaffolds/contigs for a strain, and prediction tracks."
   echo
   echo "Prediction tracks available at the moment include:"
   echo
   echo "  - Gene density (gene annotation in GFF3 format, adapted)"
   echo "  - TE density (TE annotation in GFF3 format, here using outputs from the REPET pipeline"
   echo "  - Effector density (Effector list, here using output from the effector prediction pipeline of Feurtey and Lorrain et al. 2020 (DOI: 10.1186/s12864-020-06871-w)"
   echo "  - Biosynthetic gene clusters density (BGCs list, here using output from Antismash, the geneclusters.txt file)"
   echo "  - CAZyme density (CAZyme list, here using output from dbCAN2 and considering only CAZymes predicted by >= 2 tools"
   echo "  - Species-specific gene density and strain-specific gene density (Ortholog groups table, here using output file from ProteinOrtho)"
   echo
   echo "The script can be adapted to include more tracks or remove some of the tracks (by commenting respective lines)"
   echo
   echo "The script also relies on additional software to process the files. Before running you will need to install and add to your path (or edit the script informing the paths before the commands) the following softwares (besides installing Circos):"
   echo
   echo "  - BEDTOOLS (https://bedtools.readthedocs.io/en/latest/)"
   echo "  - pyfaidx (https://pypi.org/project/pyfaidx/, easily installed with pip)"
   echo
   echo "-------------------------------------------------------------------------"
   echo
   echo "  USAGE: circosprep.sh -s [STRAIN] -t [STRAIN TAG] -p [SPECIES] -c [CUTOFF]"
   echo
   echo "  example: circosprep.sh -s PCapitalensis_CBS128856 -t Cap128_ -p PCapitalensis -c ALL"
   echo
   echo "  options:"
   echo
   echo "  -s     Strain name. For the script to work, the strain name has to be exactly the same in all the required annotation files, ortholog table headers and so on, and works as an identifier to the strain to be analysed. For example, a strain name could be PCapitalensis_CBS128856, where PCapitalensis represents the species and CBS128856 is the code for the specific strain"
   echo
   echo "  -t     Strain tag. This tag included in some output files and tracks allows the files to be re-used in synteny Circos plots, for example. A good tag is short but contains useful information to identify the strain. For example, for the strain PCapitalensis_CBS128856 a possible tag could be Cap128_ or Cap128856_. Always remember to finish the tag with an underline"
   echo
   echo "  -p     Species. Including the species code helps parsing the orthologs table later and identifying isolates belonging to the same species of the strain being analysed. For the script to work, all strains from the same species in the ortholog table need to have the same species name. For example, in the table you would have PCapitalensis_CBS128856 and PCapitalensis_CBS356_52 representing two isolates from the species PCapitalensis, so in the -p option you would state PCapitalensis"
   echo
   echo "  -c     Cutoff. Use this to limit the total number of scaffolds/contigs to be included in the plot when the assembly is too fragmented. If you do not want to set a limit, use -c ALL. Scaffolds/contigs in the final plot will be ordered from largest to smallest, so, for example, if you set -c as 10, the Circos plot will show only the 10 largest scaffolds/contigs. Note that this can mean showing scaffolds/contigs 1 to 10 if your original assembly was already ordered and numbered from largest to smallest, however, if you original assembly was not ordered this way, it will show the largest contigs, whatever their number may be."
   echo
   echo "  -*     Print this Help."
   echo
}

################################################################################

while getopts "s:t:p:c:" flag;
do
    case "${flag}" in
        s) STRAIN=${OPTARG};;
        t) STRAIN_TAG=${OPTARG};;
        p) SPECIES=${OPTARG};;
		c) SET_CUTOFF=${OPTARG};;
		*) Help
		exit 1 ;;
    esac
done

# Genome, annotation and prediction foldes
# Remember to adjust to corresponding paths in your system
GENOME_PATH="/PATH/TO/GENOME_FILE/"
GENE_ANNOTATION_PATH="/PATH/TO/GENE_ANNOTATION_FILE/"
TE_PATH="/PATH/TO/TE_ANNOTATION_FILE/"
TELOMERE_PATH="/PATH/TO/TELOMERE_ANNOTATION_FILE/"
EFFECTOR_PATH="/PATH/TO/EFFECTOR_LIST_FILE/"
CAZYME_PATH="/PATH/TO/CAZYME_LIST_FILE/"
BGC_PATH="/PATH/TO/BGC_LIST_FILE/"
ORTHOLOG_PATH="/PATH/TO/ORTHOLOGS_FILE/"
TEMPLATE_PATH="/PATH/TO/CIRCOS_CONF_TEMPLATE_FILE/"

# Path to working directory
WORK_DIR="/PATH/TO/WORK_DIR/"

# Working directories for better organization
# These directories are optional, but if you comment this section, remember to change the rest of the script accordingly

SLIDING_WINDOW_DIR="10Kb sliding windows (BED files)/"
BGC_ANNOTATION_DIR="BGC Annotation (BED files)/"
BGC_ANNOTATION_SUPPORT_DIR="BGC Annotation (BED files)/Support files/"
CAZYME_ANNOTATION_DIR="CAZymes Annotation (BED files)/"
CAZYME_ANNOTATION_SUPPORT_DIR="CAZymes Annotation (BED files)/Support files/"
CIRCOS_DIR="Circos plots/"
CIRCOS_SUPPORT_DIR="Circos plots/Configuration files/"
CIRCOS_PLOT_DIR="Circos plots/Plots/"
EFFECTOR_ANNOTATION_DIR="Effector Annotation (BED files)/"
EFFECTOR_ANNOTATION_SUPPORT_DIR="Effector Annotation (BED files)/Support files/"
FILE_PREP_DIR="File preparation/"
FILE_PREP_BED_DIR="File preparation/Basic BED file/"
FILE_PREP_KARYOTYPE_HEAD_DIR="File preparation/Karyotype head/"
FILE_PREP_KARYOTYPE_TAIL_DIR="File preparation/Karyotype tail/"
FILE_PREP_SCAFFOLD_COLUMN_DIR="File preparation/Scaffold column/"
FILE_PREP_SCAFFOLD_TAG_DIR="File preparation/Scaffold tag/"
GENE_ANNOTATION_DIR="Gene Annotation (BED files)/"
GENOME_BED_DIR="Genomes (BED files)/"
INTERSECT_BGC_DIR="Intersects of BGCs vs. 10kb sliding windows/"
INTERSECT_CAZYME_DIR="Intersects of CAZymes vs. 10kb sliding windows/"
INTERSECT_EFFECTOR_DIR="Intersects of effectors vs. 10kb sliding windows/"
INTERSECT_GENE_DIR="Intersects of genes vs. 10kb sliding windows/"
INTERSECT_SPECIES_SPECIFIC_DIR="Intersects of species specific genes vs. 10kb sliding windows/"
INTERSECT_STRAIN_SPECIFIC_DIR="Intersects of strain specific genes vs. 10kb sliding windows/"
INTERSECT_TE_DIR="Intersects of TEs vs. 10kb sliding windows/"
KARYOTYPE_DIR="Karyotypes files/"
KARYOTYPE_SUPPORT_DIR="Karyotypes files/Support files/"
SPECIES_SPECIFIC_DIR="Species-specific genes (BED files)/"
SPECIES_SPECIFIC_SUPPORT_DIR="Species-specific genes (BED files)/Support files/"
STRAIN_SPECIFIC_DIR="Strain-specific genes (BED files)/"
STRAIN_SPECIFIC_SUPPORT_DIR="Strain-specific genes (BED files)/Support files/"
TELOMERE_DIR="Telomeres (BED files)/"
TELOMERE_SUPPORT_DIR="Telomeres (BED files)/Support files/"
TE_DIR="TEs Annotation (BED files)/"

[ -d "${WORK_DIR}${SLIDING_WINDOW_DIR}" ] && echo "The folder ${SLIDING_WINDOW_DIR} already exists!" || mkdir "${WORK_DIR}${SLIDING_WINDOW_DIR}"
[ -d "${WORK_DIR}${BGC_ANNOTATION_DIR}" ] && echo "The folder ${BGC_ANNOTATION_DIR} already exists!" || mkdir "${WORK_DIR}${BGC_ANNOTATION_DIR}"
[ -d "${WORK_DIR}${BGC_ANNOTATION_SUPPORT_DIR}" ] && echo "The folder ${BGC_ANNOTATION_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${BGC_ANNOTATION_SUPPORT_DIR}"
[ -d "${WORK_DIR}${CAZYME_ANNOTATION_DIR}" ] && echo "The folder ${CAZYME_ANNOTATION_DIR} already exists!" || mkdir "${WORK_DIR}${CAZYME_ANNOTATION_DIR}"
[ -d "${WORK_DIR}${CAZYME_ANNOTATION_SUPPORT_DIR}" ] && echo "The folder ${CAZYME_ANNOTATION_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${CAZYME_ANNOTATION_SUPPORT_DIR}"
[ -d "${WORK_DIR}${CIRCOS_DIR}" ] && echo "The folder ${CIRCOS_DIR} already exists!" || mkdir "${WORK_DIR}${CIRCOS_DIR}"
[ -d "${WORK_DIR}${CIRCOS_SUPPORT_DIR}" ] && echo "The folder ${CIRCOS_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${CIRCOS_SUPPORT_DIR}"
[ -d "${WORK_DIR}${CIRCOS_PLOT_DIR}" ] && echo "The folder ${CIRCOS_PLOT_DIR} already exists!" || mkdir "${WORK_DIR}${CIRCOS_PLOT_DIR}"
[ -d "${WORK_DIR}${EFFECTOR_ANNOTATION_DIR}" ] && echo "The folder ${EFFECTOR_ANNOTATION_DIR} already exists!" || mkdir "${WORK_DIR}${EFFECTOR_ANNOTATION_DIR}"
[ -d "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}" ] && echo "The folder ${EFFECTOR_ANNOTATION_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_DIR}" ] && echo "The folder ${FILE_PREP_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_BED_DIR}" ] && echo "The folder ${FILE_PREP_BED_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_BED_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_KARYOTYPE_HEAD_DIR}" ] && echo "The folder ${FILE_PREP_KARYOTYPE_HEAD_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_KARYOTYPE_HEAD_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_KARYOTYPE_TAIL_DIR}" ] && echo "The folder ${FILE_PREP_KARYOTYPE_TAIL_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_KARYOTYPE_TAIL_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}" ] && echo "The folder ${FILE_PREP_SCAFFOLD_COLUMN_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}"
[ -d "${WORK_DIR}${FILE_PREP_SCAFFOLD_TAG_DIR}" ] && echo "The folder ${FILE_PREP_SCAFFOLD_TAG_DIR} already exists!" || mkdir "${WORK_DIR}${FILE_PREP_SCAFFOLD_TAG_DIR}"
[ -d "${WORK_DIR}${GENE_ANNOTATION_DIR}" ] && echo "The folder ${GENE_ANNOTATION_DIR} already exists!" || mkdir "${WORK_DIR}${GENE_ANNOTATION_DIR}"
[ -d "${WORK_DIR}${GENOME_BED_DIR}" ] && echo "The folder ${GENOME_BED_DIR} already exists!" || mkdir "${WORK_DIR}${GENOME_BED_DIR}"
[ -d "${WORK_DIR}${INTERSECT_BGC_DIR}" ] && echo "The folder ${INTERSECT_BGC_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_BGC_DIR}"
[ -d "${WORK_DIR}${INTERSECT_CAZYME_DIR}" ] && echo "The folder ${INTERSECT_CAZYME_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_CAZYME_DIR}"
[ -d "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}" ] && echo "The folder ${INTERSECT_EFFECTOR_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}"
[ -d "${WORK_DIR}${INTERSECT_GENE_DIR}" ] && echo "The folder ${INTERSECT_GENE_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_GENE_DIR}"
[ -d "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}" ] && echo "The folder ${INTERSECT_SPECIES_SPECIFIC_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}"
[ -d "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}" ] && echo "The folder ${INTERSECT_STRAIN_SPECIFIC_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}"
[ -d "${WORK_DIR}${INTERSECT_TE_DIR}" ] && echo "The folder ${INTERSECT_TE_DIR} already exists!" || mkdir "${WORK_DIR}${INTERSECT_TE_DIR}"
[ -d "${WORK_DIR}${KARYOTYPE_DIR}" ] && echo "The folder ${KARYOTYPE_DIR} already exists!" || mkdir "${WORK_DIR}${KARYOTYPE_DIR}"
[ -d "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}" ] && echo "The folder ${KARYOTYPE_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}"
[ -d "${WORK_DIR}${SPECIES_SPECIFIC_DIR}" ] && echo "The folder ${SPECIES_SPECIFIC_DIR} already exists!" || mkdir "${WORK_DIR}${SPECIES_SPECIFIC_DIR}"
[ -d "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}" ] && echo "The folder ${SPECIES_SPECIFIC_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}"
[ -d "${WORK_DIR}${STRAIN_SPECIFIC_DIR}" ] && echo "The folder ${STRAIN_SPECIFIC_DIR} already exists!" || mkdir "${WORK_DIR}${STRAIN_SPECIFIC_DIR}"
[ -d "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}" ] && echo "The folder ${STRAIN_SPECIFIC_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}"
[ -d "${WORK_DIR}${TELOMERE_DIR}" ] && echo "The folder ${TELOMERE_DIR} already exists!" || mkdir "${WORK_DIR}${TELOMERE_DIR}"
[ -d "${WORK_DIR}${TELOMERE_SUPPORT_DIR}" ] && echo "The folder ${TELOMERE_SUPPORT_DIR} already exists!" || mkdir "${WORK_DIR}${TELOMERE_SUPPORT_DIR}"
[ -d "${WORK_DIR}${TE_DIR}" ] && echo "The folder ${TE_DIR} already exists!" || mkdir "${WORK_DIR}${TE_DIR}"

# Neccessary files

GENOME="${GENOME_PATH}${STRAIN}.fasta"
GENE_ANNOTATION="${GENE_ANNOTATION_PATH}${STRAIN}.gff3"
TES="${TE_PATH}${STRAIN}.gff3"
TELOMERE="${TELOMERE_PATH}${STRAIN}.bed"
EFFECTOR="${EFFECTOR_PATH}${STRAIN}.effectors.tab"
CAZYME="${CAZYME_PATH}${STRAIN}.txt"
BGC="${BGC_PATH}${STRAIN}.txt"
ORTHOLOG="${ORTHOLOG_PATH}Orthology_groups.tabular"

TEMPLATE_CONF="${TEMPLATE_PATH}Template.conf"
TEMPLATE_TELOMERE_CONF="${TEMPLATE_PATH}Template_for_telomeres.conf"

CIRCOS_PATH=/path/to/circos-0.69-9/bin/circos
CIRCOS_CONF="${CIRCOS_SUPPORT_DIR}${STRAIN}.conf"
CIRCOS_TELOMERE_CONF="${CIRCOS_SUPPORT_DIR}${STRAIN}_telomeres.conf"
CIRCOS_OUTPUT_FILE="${STRAIN}"
CIRCOS_OUTPUT_TELOMERE_FILE="${STRAIN}_telomeres"

# Adjusting cutoff to work on sed
CUTOFF=$SET_CUTOFF

#############
#Basic files#
#############

#Using pyfaidx to create a pre-genome file and then format it to the GENOME format from BEDTOOLS

faidx --transform bed "${GENOME}" | sort -Vk1 > "${WORK_DIR}${FILE_PREP_BED_DIR}${STRAIN}_intermediate.bed"

# Creating GENOME file using pyfaidx file

cut -f1,3 "${WORK_DIR}${FILE_PREP_BED_DIR}${STRAIN}_intermediate.bed" > "${WORK_DIR}${GENOME_BED_DIR}${STRAIN}.bed"

# Creating KARYOTYPE file using the pyfaidx file

sed "s/^/$STRAIN_TAG/" "${WORK_DIR}${FILE_PREP_BED_DIR}${STRAIN}_intermediate.bed" | sed "s/^/-\t/" | sed "s/^/chr\t/" | cut -f1-3 > "${WORK_DIR}${FILE_PREP_KARYOTYPE_HEAD_DIR}${STRAIN}_karyotype_head.bed"

wc -l < "${WORK_DIR}${FILE_PREP_BED_DIR}${STRAIN}_intermediate.bed" >  "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}${STRAIN}_scaffold_count.bed"

CHR_COUNT=$(cat "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}${STRAIN}_scaffold_count.bed")

seq 1 "$CHR_COUNT" > "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}${STRAIN}_scaffold_sequence.bed"

cut -f2-3 "${WORK_DIR}${FILE_PREP_BED_DIR}${STRAIN}_intermediate.bed" | sed "s/$/\tvvdrgrey/" > "${WORK_DIR}${FILE_PREP_KARYOTYPE_TAIL_DIR}${STRAIN}_karyotype_tail.bed"

paste "${WORK_DIR}${FILE_PREP_KARYOTYPE_HEAD_DIR}${STRAIN}_karyotype_head.bed" "${WORK_DIR}${FILE_PREP_SCAFFOLD_COLUMN_DIR}${STRAIN}_scaffold_sequence.bed" "${WORK_DIR}${FILE_PREP_KARYOTYPE_TAIL_DIR}${STRAIN}_karyotype_tail.bed" | sort -r -Vk6 > "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes"

#Adjusting karyotype file for the Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$8 = $6+100000; $9 = $7; print}' "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes" | cut -f1-5,8-9 > "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}_for_telomeres.karyotypes"

#Creating file for chromosome order to be included in Circos conf file

cut -f3 "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes" > "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt"

tr '\n' ',' <  "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt" > "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt"

sed -i 's/,/, /g' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt"
sed -i '$ s/.$//' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt"
sed -i '$ s/.$//' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt"

###########
#Telomeres#
###########

# Creating telomere files to be included as highlights in the Circos plot and adding strain tag

cut -f1-3 "${TELOMERE}" | sed "s/^/$STRAIN_TAG/" > "${WORK_DIR}${TELOMERE_DIR}${STRAIN}_telomeres.bed"

awk 'BEGIN {OFS="\t"}; { $2 = ($2 <= 1000 ? 1 : "TAIL"); $3 = ($3 <= 1000 ? 50000 : "TAIL")} 1' "${WORK_DIR}${TELOMERE_DIR}${STRAIN}_telomeres.bed" | sed '/TAIL/!d' > "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_for_tail.bed"

sort -Vk1 "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_for_tail.bed" > "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_for_tail_sorted.bed"

awk 'BEGIN {OFS="\t"}; { $2 = ($2 <= 1000 ? 1 : "TAIL"); $3 = ($3 <= 1000 ? 50000 : "TAIL")} 1' "${WORK_DIR}${TELOMERE_DIR}${STRAIN}_telomeres.bed" | sed '/TAIL/d' > "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_head.bed"

# Sorting the karyotpe file for join

sort -k3 "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes" > "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_sorted.karyotypes"

join -e- -a1 -a2 -1 1 -2 3 -o 1.1,2.5,2.6 -t "$(printf '\t')" --nocheck-order "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_for_tail_sorted.bed" "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes" | sed '/-/d'| awk 'BEGIN {OFS="\t"}; {$4 = $3+50000; $5 = $3+100000; print}' | cut -f1,4-5 > "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_tail.bed"

cat "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_head.bed" "${WORK_DIR}${TELOMERE_SUPPORT_DIR}${STRAIN}_telomeres_tail.bed" > "${WORK_DIR}${TELOMERE_DIR}${STRAIN}_telomeres_adjusted.bed"

#################
#Sliding windows#
#################

# Making 10Kb windows with BEDTOOLS

bedtools makewindows -g "${WORK_DIR}${GENOME_BED_DIR}${STRAIN}.bed" -w 10000 > "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed"

sort -k1,1 -k2,2n "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" > "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb_sorted.bed"

###############################################
#    Creating feature files and intersects    #
###############################################

# Converting gene annotation from GFF3 to BED format

convert2bed --input=GFF < "${GENE_ANNOTATION}" > "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}.bed"

# Sorting (lexical)

sort -Vk4 -t "$(printf '\t')" "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}.bed" > "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_sorted.bed"

# Sorting (numerical)

sort -k4 -t "$(printf '\t')" "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}.bed" > "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_regular_sorted.bed"

# Sorting for windows

sort -Vk1 -Vk3 "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}.bed" > "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_windows_sorted.bed"

# Intersect gene annotation in BED format with the sliding windows and calculate gene density

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" -b "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_windows_sorted.bed" -c -sorted > "${WORK_DIR}${INTERSECT_GENE_DIR}${STRAIN}_gene_density.bed"

# Adding karyotype tag to the gene density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_GENE_DIR}${STRAIN}_gene_density.bed"

# Adjusting gene density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_GENE_DIR}${STRAIN}_gene_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_GENE_DIR}${STRAIN}_gene_density_for_telomeres.bed"


#############
# Effectors #
#############

sed '/#/d' "${EFFECTOR}" > "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}${STRAIN}_effector_without_header.tab"

sort -Vk1 "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}${STRAIN}_effector_without_header.tab" > "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}${STRAIN}_effector_without_header_sorted.tab"

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2 -t "$(printf '\t')" "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_sorted.bed" "${WORK_DIR}${EFFECTOR_ANNOTATION_SUPPORT_DIR}${STRAIN}_effector_without_header_sorted.tab" | sed '/Effector/!d' | cut -f1-10 > "${WORK_DIR}${EFFECTOR_ANNOTATION_DIR}${STRAIN}_effectors.bed"

# Intersect effectors with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" -b "${WORK_DIR}${EFFECTOR_ANNOTATION_DIR}${STRAIN}_effectors.bed" -c -sorted > "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}${STRAIN}_effector_density.bed"

# Adding karyotype tag to the effector density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}${STRAIN}_effector_density.bed"

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}${STRAIN}_effector_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_EFFECTOR_DIR}${STRAIN}_effector_density_for_telomeres.bed"


#########################
# Transposable elements #
#########################

# Converting TE annotation from GFF3 to BED format

convert2bed --input=GFF < "${TES}" > "${WORK_DIR}${TE_DIR}${STRAIN}.bed"

# Removing SSRs from the annotation file

sed '/REPET_SSRs/d' "${WORK_DIR}${TE_DIR}${STRAIN}.bed" > "${WORK_DIR}${TE_DIR}${STRAIN}_no_SSR.bed"

# Intersect TEs with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb_sorted.bed" -b "${WORK_DIR}${TE_DIR}${STRAIN}.bed" -c -sorted > "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density.bed"

# Intersect TEs with sliding windows after removing SSRs

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb_sorted.bed" -b "${WORK_DIR}${TE_DIR}${STRAIN}_no_SSR.bed" -c -sorted > "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density_no_SSR.bed"

# Adding karyotype tag to the TE density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density.bed"

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density_no_SSR.bed"

# Adjusting TE density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density_for_telomeres.bed"

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density_no_SSR.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_TE_DIR}${STRAIN}_TE_density_for_telomeres_no_SSR.bed"

#####################################
# Biosynthetic Gene Clusters (BGCs) #
#####################################

# Creating BGCs file in BED format to intersect and calculate BGC density
 
cut -f4 "${BGC}" | sed -e 's/;/\n&/g;s/;//g' | sort -k1 > "${WORK_DIR}${BGC_ANNOTATION_SUPPORT_DIR}${STRAIN}_BGCs_one_per_line_sorted.txt"

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 -t "$(printf '\t')" -i "${GENE_ANNOTATION_DIR}${STRAIN}_regular_sorted.bed" "${WORK_DIR}${BGC_ANNOTATION_SUPPORT_DIR}${STRAIN}_BGCs_one_per_line_sorted.txt" | sort -Vk4 -t "$(printf '\t')" > "${WORK_DIR}${BGC_ANNOTATION_DIR}${STRAIN}_BGCs.bed"

# Intersect BGCs with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" -b "${WORK_DIR}${BGC_ANNOTATION_DIR}${STRAIN}_BGCs.bed" -c -sorted > "${WORK_DIR}${INTERSECT_BGC_DIR}${STRAIN}_BGC_density.bed"

# Adding karyotype tag to the BGCs density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_BGC_DIR}${STRAIN}_BGC_density.bed"

# Adjusting BGC density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_BGC_DIR}${STRAIN}_BGC_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_BGC_DIR}${STRAIN}_BGC_density_for_telomeres.bed"

###########
# CAZymes #
###########

# Creating CAZymes file in BED format to intersect and calculate CAZymes density

sed 's/$/#/g;/1#/d;/HMMER/d' "${CAZYME}" | cut -f1 | sort -k1 > "${WORK_DIR}${CAZYME_ANNOTATION_SUPPORT_DIR}${STRAIN}_cazymes_without_headers_sorted.tab"

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 -t "$(printf '\t')" -i "${GENE_ANNOTATION_DIR}${STRAIN}_regular_sorted.bed" "${WORK_DIR}${CAZYME_ANNOTATION_SUPPORT_DIR}${STRAIN}_cazymes_without_headers_sorted.tab" | sort -Vk4 > "${WORK_DIR}${CAZYME_ANNOTATION_DIR}${STRAIN}_cazymes.bed"

# Intersect CAZymes with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" -b "${WORK_DIR}${CAZYME_ANNOTATION_DIR}${STRAIN}_cazymes.bed" -c -sorted > "${WORK_DIR}${INTERSECT_CAZYME_DIR}${STRAIN}_cazymes_density.bed"

# Adding karyotype tag to the CAZymes density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_CAZYME_DIR}${STRAIN}_cazymes_density.bed"

# Adjusting CAZyme density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_CAZYME_DIR}${STRAIN}_cazymes_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_CAZYME_DIR}${STRAIN}_cazymes_density_for_telomeres.bed"

##########################
# Species-specific genes #
##########################

# Retrieving orthogroups for the strain of interest from ProteinOrtho output file. We need to check the number of the column for the respective strain. Moreover, here we will extract genes that are found in all isolates from the same species of the strain of interest, but are not found in other species (species-specific). For example, if you use 5 isolates of the same species, you will need genes with the gene count equal to or lower than five. 

# Discovering number of column for the strain in original orthologs file

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed "/${STRAIN}.faa/!d" | cut -f1 | sed 's/ //g' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_column_ortholog.txt"

STRAIN_COLUMN_ORTHOLOG=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_column_ortholog.txt")

# Discovering total number of strains in original orthologs file

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed '/Species/d;/Genes/d;/Alg.-Conn./d' | wc -l > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_number_ortholog.txt"

STRAIN_NUMBER_ORTHOLOG=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_number_ortholog.txt")

# Discovering number of strains from the same species of the strain of interest in original orthologs file

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed "/${SPECIES}/!d" | wc -l > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_number_ortholog.txt"

SPECIES_NUMBER_ORTHOLOG=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_number_ortholog.txt")

OCURRENCES_ORTHOLOG=`expr $SPECIES_NUMBER_ORTHOLOG + 2`

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed "/${SPECIES}/!d" | cut -f1 | sed -n 1p | sed 's/ //g' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_first_strain_position_ortholog.txt"

FIRST_STRAIN_POSITION_ORTHOLOG=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_first_strain_position_ortholog.txt")

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed "/${SPECIES}/!d" | cut -f1 | sed -n '$p' | sed 's/ //g' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_last_strain_position_ortholog.txt"

LAST_STRAIN_POSITION_ORTHOLOG=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_last_strain_position_ortholog.txt")

# Extracting genes with gene count equal to or lower than the number of strains from the species of interest

awk -v species_number_ortholog="$SPECIES_NUMBER_ORTHOLOG" 'BEGIN {OFS="\t"}; (NR==1) || ($1 <= species_number_ortholog)' "${ORTHOLOG}" | cut -f1,$FIRST_STRAIN_POSITION_ORTHOLOG-$LAST_STRAIN_POSITION_ORTHOLOG > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes.tab"

# Saving the header for later

sed '/#/!d' "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes.tab" > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes_header.tab"

# Filtering by occurences to remove genes shared with other species

for ((i=1;i<=SPECIES_NUMBER_ORTHOLOG;i++));
do
awk -v var01="$i" 'BEGIN {OFS="\t"}; {(NR>1) && ($1 = var01); print}' "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes.tab" | awk -v var02="$OCURRENCES_ORTHOLOG" 'BEGIN {OFS="\t"}; {$var02 = gsub(/*/,"*"); print}' | awk -v var03="$OCURRENCES_ORTHOLOG" -v var04="$SPECIES_NUMBER_ORTHOLOG" -v var05="$i" -v var06="$(("$SPECIES_NUMBER_ORTHOLOG-$i"))" '($var03 == var06)' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes_only_${i}_occurrences.tab"

done

ls "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes"_only_*| sed 's/^/cat "/;s/$/"/' | sh > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_without_header.tab"

cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes_header.tab" "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_without_header.tab" | sed 's/# Species/#Species/g' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes.tab"

# Getting the list for strain of interest in species specific genes file

sed -e 's/\t/\n/g;' "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_for_specific_genes_header.tab" | nl -nln | sed "/$STRAIN/!d" | cut -f1 >  "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_strain_column_species_specific.txt"

STRAIN_COLUMN_SPECIES_SPECIFIC=$(cat "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${STRAIN}_strain_column_species_specific.txt")

cut -f$STRAIN_COLUMN_SPECIES_SPECIFIC -d "$(printf '\t')" "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes.tab" | sed '/*/d;/faa/d;s/,/\n&/g;s/,//g' | sed '/Species/d' > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_for_${STRAIN}.tab"

sort -k1 "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_for_${STRAIN}.tab" > "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_for_${STRAIN}_sorted.tab"

# Creating BED file

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 -t "$(printf '\t')" -i "${GENE_ANNOTATION_DIR}${STRAIN}_regular_sorted.bed" "${WORK_DIR}${SPECIES_SPECIFIC_SUPPORT_DIR}${SPECIES}_specific_genes_for_${STRAIN}_sorted.tab" | sort -Vk4 -t "$(printf '\t')" > "${WORK_DIR}${SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific.bed"

# Intersect species specific genes with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb.bed" -b "${WORK_DIR}${SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific.bed" -c -sorted > "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific_density.bed"

# Adding karyotype tag to the species specific density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific_density.bed"

# Adjusting species specific genes density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_SPECIES_SPECIFIC_DIR}${STRAIN}_species_specific_density_for_telomeres.bed"


#########################
# Strain-specific genes #
#########################

# Retrieving orthogroups for the strain of interest from ProteinOrtho output file. Here we need to check only the number of the column for the respective strain. 

# Discovering number of column for the strain in original orthologs file
# It would be perfectly possible to use the same variable assigned in the previous section, but I decided not to

sed -e 's/\t/\n&/g;q' "${ORTHOLOG}" | nl | sed "/${STRAIN}.faa/!d" | cut -f1 | sed 's/ //g' > "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_column_ortholog.txt"

STRAIN_SPECIFIC_COLUMN_ORTHOLOG=$(cat "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_column_ortholog.txt")

# Extracting genes with gene count equal to one, as they are supposed to be strain specific

awk 'BEGIN {OFS="\t"}; (NR==1) || ($1 = 1)' "${ORTHOLOG}" | cut -f$STRAIN_SPECIFIC_COLUMN_ORTHOLOG > "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_for_specific_genes.tab"

sort -k1 "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_for_specific_genes.tab" > "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_for_specific_genes_sorted.tab"

join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10 -t "$(printf '\t')" -i "${WORK_DIR}${GENE_ANNOTATION_DIR}${STRAIN}_regular_sorted.bed" "${WORK_DIR}${STRAIN_SPECIFIC_SUPPORT_DIR}${STRAIN}_for_specific_genes_sorted.tab" > "${WORK_DIR}${STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_genes.tab"

# Intersect strain specific genes with sliding windows

bedtools intersect -a "${WORK_DIR}${SLIDING_WINDOW_DIR}${STRAIN}_10kb_sorted.bed" -b "${WORK_DIR}${STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_genes.tab" -c > "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_density.bed"

# Adding karyotype tag to the strain specific density files

sed -i "s/^/$STRAIN_TAG/" "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_density.bed"

# Adjusting strain specific genes density file for Circos plot with telomeres

awk 'BEGIN {OFS="\t"}; {$5 = $2+50000; $6 = $3+50000; $7 = $4; print}' "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_density.bed" | cut -f1,5-7 > "${WORK_DIR}${INTERSECT_STRAIN_SPECIFIC_DIR}${STRAIN}_strain_specific_density_for_telomeres.bed"

##############################
# Circos configuration files #
##############################

if  [[ ${SET_CUTOFF} == "ALL" ]]
then
	echo "Building Circos plots configuration files that include all scaffolds/contigs"

# Without telomeres

	sed "s|STRAIN|${STRAIN}|g;s|_CUTOFF||g;s|CHR_ORDER|$(cat "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt")|;s|FIRST_CHR|$(sed -n '1p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt")|;s|LAST_CHR|$(sed -n '$p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt")|" "${TEMPLATE_CONF}" | sed -e "s|WORKDIR/|$WORK_DIR|" > "${CIRCOS_CONF}"
 
# With telomeres

	sed "s|STRAIN|${STRAIN}|g;s|_CUTOFF||g;s|CHR_ORDER|$(cat "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order.txt")|;s|FIRST_CHR|$(sed -n '1p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt")|;s|LAST_CHR|$(sed -n '$p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt")|" "${TEMPLATE_TELOMERE_CONF}" | sed -e "s|WORKDIR/|$WORK_DIR|" > "${CIRCOS_TELOMERE_CONF}"

else 

	echo "Building Circos plot configuration files that include only the ${CUTOFF} largest contigs"

# Without telomeres
	head -n "$CUTOFF" "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}.karyotypes" > "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}_cutoff.karyotypes"

	head -n "$CUTOFF" "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order.txt" > "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt"

	tr '\n' ',' <  "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt" > "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt"

	sed -i 's/,/, /g' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt"
	sed -i '$ s/.$//' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt"
	sed -i '$ s/.$//' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt"

	sed "s|STRAIN|${STRAIN}|g;s|_CUTOFF|_cutoff|g;s|CHR_ORDER|$(cat "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt")|;s|FIRST_CHR|$(sed -n '1p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt")|;s|LAST_CHR|$(sed -n '$p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt")|" "${TEMPLATE_CONF}" | sed -e "s|WORKDIR/|$WORK_DIR|" > "${CIRCOS_CONF}"

# With telomeres

	head -n "$CUTOFF" "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}_for_telomeres.karyotypes" > "${WORK_DIR}${KARYOTYPE_DIR}${STRAIN}_for_telomeres_cutoff.karyotypes"

	sed "s|STRAIN|${STRAIN}|g;s|_CUTOFF|_cutoff|g;s|CHR_ORDER|$(cat "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_chr_order_cutoff.txt")|;s|FIRST_CHR|$(sed -n '1p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt")|;s|LAST_CHR|$(sed -n '$p' "${WORK_DIR}${KARYOTYPE_SUPPORT_DIR}${STRAIN}_for_chr_order_cutoff.txt")|" "${TEMPLATE_TELOMERE_CONF}" | sed -e "s|WORKDIR/|$WORK_DIR|" > "${CIRCOS_TELOMERE_CONF}"

fi

################
#              #
# Circos plots #
#              # 
################

# Without telomeres

perl ${CIRCOS_PATH} -conf "${CIRCOS_CONF}" -outputdir "${CIRCOS_PLOT_DIR}" -outputfile "${CIRCOS_OUTPUT_FILE}"

# With telomeres

perl ${CIRCOS_PATH}  -conf "${CIRCOS_TELOMERE_CONF}" -outputdir "${CIRCOS_PLOT_DIR}" -outputfile "${CIRCOS_OUTPUT_TELOMERE_FILE}"



