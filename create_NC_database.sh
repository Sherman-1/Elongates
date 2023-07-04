#!/bin/bash

# Create non-coding databases for the different species 
# From GTF and Fasta, create a databaase from everything not annotated as CDS

filter_CDS() {
 
    local input_gff="$1"

    # Check if the input file exists
    if [[ ! -f "$input_gff" ]]; then
    echo "Error: Input GFF file not found."
        return 1
    fi

    # Filter out only the CDS features using Awk
    awk '$3 == "gene"' "$input_gff" > "gene.gff"

}

sizes_gff() {

    local fasta="$1"
    local flag="$2"

    if [[ ! -f "$fasta" ]]; then
    echo "Error: Input fasta file not found."
        return 1
    fi

    faSize -detailed "$fasta" > "${fasta%.*}.sizes"

    # Create a GFF file with the coordinates of the chromosomes
    awk -v OFS='\t' '{print $1, "noncoding", "exon", 1, $2, ".", "+", ".", "ID=" $1}' "${fasta%.*}.sizes" > "genome.gff"
    rm "${fasta%.*}.sizes"

} 

species=$(yq '.Species_order.Scer_for_bash' env.yaml)
cd input || exit
for specie in $species; do

    filter_CDS ${specie}.gff 
    sizes_gff ${specie}.fna 


    bedtools subtract -a "genome.gff" -b "CDS.gff" > "tmp.gff" && rm "gene.gff" "genome.gff"
    bedtools getfasta -fi ${specie}.fna -bed "tmp.gff" -fo "nc.fna" && rm "tmp.gff"
    faTrans nc.fna ${specie}_nc.faa && rm "nc.fna"

    # Add the specie name to the fasta header
    awk -i inplace -v var="$specie" '/^>/ {$0=$0" "var} 1' ${specie}_nc.faa

done

cat *_nc.faa > full_nc.faa 
makeblastdb -in full_nc.faa -dbtype prot -out interCDS_DB && rm *_nc.faa