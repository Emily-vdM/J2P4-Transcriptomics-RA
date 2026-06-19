# Reference files

Deze map is bedoeld voor de referentiebestanden die nodig zijn voor het indexeren en mappen van de RNA-seq reads.

In deze analyse is het humane NCBI RefSeq referentiegenoom gebruikt:

* Referentiegenoom: NCBI RefSeq GRCh38.p14
* Assembly accession: GCF_000001405.40_GRCh38.p14
* Bron: NCBI Datasets
* Link: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

Voor het mappen en tellen van reads moeten het referentiegenoom en de annotatie van dezelfde bron en versie afkomstig zijn. Daarom zijn in deze analyse zowel het FASTA-bestand als het GTF-annotatiebestand van NCBI RefSeq GRCh38.p14 gebruikt.

De daadwerkelijke FASTA-, GTF- en indexbestanden zijn niet toegevoegd aan GitHub, omdat deze bestanden groot zijn. Ze zijn lokaal bewaard en gebruikt tijdens de analyse. In het R-script staat beschreven hoe de referentie is gebruikt voor het indexeren en mappen.
