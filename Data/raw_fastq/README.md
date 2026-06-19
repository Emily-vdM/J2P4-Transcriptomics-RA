# Raw FASTQ data

Deze map is bedoeld voor de ruwe paired-end FASTQ-bestanden die zijn gebruikt voor de RNA-seq analyse van synoviumbiopten van controlepersonen en RA-patiënten.

In deze analyse zijn acht samples gebruikt:

| Sample ID  | Groep          | Leeftijd | Geslacht |
| ---------- | -------------- | -------: | -------- |
| SRR4785819 | Controle       |       31 | female   |
| SRR4785820 | Controle       |       15 | female   |
| SRR4785828 | Controle       |       31 | female   |
| SRR4785831 | Controle       |       42 | female   |
| SRR4785979 | RA established |       54 | female   |
| SRR4785980 | RA established |       66 | female   |
| SRR4785986 | RA established |       60 | female   |
| SRR4785988 | RA established |       59 | female   |

De FASTQ-bestanden zijn paired-end bestanden. Per sample hoort er daarom een read 1- en read 2-bestand bij, bijvoorbeeld:

* `SRR4785819_1_subset40k.fastq`
* `SRR4785819_2_subset40k.fastq`

De daadwerkelijke FASTQ-bestanden zijn niet toegevoegd aan GitHub, omdat sequencingbestanden groot zijn en niet geschikt zijn om rechtstreeks in de repository op te slaan. De bestanden zijn lokaal bewaard en gebruikt als input voor het mappen van reads met Rsubread.

De gebruikte sample-informatie staat opgeslagen in `Data/Metadata/sample_metadata_RA.csv`.
