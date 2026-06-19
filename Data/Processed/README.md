# Processed data

Deze map bevat verwerkte databestanden die zijn gebruikt of gemaakt tijdens de transcriptomics-analyse van reumatoïde artritis versus controle.

## Bestanden

### `count_matrix_RA.txt`

Dit is de officiële count matrix die is gebruikt als input voor de differentiële genexpressieanalyse met `DESeq2`.

In deze matrix staan genen als rijen en samples als kolommen. De kolommen bestaan uit vier controle-samples en vier RA-samples:

* `SRR4785819`
* `SRR4785820`
* `SRR4785828`
* `SRR4785831`
* `SRR4785979`
* `SRR4785980`
* `SRR4785986`
* `SRR4785988`

De sample-informatie die bij deze kolommen hoort staat in `Data/Metadata/sample_metadata_RA.csv`.

### `RA_countmatrix.csv`

Dit is een CSV-versie van `count_matrix_RA.txt`. Dit bestand is opgeslagen voor overzicht en controle, zodat de count matrix makkelijker geopend kan worden in bijvoorbeeld Excel.

### `RA_countmatrix_NCBI_GRCh38p14_subset40k.csv`

Dit is een zelfgemaakte count matrix die is verkregen na het mappen van de subset40k FASTQ-bestanden tegen het humane NCBI RefSeq referentiegenoom GRCh38.p14 en het tellen van reads met `featureCounts`.

Deze count matrix is gebruikt om te laten zien hoe de count matrix vanuit gemapte reads gemaakt kan worden. Voor de uiteindelijke DESeq2-analyse is `count_matrix_RA.txt` gebruikt.

### `RA_countmatrix_NCBI_GRCh38p14_subset40k_annotation.csv`

Dit bestand bevat de annotatie-output van `featureCounts` die hoort bij de zelfgemaakte NCBI count matrix. Hierin staat aanvullende informatie over de features/genen die zijn gebruikt bij het tellen van reads.

### `featureCounts_NCBI_GRCh38p14_subset40k.rds`

Dit is het volledige R-object dat door `featureCounts` is aangemaakt. Dit bestand is opgeslagen als RDS-bestand, zodat de volledige featureCounts-output later opnieuw in R kan worden ingeladen zonder de telling opnieuw uit te voeren.

## Belang van deze map

De bestanden in deze map vormen de overgang tussen de ruwe/gemapte sequencingdata en de statistische analyse. De ruwe FASTQ-bestanden en BAM-bestanden zijn grote bestanden en staan niet op GitHub. De verwerkte count matrices zijn wel toegevoegd, omdat deze nodig zijn om de differentiële genexpressieanalyse reproduceerbaar te maken.

## Gebruik in de analyse

Voor de uiteindelijke differentiële genexpressieanalyse is het bestand `count_matrix_RA.txt` gebruikt. Dit bestand is ingelezen in R, omgezet naar een integer count matrix en gebruikt als input voor `DESeq2`.

De overige bestanden zijn opgeslagen als extra verwerkte output en documentatie van de mapping- en featureCounts-stap.
