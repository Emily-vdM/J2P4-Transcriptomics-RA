# Metadata

Deze map bevat de metadata die hoort bij de RNA-seq samples die zijn gebruikt in de transcriptomics-analyse van reumatoïde artritis versus controle.

## Bestanden

### `sample_metadata_RA.csv`

Dit bestand bevat de biologische sample-informatie per sample. Hierin staat onder andere:

* sample ID
* conditie
* diagnose
* leeftijd
* geslacht

De gebruikte samples zijn vier controlepersonen en vier personen met gevestigde reumatoïde artritis. De controlegroep is gebruikt als referentiegroep in de DESeq2-analyse. Hierdoor betekent een positieve `log2FoldChange` dat een gen hoger tot expressie komt in RA ten opzichte van controle.

### `sample_info_RA_subset40k.csv`

Dit bestand bevat de technische sample-informatie die is gebruikt voor de subset40k FASTQ-bestanden. Hierin staat per sample:

* SRR-runnummer
* sample name
* conditie

Dit bestand is gebruikt om de juiste FASTQ-bestanden te koppelen aan de juiste sample-namen en condities tijdens het mappen van reads en het maken van de count matrix.

## Belang van metadata

Metadata is belangrijk omdat de count matrix alleen read counts per sample bevat. Zonder metadata is niet duidelijk welke kolommen bij controlepersonen horen en welke kolommen bij RA-patiënten horen. De metadata zorgt er dus voor dat de differentiële genexpressieanalyse correct kan worden uitgevoerd en dat de analyse reproduceerbaar blijft.

## Gebruikte groepen

In deze analyse zijn de samples als volgt ingedeeld:

| Sample ID  | Conditie |
| ---------- | -------- |
| SRR4785819 | control  |
| SRR4785820 | control  |
| SRR4785828 | control  |
| SRR4785831 | control  |
| SRR4785979 | RA       |
| SRR4785980 | RA       |
| SRR4785986 | RA       |
| SRR4785988 | RA       |
