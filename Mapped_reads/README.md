# Mapped reads

Deze map is bedoeld voor de gemapte reads die ontstaan na het alignen van de paired-end FASTQ-bestanden tegen het humane referentiegenoom.

In deze analyse zijn de reads gemapt met het R-package `Rsubread`. De output van deze stap bestaat uit BAM-bestanden. Deze BAM-bestanden geven weer waar de reads op het referentiegenoom zijn uitgelijnd.

De BAM-bestanden zijn tussenbestanden in de RNA-seq workflow. Ze zijn nodig om vervolgens met `featureCounts` een count matrix te maken. De count matrix is daarna gebruikt als input voor de differentiële genexpressieanalyse met `DESeq2`.

De daadwerkelijke BAM-, BAI- en eventuele SAM-bestanden zijn niet toegevoegd aan GitHub, omdat dit grote bestanden zijn. Ze zijn lokaal bewaard. De stappen waarmee deze bestanden zijn gemaakt, zijn vastgelegd in het R-script in de map `Scripts`.

Voorbeelden van bestanden die lokaal in deze map kunnen staan:

* `SRR4785819.bam`
* `SRR4785819.bam.bai`
* `SRR4785820.bam`
* `SRR4785820.bam.bai`

Deze bestanden worden uitgesloten van GitHub via `.gitignore`.
