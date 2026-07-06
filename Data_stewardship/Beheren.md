# Beheren

Binnen dit transcriptomics-project is gewerkt met RNA-seq data van synoviumbiopten van controlepersonen en personen met reumatoïde artritis. Als data steward heb ik de projectgegevens zo beheerd dat duidelijk is welke bestanden input, tussenbestand, verwerkte data of resultaat zijn. Hierdoor blijft de analyse beter controleerbaar en reproduceerbaar.

## Mappenstructuur

De repository is ingericht met een vaste mappenstructuur:

```text
Data/
├── Metadata/
├── Processed/
└── raw_fastq/

Mapped_reads/
Reference/
Reference_articles/
Results/
Scripts/
Data_stewardship/
```

`Data/raw_fastq/` is bedoeld voor de ruwe paired-end FASTQ-bestanden. `Mapped_reads/` is bedoeld voor BAM- en indexbestanden na mapping. `Reference/` is bedoeld voor het referentiegenoom, annotatiebestand en Rsubread-index. Deze grote bestanden zijn niet op GitHub geplaatst, maar lokaal bewaard. De mappen blijven wel zichtbaar door een `README.md`-bestand.

`Data/Processed/` bevat verwerkte data zoals count matrices. `Results/` bevat tabellen, figuren, pathway-output en `sessionInfo`. `Scripts/` bevat het R-script waarmee de analyse is uitgevoerd.

## Metadata en bestandsnamen

De sample-informatie staat in `Data/Metadata/sample_metadata_RA.csv` en `Data/Metadata/sample_info_RA_subset40k.csv`. Deze metadata koppelt de sample IDs aan de juiste conditie, waardoor DESeq2 RA correct met controle kan vergelijken.

Bestandsnamen zijn zo gekozen dat de inhoud en vergelijking herkenbaar zijn, bijvoorbeeld:

```text
DESeq2_resultaten_RA_vs_control.csv
GO_Biological_Process_RA_vs_control.csv
Volcanoplot_RA_vs_control.png
```

## Grote bestanden en .gitignore

Grote bestanden zoals FASTQ-, BAM-, BAI-, FASTA-, GTF- en indexbestanden zijn uitgesloten via `.gitignore`. Hierdoor blijft de repository overzichtelijk en geschikt om te delen. Alleen scripts, metadata, verwerkte tabellen, figuren en documentatie zijn toegevoegd aan GitHub.

## Reproduceerbaarheid

De analyse is vastgelegd in `Scripts/transcriptomics_RA_workflow.R`. In dit script staan de gebruikte stappen, packages, bestandslocaties en outputbestanden. Daarnaast is `Results/sessionInfo_transcriptomics_RA.txt` opgeslagen, zodat zichtbaar is met welke R-versie en packageversies de analyse is uitgevoerd.

## Zorgvuldig databeheer

De gebruikte data is afkomstig uit een openbare dataset en wordt aangeduid met SRR-accessienummers. Er zijn geen direct identificeerbare persoonsgegevens in de repository geplaatst. Door ruwe data, scripts, metadata en resultaten gescheiden te beheren, blijft het project overzichtelijk, controleerbaar en reproduceerbaar.
