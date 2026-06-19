# Beheren

## Inleiding

Binnen dit transcriptomics-project is gewerkt met RNA-seq data van synoviumbiopten van controlepersonen en personen met reumatoïde artritis. Bij dit soort projecten is goed databeheer belangrijk, omdat er meerdere soorten bestanden ontstaan: ruwe sequencingdata, gemapte reads, referentiebestanden, count matrices, scripts, figuren en resultaat-tabellen. Wanneer deze bestanden niet duidelijk worden georganiseerd, is het later moeilijk om de analyse terug te vinden, te controleren of opnieuw uit te voeren.

Data stewardship betekent in dit project dat de onderzoeksgegevens op een gestructureerde, begrijpelijke en reproduceerbare manier zijn beheerd. Daarbij is onderscheid gemaakt tussen ruwe data, tussenbestanden, verwerkte data en eindresultaten.

## Mappenstructuur

Voor dit project is een vaste mappenstructuur gebruikt:

```text
Data/
├── Metadata/
├── Processed/
└── Raw_FASTQ/

Mapped_reads/

Reference/

Reference_articles/

Results/
├── Tables/
├── Figures/
└── Pathways/

Scripts/

Documentation/
```

De map `Data/Raw_FASTQ/` is bedoeld voor de ruwe paired-end FASTQ-bestanden. Deze bestanden zijn lokaal gebruikt als input voor het mappen van reads. De map `Mapped_reads/` is bedoeld voor de BAM-bestanden die ontstaan na het alignen van reads tegen het referentiegenoom. De map `Reference/` is bedoeld voor het referentiegenoom, annotatiebestanden en indexbestanden. De map `Data/Processed/` bevat verwerkte data, zoals de count matrix. De map `Results/` bevat de output van de analyse, zoals tabellen, figuren en pathway-resultaten.

Door deze structuur is duidelijk welke bestanden input, tussenbestand of resultaat zijn. Dit maakt het project overzichtelijker en beter reproduceerbaar.

## Bestandsnamen

Voor de naamgeving van bestanden is geprobeerd om korte en duidelijke namen te gebruiken. Lange bestandsnamen zijn vermeden, omdat dit problemen kan geven bij GitHub Desktop en omdat korte namen overzichtelijker zijn. Voor bronartikelen zijn bijvoorbeeld namen gebruikt zoals:

```text
Platzer_2019_RA_gene_expression.pdf
Gabriel_2001_RA_epidemiology.pdf
```

Bij resultaten is in de bestandsnaam opgenomen welke analyse is uitgevoerd en welke vergelijking is gemaakt, bijvoorbeeld:

```text
DESeq2_resultaten_RA_vs_control.csv
significante_genen_RA_vs_control.csv
Volcanoplot_RA_vs_control.png
GO_Biological_Process_RA_vs_control.csv
```

Hierdoor is zonder het bestand te openen al zichtbaar wat de inhoud ongeveer is.

## Metadata

De sample-informatie is apart opgeslagen in `Data/Metadata/sample_metadata_RA.csv`. Hierin staat per sample tot welke groep het sample behoort, bijvoorbeeld controle of RA. Metadata is belangrijk omdat de count matrix alleen de read counts bevat, maar niet automatisch uitlegt welke kolom bij welke conditie hoort.

In dit project is gewerkt met vier controles en vier RA-samples. De metadata is gebruikt om in DESeq2 de vergelijking tussen RA en controle correct uit te voeren.

## Grote bestanden en .gitignore

Niet alle bestanden zijn geschikt om naar GitHub te uploaden. Grote bestanden zoals FASTQ-, BAM-, BAI-, FASTA- en GTF-bestanden zijn lokaal bewaard en uitgesloten via `.gitignore`. Dit voorkomt dat de repository onnodig groot wordt en houdt de GitHub-pagina overzichtelijk.

De mappen `Data/Raw_FASTQ/`, `Mapped_reads/` en `Reference/` zijn wel zichtbaar gemaakt op GitHub door er een `README.md`-bestand in te plaatsen. In deze README-bestanden staat welke bestanden lokaal in deze mappen horen en waarom deze niet op GitHub zijn geplaatst.

## Reproduceerbaarheid

De reproduceerbaarheid is verbeterd door de analyse in R-scripts vast te leggen. In de scripts staan de gebruikte stappen, package calls, bestandslocaties en commentaarregels. Hierdoor kan iemand anders de analyse beter volgen en eventueel opnieuw uitvoeren wanneer dezelfde inputbestanden beschikbaar zijn.

Ook zijn resultaten opgeslagen in vaste mappen, zoals `Results/Tables/` voor tabellen en `Results/Figures/` voor figuren. Hierdoor zijn de outputbestanden makkelijk terug te vinden.

## Omgaan met gevoelige data

De gebruikte data is afkomstig uit een openbaar beschikbare dataset en wordt aangeduid met SRR-accessienummers. Er zijn geen direct identificeerbare persoonsgegevens, zoals namen of adressen, in de repository geplaatst. Wanneer in toekomstig onderzoek wel gevoelige of herleidbare data gebruikt zou worden, zou deze niet openbaar gedeeld mogen worden. In dat geval zouden bestanden alleen op een beveiligde opslaglocatie geplaatst worden, met beperkte toegang voor personen die de data nodig hebben.

## Open data

Open data is belangrijk omdat het onderzoek controleerbaar en herbruikbaar maakt. Bij transcriptomicsonderzoek is het waardevol als ruwe data, metadata en analysecode beschikbaar zijn, zodat andere onderzoekers analyses kunnen controleren of uitbreiden. Tegelijk moet er zorgvuldig worden omgegaan met grote bestanden en eventuele gevoelige data. In dit project is daarom gekozen om de analysecode, metadata, resultaten en documentatie op GitHub te plaatsen, terwijl grote ruwe databestanden lokaal zijn bewaard en via documentatie worden beschreven.
