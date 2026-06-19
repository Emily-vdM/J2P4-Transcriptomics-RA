# J2P4-Transcriptomics-RA

RNA-seq analyse van synoviumbiopten van personen met reumatoïde artritis en controlepersonen met behulp van `Rsubread`, `DESeq2`, KEGG en Gene Ontology analyse.

---

## Inhoud/structuur

Deze repository bevat de uitwerking van de transcriptomics-casus voor J2P4. De analyse is uitgevoerd in R en de repository is gestructureerd zodat data, scripts, resultaten, referenties en beheerdocumentatie gescheiden zijn opgeslagen.

* `Data/Metadata/` – metadata over de gebruikte samples. In deze map staat ook een README met uitleg over de metadata-bestanden.
* `Data/Processed/` – verwerkte data, zoals de officiële count matrix en zelfgemaakte featureCounts-output. In deze map staat ook een README met uitleg over de bestanden.
* `Data/Raw_FASTQ/` – map voor ruwe FASTQ-bestanden. De bestanden zelf zijn niet toegevoegd aan GitHub.
* `Mapped_reads/` – map voor lokaal opgeslagen BAM-bestanden. De bestanden zelf zijn niet toegevoegd aan GitHub.
* `Reference/` – map voor referentiegenoom-, annotatie- en indexbestanden. De bestanden zelf zijn niet toegevoegd aan GitHub.
* `Reference_articles/` – gebruikte literatuurbronnen.
* `Results/Tables/` – resultaat-tabellen uit de analyse.
* `Results/Figures/` – figuren, zoals de volcano plot en GO-dotplot.
* `Results/Pathways/` – KEGG/pathview-output.
* `Scripts/` – R-scripts voor de transcriptomics-analyse.
* `Data_stewardship/` – uitleg over data stewardship en GitHub-beheer.

---

## Inleiding

Reumatoïde artritis, afgekort RA, is een systemische auto-immuunziekte waarbij ontsteking van het gewrichtsslijmvlies een belangrijke rol speelt. Deze ontsteking wordt ook wel synovitis genoemd en kan uiteindelijk leiden tot gewrichtsschade. De precieze oorzaak van RA is niet volledig bekend, maar genetische aanleg, omgevingsfactoren en een verstoord immuunsysteem spelen waarschijnlijk samen een rol.

In deze casus wordt RNA-seq data gebruikt om verschillen in genexpressie te onderzoeken tussen personen met RA en controlepersonen. Met RNA-seq kan op grote schaal worden bepaald welke genen meer of minder tot expressie komen. Hierdoor kan inzicht worden verkregen in biologische processen en pathways die betrokken zijn bij het ziekteproces.

De gebruikte dataset bestaat uit synoviumbiopten van vier controlepersonen en vier personen met gevestigde RA. De personen met RA zijn ACPA-positief en de controlepersonen zijn ACPA-negatief. Het doel van deze analyse is om met behulp van R te onderzoeken welke genen differentieel tot expressie komen in RA ten opzichte van controle, en welke KEGG-pathways en Gene Ontology-termen hierbij betrokken zijn.

### Onderzoeksvraag

Welke genen, biologische processen en pathways verschillen tussen synoviumbiopten van personen met reumatoïde artritis en controlepersonen?

---

## Methode

De analyse is uitgevoerd in R. De ruwe sequencingdata bestond uit paired-end FASTQ-bestanden van acht samples. Eerst zijn de reads gemapt tegen het humane referentiegenoom met het R-package `Rsubread`. Voor deze analyse is het humane NCBI RefSeq referentiegenoom GRCh38.p14 gebruikt, met assembly accession `GCF_000001405.40_GRCh38.p14`. De FASTA- en GTF-bestanden zijn van dezelfde bron gebruikt, zodat het referentiegenoom en de annotatie met elkaar overeenkomen.

Na het mappen zijn de BAM-bestanden gebruikt om met `featureCounts` een count matrix te maken. Deze zelfgemaakte count matrix is opgeslagen in `Data/Processed/RA_countmatrix_NCBI_GRCh38p14_subset40k.csv` en laat zien hoe de count matrix vanuit gemapte reads kan worden opgebouwd.

Voor de uiteindelijke differentiële genexpressieanalyse is het bestand `Data/Processed/count_matrix_RA.txt` gebruikt als officiële input voor `DESeq2`. Deze count matrix bevat per gen het aantal reads per sample. De count matrix is ingelezen in R, omgezet naar een integer matrix en gekoppeld aan de samplemetadata. In DESeq2 is de conditie `control` ingesteld als referentiegroep. Hierdoor betekent een positieve `log2FoldChange` dat een gen hoger tot expressie komt in RA ten opzichte van controle.

Significante genen zijn geselecteerd op basis van `padj < 0.05` en `|log2FoldChange| > 1`. Daarna zijn de resultaten geannoteerd met Entrez ID’s met behulp van `org.Hs.eg.db` en `AnnotationDbi`. Voor visualisatie is een volcano plot gemaakt met `EnhancedVolcano`. Daarnaast is met `pathview` de KEGG-pathway voor reumatoïde artritis gevisualiseerd.

Voor de Gene Ontology analyse is gebruikgemaakt van `clusterProfiler::enrichGO()`. Hiermee is een over-representation analysis uitgevoerd voor de drie GO-ontologieën: Biological Process (BP), Molecular Function (MF) en Cellular Component (CC). Bij deze methode wordt onderzocht of bepaalde GO-termen vaker voorkomen in de lijst met significant differentieel tot expressie komende genen dan verwacht op basis van een achtergrondlijst.

Als input voor de GO-analyse zijn de significant veranderde genen gebruikt met `padj < 0.05` en `|log2FoldChange| > 1`. Deze genen zijn omgezet van gene symbols naar Entrez ID’s met `org.Hs.eg.db` en `AnnotationDbi`. Als achtergrond zijn alle betrouwbaar geteste genen gebruikt, dus de genen zonder `NA` in `padj`. Hierdoor wordt de GO-analyse niet vergeleken met het volledige humane genoom, maar met de genen die in deze RNA-seq analyse daadwerkelijk betrouwbaar getest konden worden.

Er is gekozen voor `clusterProfiler::enrichGO()` omdat dit package geschikt is voor functionele enrichmentanalyse en visualisatie van genlijsten, waaronder GO-analyse (Yu et al., 2012). Een beperking van deze methode is dat er geen specifieke correctie wordt uitgevoerd voor gene-length bias. Bij RNA-seq data kan deze bias een rol spelen, omdat langere genen een grotere kans kunnen hebben om als differentieel tot expressie komend gen gedetecteerd te worden. Een alternatief hiervoor is `goseq`, dat specifiek is ontwikkeld voor GO-analyse van RNA-seq data waarbij rekening wordt gehouden met deze selectie-bias (Young et al., 2010). Omdat in deze analyse is gewerkt met gene symbols en geen genlengtes zijn meegenomen, is gekozen voor `clusterProfiler::enrichGO()` met alle betrouwbaar geteste genen als achtergrond. Deze beperking is meegenomen bij de interpretatie van de GO-resultaten.

Het volledige R-script waarmee de analyse is uitgevoerd staat in:

* [Transcriptomics RA script](Scripts/transcriptomics_RA_script_georganiseerd.R)

### Globale workflow

```mermaid
flowchart TD
    A[Raw paired-end FASTQ] --> B[Reads mappen met Rsubread]
    B --> C[BAM-bestanden]
    C --> D[Count matrix maken met featureCounts]
    D --> E[Differentiële genexpressie met DESeq2]
    E --> F[Annotatie met SYMBOL en ENTREZID]
    F --> G[Volcano plot]
    F --> H[KEGG/pathview analyse]
    F --> I[Gene Ontology analyse met clusterProfiler]
```

---

## Resultaten

Na het uitvoeren van de DESeq2-analyse zijn de RA-samples vergeleken met de controles. Hierbij was de controlegroep ingesteld als referentiegroep. Een positieve `log2FoldChange` betekent daarom dat een gen hoger tot expressie komt in RA ten opzichte van controle, terwijl een negatieve `log2FoldChange` betekent dat een gen lager tot expressie komt in RA.

Op basis van de criteria `padj < 0.05` en `|log2FoldChange| > 1` zijn in totaal **4528 differentieel tot expressie komende genen** gevonden. Hiervan kwamen **2057 genen hoger tot expressie in RA** en **2471 genen lager tot expressie in RA**. De volledige DESeq2-resultaten en de lijst met significante genen zijn opgeslagen in de map `Results/Tables/`.

Belangrijke resultaatbestanden zijn:

* [DESeq2 resultaten](Results/Tables/DESeq2_resultaten_RA_vs_control.csv)
* [Aantal differentieel tot expressie komende genen](Results/Tables/DEG_aantallen_RA_vs_control.csv)
* [Significante genen met Entrez ID](Results/Tables/significante_genen_RA_vs_control_met_EntrezID.csv)

### Volcano plot

De volcano plot laat per gen de `log2FoldChange` en aangepaste p-waarde zien. Genen die hoger tot expressie komen in RA liggen rechts van het midden. Genen die lager tot expressie komen in RA liggen links van het midden. De rode punten voldoen aan zowel de p-waardegrens als de fold change-grens.

![Volcano plot RA vs controle](Results/Figures/Volcanoplot_RA_vs_control.png)

In de volcano plot is te zien dat er veel genen significant verschillen tussen RA en controle. Sterk significant lager tot expressie komende genen zijn onder andere `ANKRD30BL`, `MT-ND6`, `RAB3IL1`, `SLC9A3R2` en `ZNF598`. Aan de rechterkant van de plot is onder andere `SRGN` zichtbaar als gen met hogere expressie in RA.

Daarnaast kwamen in de lijst met sterk hoger tot expressie komende genen meerdere immunoglobuline-gerelateerde genen naar voren, zoals `IGHV3-53`, `IGKV1-39`, `IGKV3D-15`, `IGHV6-1` en `IGHV1-69`. Dit wijst op betrokkenheid van B-cellen, antistoffen en adaptieve immuunprocessen in het RA-synovium.

### Gene Ontology analyse

De Gene Ontology analyse is uitgevoerd om te onderzoeken welke biologische processen, moleculaire functies en cellulaire componenten oververtegenwoordigd zijn in de lijst met significante genen. Hierbij zijn de significant veranderde genen vergeleken met alle betrouwbaar geteste genen als achtergrond.

In totaal zijn **138 significante GO-termen** gevonden voor Biological Process, **9 GO-termen** voor Molecular Function en **4 GO-termen** voor Cellular Component.

De GO Biological Process-resultaten laten vooral termen zien die te maken hebben met immuunactiviteit. Belangrijke termen zijn onder andere:

* `B cell mediated immunity`
* `adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains`
* `immunoglobulin mediated immune response`
* `lymphocyte differentiation`
* `lymphocyte mediated immunity`
* `immune response-regulating cell surface receptor signaling pathway`
* `antigen receptor-mediated signaling pathway`

De GO-dotplot laat zien dat vooral processen rond lymfocyten, B-cellen, antigeenreceptor-signalering, immunoglobuline-gemedieerde immuunrespons en adaptieve immuunrespons verrijkt zijn.

![GO Biological Process dotplot](Results/Figures/GO_BP_dotplot_RA_vs_control.png)

De GO Molecular Function-resultaten bevatten onder andere `antigen binding` en `immunoglobulin receptor binding`. De GO Cellular Component-resultaten bevatten vooral termen zoals `immunoglobulin complex`, `IgG immunoglobulin complex`, `immunoglobulin complex, circulating` en `plasma membrane signaling receptor complex`. Samen ondersteunen deze resultaten dat antistof-gerelateerde processen en immuuncelactiviteit sterk aanwezig zijn in de differentieel tot expressie komende genen.

De GO-resultaten zijn opgeslagen in:

* [GO Biological Process](Results/Tables/GO_Biological_Process_RA_vs_control.csv)
* [GO Molecular Function](Results/Tables/GO_Molecular_Function_RA_vs_control.csv)
* [GO Cellular Component](Results/Tables/GO_Cellular_Component_RA_vs_control.csv)

### KEGG-pathwayanalyse

Voor de KEGG-analyse is gekeken welke pathways overlap hebben met de significante genen. Daarnaast is de KEGG-pathway voor reumatoïde artritis (`hsa05323`) gevisualiseerd met `pathview`.

![KEGG rheumatoid arthritis pathway](Results/Pathways/hsa05323.RA_vs_control.png)

In de KEGG-visualisatie geeft rood een hogere expressie in RA aan en groen een lagere expressie in RA. In de RA-pathway zijn meerdere ontstekings- en immuungerelateerde genen hoger tot expressie in RA, waaronder `IFNG`, `CD80/86`, `MHCII`, `CD28`, `CTLA4`, `IL6`, `IL1B`, `IL8`, `CCL2`, `CCL3`, `CXCL1`, `CXCL5`, `TLR2/4` en `MMP1/3`. Dit wijst op verhoogde ontstekings- en immuunsignalering in de RA-samples.

De KEGG-overlapresultaten bevatten meerdere pathways die passen bij ontsteking en auto-immuniteit, waaronder:

* `Cytokine-cytokine receptor interaction`
* `Chemokine signaling pathway`
* `TNF signaling pathway`
* `NF-kappa B signaling pathway`
* `Th17 cell differentiation`
* `T cell receptor signaling pathway`
* `Toll-like receptor signaling pathway`
* `IL-17 signaling pathway`
* `Rheumatoid arthritis`
* `B cell receptor signaling pathway`
* `Antigen processing and presentation`

Deze resultaten moeten worden geïnterpreteerd als overlap met KEGG-pathways. Het bestand `KEGG_pathways_met_overlap.csv` bevat aantallen overlappende genen per pathway, maar dit is geen volledige statistische KEGG-enrichmentanalyse met p-waarden.

Het KEGG-resultaatbestand staat hier:

* [KEGG pathways met overlap](Results/Tables/KEGG_pathways_met_overlap.csv)

---

## Conclusie

Met behulp van RNA-seq analyse zijn duidelijke verschillen in genexpressie gevonden tussen synoviumbiopten van personen met reumatoïde artritis en controlepersonen. In totaal zijn **4528 differentieel tot expressie komende genen** gevonden, waarvan **2057 genen hoger tot expressie kwamen in RA** en **2471 genen lager tot expressie kwamen in RA**.

De biologische interpretatie van de resultaten wijst vooral richting immuunactiviteit. De GO-analyse liet verrijking zien van processen rond B-cellen, lymfocyten, antigeenreceptor-signalering, immunoglobuline-gemedieerde immuunrespons en adaptieve immuunrespons. Ook de Molecular Function- en Cellular Component-resultaten ondersteunen dit beeld, doordat termen zoals `antigen binding`, `immunoglobulin receptor binding` en `immunoglobulin complex` naar voren kwamen.

De KEGG-resultaten sluiten hierbij aan. Er was overlap met meerdere immuun- en ontstekingsgerelateerde pathways, zoals `cytokine-cytokine receptor interaction`, `chemokine signaling pathway`, `TNF signaling pathway`, `NF-kappa B signaling pathway`, `IL-17 signaling pathway` en de `rheumatoid arthritis` pathway. In de KEGG RA-pathway waren meerdere ontstekingsgerelateerde genen hoger tot expressie in RA, waaronder cytokines, chemokines en genen die betrokken zijn bij T-celactivatie en ontstekingssignalering.

Samen passen deze resultaten bij RA als auto-immuunziekte waarbij synovitis, immuunactivatie, B-celactiviteit, T-celactiviteit en ontstekingsprocessen een belangrijke rol spelen. Een beperking van deze analyse is dat er maar acht samples zijn gebruikt en dat er gewerkt is met subset40k FASTQ-bestanden. De resultaten zijn daarom vooral geschikt om de transcriptomics-workflow te demonstreren en moeten voorzichtig biologisch geïnterpreteerd worden.

Daarnaast is de GO-analyse uitgevoerd met `clusterProfiler::enrichGO()`. Deze methode is geschikt voor over-representation analysis van genlijsten, maar corrigeert niet specifiek voor gene-length bias. Omdat deze bias bij RNA-seq data invloed kan hebben op GO-resultaten, moeten de GO-termen voorzichtig worden geïnterpreteerd.

Vervolgonderzoek zou kunnen bestaan uit analyse van een grotere dataset, validatie van interessante genen met qPCR en een volledige statistische pathway-enrichmentanalyse. Ook zou een GO-analyse met `goseq` uitgevoerd kunnen worden om beter rekening te houden met gene-length bias. Ondanks deze beperkingen laat deze analyse zien hoe RNA-seq gebruikt kan worden om ziekte-gerelateerde genexpressiepatronen bij reumatoïde artritis te onderzoeken.

---

## Data stewardship en beheer

Voor de competentie beheren zijn aparte documenten toegevoegd waarin wordt uitgelegd hoe de projectgegevens en scripts zijn beheerd:

* [Data stewardship / beheren](Data_stewardship/Beheren.md)
* [GitHub beheer](Data_stewardship/GitHub_beheer.md)

In dit project is onderscheid gemaakt tussen ruwe data, tussenbestanden, verwerkte data, resultaten, scripts en documentatie. Grote bestanden zoals FASTQ-, BAM- en referentiebestanden zijn niet toegevoegd aan GitHub en worden uitgesloten via `.gitignore`. De mappen waarin deze bestanden lokaal horen te staan bevatten wel een README-bestand met uitleg.

---

## Gebruikte software en packages

De analyse is uitgevoerd in R/RStudio. De belangrijkste gebruikte packages zijn:

* `Rsubread`
* `Rsamtools`
* `DESeq2`
* `KEGGREST`
* `EnhancedVolcano`
* `pathview`
* `org.Hs.eg.db`
* `AnnotationDbi`
* `clusterProfiler`
* `enrichplot`
* `ggplot2`

De volledige R-sessie, RStudio-versie en gebruikte package-versies zijn opgeslagen in:

* [sessionInfo_transcriptomics_RA.txt](Results/sessionInfo_transcriptomics_RA.txt)

Deze informatie is toegevoegd om de analyse beter reproduceerbaar te maken.

---

## Bronnen

Gabriel, S. E. (2001). *The epidemiology of rheumatoid arthritis*. Rheumatic Disease Clinics of North America, 27(2), 269–281. doi:10.1016/S0889-857X(05)70201-5

Majithia, V., & Geraci, S. A. (2007). *Rheumatoid arthritis: Diagnosis and management*. The American Journal of Medicine, 120(11), 936–939. doi:10.1016/j.amjmed.2007.04.005

Platzer, A., Nussbaumer, T., Karonitsch, T., Smolen, J. S., & Aletaha, D. (2019). *Analysis of gene expression in rheumatoid arthritis and related conditions offers insights into sex-bias, gene biotypes and co-expression patterns*. PLOS ONE, 14(7), e0219698. doi:10.1371/journal.pone.0219698

Radu, A.-F., & Bungau, S. G. (2021). *Management of rheumatoid arthritis: An overview*. Cells, 10(11), 2857. doi:10.3390/cells10112857

Young, M. D., Wakefield, M. J., Smyth, G. K., & Oshlack, A. (2010). *Gene ontology analysis for RNA-seq: Accounting for selection bias*. Genome Biology, 11, R14. https://doi.org/10.1186/gb-2010-11-2-r14

Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). *ClusterProfiler: An R package for comparing biological themes among gene clusters*. OMICS: A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118
