# Start ####

  # Working directory instellen
setwd("C:/Users/31613/OneDrive - NHL Stenden/00 jaar 2/periode 4/Transcriptomics/Github/J2P4-Transcriptomics-RA/")
getwd()

  # Mappen aanmaken
    # Data           = inputdata en verwerkte data
      # /Metadata    = sample-informatie en metadata
      # /Processed   = officiële countmatrix + zelfgemaakte NCBI countmatrix + featureCounts-output
      # /raw_fastq   = ruwe subset FASTQ-bestanden
    # Mapped_reads   = BAM, sorted BAM, BAI, summary/indel-bestanden
    # Reference      = NCBI FASTA, GTF en Rsubread-index
    # Results        = analyse-uitkomsten en sessionInfo
      # /Tables      = DESeq2, significante genen, KEGG, GO tabellen
      # /Figures     = volcano plot en GO dotplot
      # /Pathways    = KEGG/pathview figuren
    # Scripts        = R-scripts


folders <- c(
  "Data/Metadata",
  "Data/Processed",
  "Data/raw_fastq",
  "Reference",
  "Mapped_reads",
  "Results/Tables",
  "Results/Figures",
  "Results/Pathways",
  "Scripts"
)

for (folder in folders) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}




  # Packages laden [1e deel] ####
library(BiocManager)   # install.packages('BiocManager')
library(Rsubread)      # BiocManager::install('Rsubread')
library(Rsamtools)     # BiocManager::install("Rsamtools")






  # Referentiegenoom en annotatie: NCBI RefSeq GRCh38.p14 ####
    # In deze analyse wordt het humane NCBI RefSeq referentiegenoom gebruikt:
    # GCF_000001405.40_GRCh38.p14

  # Belangrijk:
    # - FASTA en GTF moeten van dezelfde bron komen.
    # - Daarom worden hier beide bestanden van NCBI gebruikt.
    # https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

  # Controleren of de bestanden aanwezig zijn 
list.files("Reference")

  # Bestandsnamen instellen 
ncbi_fasta <- "Reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
ncbi_gtf   <- "Reference/genomic.gtf"

  # Controleren of R de bestanden kan vinden
file.exists(ncbi_fasta)
file.exists(ncbi_gtf)

stopifnot(file.exists(ncbi_fasta))
stopifnot(file.exists(ncbi_gtf))



  # Referentiegenoom indexeren ####
    # De index wordt gemaakt met Rsubread::buildindex.
    # indexSplit = TRUE wordt gebruikt omdat het humane genoom groot is.
    # De basename is de naam waaronder de indexbestanden worden opgeslagen.

if (!file.exists("Reference/ref_NCBI_GRCh38p14.files")) {
  buildindex(
    basename = "Reference/ref_NCBI_GRCh38p14",
    reference = ncbi_fasta,
    memory = 4000,
    indexSplit = TRUE
  )
}

  # Controleren of de indexbestanden zijn aangemaakt
list.files("Reference", pattern = "ref_NCBI_GRCh38p14")

  # Controleren of de raw data bestanden aanwezig zijn
list.files("Data/raw_fastq")






  # Sample-informatie ####
    # In deze tabel staat per sample:
      # - de originele SRR-code
      # - de naam die gebruikt wordt voor BAM-bestanden en kolomnamen
      # - de conditie: control of RA

sample_info <- data.frame(
  run = c(
    "SRR4785819",
    "SRR4785820",
    "SRR4785828",
    "SRR4785831",
    "SRR4785979",
    "SRR4785980",
    "SRR4785986",
    "SRR4785988"
  ),
  sample_name = c(
    "SRR4785819_con1",
    "SRR4785820_con2",
    "SRR4785828_con3",
    "SRR4785831_con4",
    "SRR4785979_ra1",
    "SRR4785980_ra2",
    "SRR4785986_ra3",
    "SRR4785988_ra4"
  ),
  condition = c(
    "control",
    "control",
    "control",
    "control",
    "RA",
    "RA",
    "RA",
    "RA"
  )
)

sample_info

  # Sample informatie opslaan
write.csv(
  sample_info,
  "Data/Metadata/sample_info_RA_subset40k.csv",
  row.names = FALSE
)



  # FASTQ-bestanden instellen ####
    # De bestanden staan in Data/raw_fastq.
    # Per sample is er een forward read bestand (_1_) en reverse read bestand (_2_).

readfile1 <- file.path(
  "Data",
  "raw_fastq",
  paste0(sample_info$run, "_1_subset40k.fastq")
)

readfile2 <- file.path(
  "Data",
  "raw_fastq",
  paste0(sample_info$run, "_2_subset40k.fastq")
)

  # Controleren of alle FASTQ-bestanden aanwezig zijn
file.exists(readfile1)
file.exists(readfile2)

all(file.exists(readfile1))
all(file.exists(readfile2))

  # Script stoppen als niet alle FASTQ-bestanden gevonden worden
stopifnot(all(file.exists(readfile1)))
stopifnot(all(file.exists(readfile2)))




  # Mapping ####
    # - Onderstaande code mapt de controle en RA samples tegen het humane GRCh38.p14 referentiegenoom.
    # - index = "Reference/ref_NCBI_GRCh38p14": gebruik de eerder gemaakte index van het referentiegenoom.
    # - Let hierbij op dat je alleen de 'basename' moet gebruiken, omdat bij het indexeren meerdere files zijn gemaakt.
    # - readfile1 = ...: het bestand met forward reads (_1_) voor dit monster.
    # - readfile2 = ...: het bestand met reverse reads (_2_) voor dit monster.
    # - input_format = "FASTQ": de reads zijn FASTQ-bestanden.
    # - type = "rna": de reads zijn afkomstig van RNA-seq.
    # - output_file = ...: de naam van het resultaatbestand (.BAM)

align_results <- lapply(seq_len(nrow(sample_info)), function(i) {
  
  output_bam <- file.path(
    "Mapped_reads",
    paste0(sample_info$sample_name[i], ".BAM")
  )
  
  if (!file.exists(output_bam)) {
    
    align(
      index = "Reference/ref_NCBI_GRCh38p14",
      readfile1 = readfile1[i],
      readfile2 = readfile2[i],
      input_format = "FASTQ",
      type = "rna",
      output_file = output_bam
    )
    
  } else {
    message("BAM-bestand bestaat al en wordt overgeslagen: ", output_bam)
  }
  
})

names(align_results) <- sample_info$sample_name

  # Controleren of de BAM-bestanden zijn aangemaakt
list.files("Mapped_reads", pattern = "\\.BAM$")


  # BAM-bestanden sorteren & indexeren ####
    # - In een BAM-bestand staat: welke reads gemapt zijn, waar op het genoom ze gemapt zijn, de mappingkwaliteit en eventuele mismatches of clipping.
    # - Sorteren: zorgt ervoor dat de reads op volgorde van genoompositie staan.
    # - Indexeren: maakt een soort inhoudsopgave zodat programma's snel naar een specifieke regio kunnen springen.
    # - De gesorteerde en geïndexeerde BAM-bestanden kunnen daarna worden bekeken in bijvoorbeeld IGV.


bam_files <- file.path(
  "Mapped_reads",
  paste0(sample_info$sample_name, ".BAM")
)

  # Controleren of alle BAM-bestanden bestaan
file.exists(bam_files)
stopifnot(all(file.exists(bam_files)))


  # BAM-bestanden sorteren
lapply(seq_along(bam_files), function(i) {
  
  sorted_bam <- file.path(
    "Mapped_reads",
    paste0(sample_info$sample_name[i], ".sorted.bam")
  )
  
  if (!file.exists(sorted_bam)) {
    
    sortBam(
      file = bam_files[i],
      destination = file.path(
        "Mapped_reads",
        paste0(sample_info$sample_name[i], ".sorted")
      )
    )
    
  } else {
    message("Gesorteerde BAM bestaat al en wordt overgeslagen: ", sorted_bam)
  }
  
})

  # Gesorteerde BAM-bestanden
sorted_bams <- file.path(
  "Mapped_reads",
  paste0(sample_info$sample_name, ".sorted.bam")
)

  # Controleren of alle gesorteerde BAM-bestanden bestaan
file.exists(sorted_bams)
stopifnot(all(file.exists(sorted_bams)))


  # Gesorteerde BAM-bestanden indexeren
lapply(sorted_bams, function(bam) {
  
  bai_file <- paste0(bam, ".bai")
  
  if (!file.exists(bai_file)) {
    indexBam(bam)
  } else {
    message("BAM-index bestaat al en wordt overgeslagen: ", bai_file)
  }
  
})

  # Controleren of de gesorteerde en geïndexeerde bestanden zijn aangemaakt
list.files("Mapped_reads", pattern = "sorted")

length(list.files("Mapped_reads", pattern = "sorted.bam$"))
length(list.files("Mapped_reads", pattern = "sorted.bam.bai$"))





  # (Count matrix) Uitgelijnde reads (uit de BAM-bestanden) vergelijken met de genannotaties in het GTF-bestand ####
    # - Deze vergelijking bepaalt met de functie featureCounts hoeveel reads er op elk gen vallen.
    # - files = ...: de BAM-bestanden die je wilt tellen.
    # - annot.ext = ...: het GTF-bestand dat de humane genannotaties bevat.
    # - isPairedEnd = TRUE: de reads zijn paired-end, dus elk sample heeft een _1_ en _2_ bestand.
    # - isGTFAnnotationFile = TRUE: het annotatiebestand is een GTF-bestand.
    # - GTF.featureType = "exon": reads worden geteld op exon-regio's.
    # - GTF.attrType = "gene_id": exonen worden gekoppeld aan genen via gene_id.
    # - useMetaFeatures = TRUE: counts worden samengevat op genniveau.

  # Controleren of alle gesorteerde BAM-bestanden en het GTF-bestand aanwezig zijn
file.exists(sorted_bams)
file.exists(ncbi_gtf)

  # Script stoppen als bestanden ontbreken
stopifnot(all(file.exists(sorted_bams)))
stopifnot(file.exists(ncbi_gtf))

count_matrix <- featureCounts(
  files = sorted_bams,
  annot.ext = ncbi_gtf,
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

  # Count matrix bekijken
str(count_matrix)

  # Alleen de tellingen ophalen
counts <- count_matrix$counts
head(counts)

  # Kolomnamen overzichtelijk maken
colnames(counts) <- sample_info$sample_name

head(counts)

  # Extra controles
dim(counts)
colSums(counts)
sum(rowSums(counts) > 0)


  # Zelfgemaakte NCBI subset40k count matrix
write.csv(
  counts,
  "Data/Processed/RA_countmatrix_NCBI_GRCh38p14_subset40k.csv"
)

  # FeatureCounts annotatie opslaan
write.csv(
  count_matrix$annotation,
  "Data/Processed/RA_countmatrix_NCBI_GRCh38p14_subset40k_annotation.csv",
  row.names = FALSE
)

  # Volledig featureCounts object opslaan
saveRDS(
  count_matrix,
  "Data/Processed/featureCounts_NCBI_GRCh38p14_subset40k.rds"
)

# Count matrix inladen ####
   # Voor de differentiële genexpressieanalyse wordt de officiële count matrix gebruikt:
   # Data/Processed/count_matrix_RA.txt
    # De zelfgemaakte NCBI count matrix is opgeslagen als extra output van de mapping/featureCounts-stap.

counts_raw <- read.delim(
  "Data/Processed/count_matrix_RA.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

  # Bekijken hoe het bestand eruitziet
dim(counts_raw)
head(counts_raw)
colnames(counts_raw)

  # Definitieve count matrix maken 
counts <- counts_raw

  # Omzetten naar matrix met integers
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

  # Controles
dim(counts)
head(counts)
colSums(counts)
sum(rowSums(counts) > 0)

  # Count matrix opslaan als CSV voor overzicht
write.csv(counts, "Data/Processed/RA_countmatrix.csv")


  # Packages laden [2e deel] ####
    # Let op: installeren hoeft maar één keer.
    # Als de packages al geïnstalleerd zijn, kun je de BiocManager::install regels overslaan.

library(DESeq2)            # BiocManager::install("DESeq2")
library(KEGGREST)          # BiocManager::install("KEGGREST")
library(EnhancedVolcano)   # BiocManager::install("EnhancedVolcano")
library(pathview)          # BiocManager::install("pathview")
library(org.Hs.eg.db)      # BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)     # BiocManager::install("AnnotationDbi")





# Metadata ####
    # Tabel die beschrijft welke monsters controle zijn en welke afkomstig zijn van RA-patiënten.
    # De volgorde van de metadata moet overeenkomen met de kolommen van de count matrix.

      # SRR4785819, SRR4785820, SRR4785828 en SRR4785831 = controle
      # SRR4785979, SRR4785980, SRR4785986 en SRR4785988 = RA
sample_table <- data.frame(
  condition = c("control", "control", "control", "control",
                "RA", "RA", "RA", "RA")
)

  # Rijnamen gelijk maken aan de kolomnamen van de count matrix
rownames(sample_table) <- colnames(counts)

  # Zorgt ervoor dat control de referentiegroep is.
  # Hierdoor betekent een positieve log2FoldChange: hogere expressie in RA t.o.v. controle.
sample_table$condition <- factor(
  sample_table$condition,
  levels = c("control", "RA")
)

  # Controleren of de kolomnamen van counts overeenkomen met de rijnamen van sample_table
all(colnames(counts) == rownames(sample_table))
      # !!! MOET TRUE GEVEN !!!
stopifnot(all(colnames(counts) == rownames(sample_table)))

sample_table

  # Metadata opslaan
sample_metadata_export <- data.frame(
  sample_id = c(
    "SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831",
    "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988"
  ),
  condition = c(
    "control", "control", "control", "control",
    "RA", "RA", "RA", "RA"
  ),
  diagnosis = c(
    "Normal", "Normal", "Normal", "Normal",
    "Rheumatoid arthritis established",
    "Rheumatoid arthritis established",
    "Rheumatoid arthritis established",
    "Rheumatoid arthritis established"
  ),
  age = c(31, 15, 31, 42, 54, 66, 60, 59),
  sex = rep("female", 8)
)

write.csv(
  sample_metadata_export,
  file = "Data/Metadata/sample_metadata_RA.csv",
  row.names = FALSE
)






# Statistiek ####
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_table,
  design = ~ condition
)

  # Genen met heel lage read counts verwijderen
dds <- dds[rowSums(counts(dds)) > 10, ]

  # DESeq2 uitvoeren
dds <- DESeq(dds)

  # Resultaten: RA vergeleken met controle
resultaten <- results(
  dds,
  contrast = c("condition", "RA", "control")
)

  # Resultaten bekijken & Opslaan
head(resultaten)
                     # Rijnamen: gene symbols.
                     # baseMean: Gemiddelde expressie van het gen over alle samples.
                     # log2FoldChange: verandering in expressie in RA t.o.v. controle.
                     # lfcSE: standaardfout van de log2FoldChange.
                     # stat: Wald-statistiek.
                     # pvalue: p-waarde.
                     # padj: gecorrigeerde p-waarde volgens Benjamini-Hochberg.
summary(resultaten, alpha = 0.05)            
                                

resultaten_df <- as.data.frame(resultaten)
resultaten_df$GeneSymbol <- rownames(resultaten_df)

write.csv(
  resultaten_df,
  "Results/Tables/DESeq2_resultaten_RA_vs_control.csv",
  row.names = FALSE
)






# Aantal differentieel tot expressie komende genen ####
  # Criteria:
              # padj < 0.05 = significant na correctie voor multiple testing
              # log2FoldChange > 1 = hoger tot expressie in RA
              # log2FoldChange < -1 = lager tot expressie in RA

up_RA <- sum(resultaten_df$padj < 0.05 & resultaten_df$log2FoldChange > 1, na.rm = TRUE)
down_RA <- sum(resultaten_df$padj < 0.05 & resultaten_df$log2FoldChange < -1, na.rm = TRUE)
totaal_DEG <- up_RA + down_RA

  # Aantallen overzichtelijk tonen
cat("Aantal hoger in RA:", up_RA, "\n")
cat("Aantal lager in RA:", down_RA, "\n")
cat("Totaal aantal DEG's:", totaal_DEG, "\n")

  # Samenvatting met aantallen opslaan
deg_summary <- data.frame(
  categorie = c("Hoger in RA", "Lager in RA", "Totaal DEG"),
  aantal = c(up_RA, down_RA, totaal_DEG)
)

write.csv(
  deg_summary,
  file = "Results/Tables/DEG_aantallen_RA_vs_control.csv",
  row.names = FALSE
)




# Resultaten sorteren ####
  # Voor het sorteren worden genen met NA in padj eerst verwijderd.
  # Genen met padj = NA zijn door DESeq2 niet betrouwbaar getest,
  # bijvoorbeeld door outliers of te lage counts.

resultaten_zonder_NA <- resultaten_df[
  !is.na(resultaten_df$padj),
]

  # hoogste_fold_change: betrouwbaar geteste genen met de sterkste opregulatie in RA bovenaan.
  # laagste_fold_change: betrouwbaar geteste genen met de sterkste neerregulatie in RA bovenaan.
  # laagste_p_waarde: genen met de laagste gecorrigeerde p-waarde bovenaan.

hoogste_fold_change <- resultaten_zonder_NA[order(resultaten_zonder_NA$log2FoldChange, decreasing = TRUE),]

laagste_fold_change <- resultaten_zonder_NA[order(resultaten_zonder_NA$log2FoldChange, decreasing = FALSE),]

laagste_p_waarde <- resultaten_zonder_NA[order(resultaten_zonder_NA$padj, decreasing = FALSE),]

head(hoogste_fold_change)
head(laagste_fold_change)
head(laagste_p_waarde)

  # Gesorteerde resultaten opslaan
write.csv(hoogste_fold_change,"Results/Tables/hoogste_fold_change_RA.csv", row.names = FALSE)

write.csv(laagste_fold_change,"Results/Tables/laagste_fold_change_RA.csv", row.names = FALSE)

write.csv(laagste_p_waarde,"Results/Tables/laagste_p_waarde_RA.csv", row.names = FALSE)






# Significante genen selecteren ####
   # Dit zijn genen met padj < 0.05 en minimaal een verdubbeling of halvering in expressie.

significante_genen <- resultaten_df[
  !is.na(resultaten_df$padj) &
    resultaten_df$padj < 0.05 &
    abs(resultaten_df$log2FoldChange) > 1,
]

significante_genen <- significante_genen[
  order(significante_genen$padj, decreasing = FALSE),
]

head(significante_genen)
nrow(significante_genen)

write.csv(
  significante_genen,
  file = "Results/Tables/significante_genen_RA_vs_control.csv",
  row.names = FALSE
)




# Volcano plot ####
      # Een volcano plot laat per gen de log2FoldChange en aangepaste p-waarde zien.
      # Op de x-as staat de log2FoldChange:
      #   log2FoldChange > 0 = hoger tot expressie in RA
      #   log2FoldChange < 0 = lager tot expressie in RA
      # Op de y-as staat de -log10 van de aangepaste p-waarde.
      # Hoger op de y-as betekent sterker significant.

  # Genen met NA in padj of log2FoldChange verwijderen voor de plot
volcano_data <- resultaten_df[
  !is.na(resultaten_df$padj) &
    !is.na(resultaten_df$log2FoldChange),
]

  # Top 10 meest significante genen labelen
top_labels <- volcano_data$GeneSymbol[
  order(volcano_data$padj)
][1:10]

  # Volcano plot maken en opslaan
png(
  filename = "Results/Figures/Volcanoplot_RA_vs_control.png",
  width = 8,
  height = 10,
  units = "in",
  res = 500
)

EnhancedVolcano(
  volcano_data,
  lab = volcano_data$GeneSymbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = top_labels,
  title = "RA vs controle",
  subtitle = "Differentiële genexpressie",
  xlab = "log2 fold change",
  ylab = "-log10 aangepaste p-waarde",
  caption = "Criteria: padj < 0.05 en |log2FC| > 1"
)

dev.off()

  # Controleren en bekijken van Volcano plot
file.exists("Results/Figures/Volcanoplot_RA_vs_control.png")
browseURL("Results/Figures/Volcanoplot_RA_vs_control.png")


# KEGG Pathway-analyse ####
    # Voor pathview moeten gene symbols worden omgezet naar ENTREZID.
    # De log2FoldChange wordt gebruikt om genen in KEGG pathways rood/groen te kleuren.
    # Rood = hogere expressie in RA
    # Groen = lagere expressie in RA


  # Gene symbols omzetten naar ENTREZID
resultaten_df$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = resultaten_df$GeneSymbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

  # Alleen genen bewaren met ENTREZID en log2FoldChange
res_kegg <- resultaten_df[
  !is.na(resultaten_df$ENTREZID) &
    !is.na(resultaten_df$log2FoldChange),
]

  # Controleren hoeveel genen gemapt konden worden
nrow(res_kegg)
head(res_kegg)

  # Opslaan
write.csv(
  res_kegg,
  file = "Results/Tables/DESeq2_resultaten_met_EntrezID.csv",
  row.names = FALSE
)



# Vector maken voor pathview ####
    # De namen van de vector zijn ENTREZID's.
    # De waarden zijn de log2FoldChanges.

gene_fc_temp <- tapply(
  res_kegg$log2FoldChange,
  res_kegg$ENTREZID,
  mean,
  na.rm = TRUE
)

  # Omzetten naar gewone named numeric vector
gene_fc <- as.numeric(gene_fc_temp)
names(gene_fc) <- names(gene_fc_temp)

  # Missende of lege waarden verwijderen
gene_fc <- gene_fc[!is.na(gene_fc)]
gene_fc <- gene_fc[!is.na(names(gene_fc))]
gene_fc <- gene_fc[names(gene_fc) != ""]

  # Controleren
head(gene_fc)
length(gene_fc)
str(gene_fc)






# KEGG pathway visualiseren: Rheumatoid arthritis ####
    # hsa05323 = Rheumatoid arthritis pathway bij Homo sapiens

pathview(
  gene.data = gene_fc,
  pathway.id = "05323",
  species = "hsa",
  gene.idtype = "entrez",
  out.suffix = "RA_vs_control",
  kegg.dir = "Results/Pathways",
  limit = list(gene = 5)
)

# Controleren of de pathwayfiguur in Results/Pathways staat & openen
file.exists("Results/Pathways/hsa05323.RA_vs_control.png")
browseURL("Results/Pathways/hsa05323.RA_vs_control.png")




# KEGG pathways zoeken met overlap met significante genen ####

  # ENTREZID toevoegen aan significante genen
significante_genen$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = significante_genen$GeneSymbol,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

  # Significante genen mét EntrezID opslaan
write.csv(
  significante_genen,
  file = "Results/Tables/significante_genen_RA_vs_control_met_EntrezID.csv",
  row.names = FALSE
)

  # Alleen significante genen met ENTREZID bewaren
sig_entrez <- unique(na.omit(significante_genen$ENTREZID))

length(sig_entrez)

  # KEGG pathways ophalen waarin deze genen voorkomen
     # In blokjes om problemen met te lange queries te voorkomen
chunks <- split(sig_entrez, ceiling(seq_along(sig_entrez) / 100))

gene_pathways <- unlist(
  lapply(chunks, function(x) {
    keggLink("pathway", paste0("hsa:", x))
  })
)

  # Controleren of er pathways gevonden zijn
length(gene_pathways)
head(gene_pathways)

  # Aantal overlappende genen per pathway tellen
pathway_counts <- sort(table(gene_pathways), decreasing = TRUE)

  # Namen van humane KEGG pathways ophalen
pathway_info <- keggList("pathway", "hsa")

pathway_info_df <- data.frame(
  pathway_id = sub("^path:", "", names(pathway_info)),
  pathway_name = as.character(pathway_info)
)

  # Resultaten combineren
kegg_results <- data.frame(
  pathway_id = sub("^path:", "", names(pathway_counts)),
  count = as.numeric(pathway_counts)
)

kegg_results <- merge(
  kegg_results,
  pathway_info_df,
  by = "pathway_id",
  all.x = TRUE
)

kegg_results <- kegg_results[
  order(kegg_results$count, decreasing = TRUE),
]

head(kegg_results, 20)

  # Opslaan
write.csv(
  kegg_results,
  file = "Results/Tables/KEGG_pathways_met_overlap.csv",
  row.names = FALSE
)


# GO-analyse met goseq - packages [3e deel] ####

  # Doel:
    # Onderzoeken welke Gene Ontology-termen oververtegenwoordigd zijn in de significant differentieel tot expressie komende genen.
    # Voor RNA-seq data wordt hier goseq gebruikt, omdat goseq rekening kan houden met gene-length bias. 
    # Langere genen hebben bij RNA-seq namelijk een grotere kans om als differentieel tot expressie komend gen gedetecteerd te worden.


  # Packages laden

library(goseq)            # BiocManager::install("goseq")
library(rtracklayer)      # BiocManager::install("rtracklayer")
library(GenomicRanges)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)            # BiocManager::install("GO.db")
library(ggplot2)



# Controle: benodigde objecten aanwezig ####
    # Deze objecten worden eerder in het script gemaakt.
    # Als één van deze FALSE geeft, moet het script eerst vanaf DESeq2 opnieuw gerund worden.

exists("resultaten_df")
exists("resultaten_zonder_NA")
exists("significante_genen")

stopifnot(exists("resultaten_df"))
stopifnot(exists("resultaten_zonder_NA"))
stopifnot(exists("significante_genen"))




# Genlengtes bepalen uit GTF-bestand ####
    # goseq heeft lengte-informatie nodig om te corrigeren voor gene-length bias.
    # Hiervoor wordt hetzelfde GTF-bestand gebruikt als bij het mappen/tellen van reads.
    # Pas de bestandsnaam hieronder aan als jouw GTF-bestand anders heet.

gtf_file <- "Reference/genomic.gtf"

file.exists(gtf_file)
  # !!! MOET TRUE GEVEN !!!
stopifnot(file.exists(gtf_file))

  # GTF-bestand inladen
gtf <- rtracklayer::import(gtf_file)

  # Alleen exonen gebruiken voor genlengtes
exons <- gtf[gtf$type == "exon"]

  # Controleren welke kolommen met annotatie aanwezig zijn
colnames(mcols(exons))

  # Gennaam/symbool uit de GTF halen
    # Bij NCBI RefSeq GTF staat het gensymbool meestal in de kolom 'gene'.
if ("gene" %in% colnames(mcols(exons))) {
  gene_symbols_gtf <- mcols(exons)$gene
} else if ("gene_name" %in% colnames(mcols(exons))) {
  gene_symbols_gtf <- mcols(exons)$gene_name
} else if ("gene_id" %in% colnames(mcols(exons))) {
  gene_symbols_gtf <- mcols(exons)$gene_id
} else {
  stop("Geen geschikte gennaam-kolom gevonden in de GTF. Controleer colnames(mcols(exons)).")
}

  # Gennaam toevoegen aan exon-object
exons$gene_symbol <- as.character(gene_symbols_gtf)

  # Exonen zonder gennaam verwijderen
exons <- exons[
  !is.na(exons$gene_symbol) &
    exons$gene_symbol != ""
]

  # Exonen per gen groeperen
exons_per_gene <- split(exons, exons$gene_symbol)

  # Overlappende exonen per gen samenvoegen en totale exonlengte per gen berekenen
gene_lengths <- sum(width(reduce(exons_per_gene)))

  # Omzetten naar gewone numerieke vector
gene_lengths <- as.numeric(gene_lengths)
names(gene_lengths) <- names(exons_per_gene)

  # Controleren
head(gene_lengths)
length(gene_lengths)
summary(gene_lengths)

  # Controle of genen uit DESeq2 ook in de lengtevector zitten
sum(resultaten_zonder_NA$GeneSymbol %in% names(gene_lengths)) 
sum(significante_genen$GeneSymbol %in% names(gene_lengths))






# Binaire genenlijst maken voor goseq ####
    # goseq verwacht een named vector met:
      # 1 = significant differentieel tot expressie komend gen
      # 0 = niet significant, maar wel betrouwbaar getest
      # De achtergrond bestaat uit alle betrouwbaar geteste genen: genen zonder NA in padj.

tested_genes <- unique(resultaten_zonder_NA$GeneSymbol)

sig_genes <- unique(significante_genen$GeneSymbol)

  # Alleen genen gebruiken waarvoor ook een genlengte beschikbaar is
tested_genes_with_length <- intersect(tested_genes, names(gene_lengths))

length(tested_genes)
length(tested_genes_with_length)

  # Binaire vector maken
de_genes <- as.integer(tested_genes_with_length %in% sig_genes)
names(de_genes) <- tested_genes_with_length

  # Lengtevector in dezelfde volgorde zetten
bias_data <- gene_lengths[names(de_genes)]

  # Controles
table(de_genes)
length(de_genes)
length(bias_data)
all(names(de_genes) == names(bias_data))
  # !!! MOET TRUE GEVEN !!!
stopifnot(all(names(de_genes) == names(bias_data)))





# Probability Weighting Function berekenen en opslaan ####
    # Deze plot laat zien of genlengte invloed heeft op de kans dat een gen als differentieel tot expressie komend wordt gevonden.

png(
  filename = "Results/Figures/goseq_PWF_RA_vs_control.png",
  width = 8,
  height = 7,
  units = "in",
  res = 300
)

pwf <- nullp(
  de_genes,
  bias.data = bias_data,
  plot.fit = TRUE
)

dev.off()

  # Controleren en openen
file.exists("Results/Figures/goseq_PWF_RA_vs_control.png")
browseURL("Results/Figures/goseq_PWF_RA_vs_control.png")

head(pwf)
summary(pwf$bias.data)





# GO-annotatie ophalen ####
    # De genen in deze analyse zijn gene symbols.
    # Met org.Hs.eg.db worden GO-termen gekoppeld aan deze gene symbols.

gene2go_all <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = names(de_genes),
  keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
)

  # Dubbele regels en NA's verwijderen
gene2go_all <- unique(gene2go_all)

gene2go_all <- gene2go_all[
  !is.na(gene2go_all$GOALL) &
    !is.na(gene2go_all$ONTOLOGYALL),
]

  # Controles 
head(gene2go_all) 
table(gene2go_all$ONTOLOGYALL) 
length(unique(gene2go_all$SYMBOL))


# Functie om goseq per GO-ontologie uit te voeren ####

run_goseq_ontology <- function(ontology_code) {
  
  # GO-termen selecteren voor BP, MF of CC
  gene2go_ont <- gene2go_all[
    gene2go_all$ONTOLOGYALL == ontology_code,
    c("SYMBOL", "GOALL")
  ]
  
  # Omzetten naar gene2cat-formaat:
  # per gen een lijst met bijbehorende GO-termen
  gene2cat <- split(gene2go_ont$GOALL, gene2go_ont$SYMBOL)
  
  # Alleen genen bewaren die ook in de goseq-genenlijst zitten
  gene2cat <- gene2cat[names(gene2cat) %in% names(de_genes)]
  
  # goseq uitvoeren
  go_result <- goseq(
    pwf,
    gene2cat = gene2cat,
    method = "Wallenius"
  )
  
  # GO-termen en ontologie toevoegen
  go_result$term <- Term(GOTERM[go_result$category])
  go_result$ontology <- Ontology(GOTERM[go_result$category])
  
  # BH-correctie uitvoeren op de over-representation p-waarde
  go_result$padj <- p.adjust(
    go_result$over_represented_pvalue,
    method = "BH"
  )
  
  # GeneRatio toevoegen: aantal significante genen in GO-term / totaal aantal significante genen
  go_result$GeneRatio <- go_result$numDEInCat / sum(de_genes == 1)
  
  # Sorteren op aangepaste p-waarde
  go_result <- go_result[
    order(go_result$padj, decreasing = FALSE),
  ]
  
  return(go_result)
}






# GO-analyse uitvoeren per ontologie ####

go_bp_all <- run_goseq_ontology("BP")
go_mf_all <- run_goseq_ontology("MF")
go_cc_all <- run_goseq_ontology("CC")

  # Significante GO-termen selecteren
go_bp_sig <- go_bp_all[
  !is.na(go_bp_all$padj) &
    go_bp_all$padj < 0.05,
]

go_mf_sig <- go_mf_all[
  !is.na(go_mf_all$padj) &
    go_mf_all$padj < 0.05,
]

go_cc_sig <- go_cc_all[
  !is.na(go_cc_all$padj) &
    go_cc_all$padj < 0.05,
]

  # Resultaten bekijken
head(go_bp_sig)
head(go_mf_sig)
head(go_cc_sig)

  # Aantal significante GO-termen
nrow(go_bp_sig)
nrow(go_mf_sig)
nrow(go_cc_sig)




# GO-resultaten opslaan ####
    # Dezelfde bestandsnamen worden gebruikt als eerder, zodat de links in de README blijven werken.

write.csv(
  go_bp_sig,
  "Results/Tables/GO_Biological_Process_RA_vs_control.csv",
  row.names = FALSE
)

write.csv(
  go_mf_sig,
  "Results/Tables/GO_Molecular_Function_RA_vs_control.csv",
  row.names = FALSE
)

write.csv(
  go_cc_sig,
  "Results/Tables/GO_Cellular_Component_RA_vs_control.csv",
  row.names = FALSE
)

  # Optioneel: ook alle GO-resultaten opslaan, inclusief niet-significante termen
write.csv(
  go_bp_all,
  "Results/Tables/GO_Biological_Process_RA_vs_control_all_goseq.csv",
  row.names = FALSE
)

write.csv(
  go_mf_all,
  "Results/Tables/GO_Molecular_Function_RA_vs_control_all_goseq.csv",
  row.names = FALSE
)

write.csv(
  go_cc_all,
  "Results/Tables/GO_Cellular_Component_RA_vs_control_all_goseq.csv",
  row.names = FALSE
)

  # Controleren of bestanden zijn opgeslagen 
file.exists("Results/Tables/GO_Biological_Process_RA_vs_control.csv") 
file.exists("Results/Tables/GO_Molecular_Function_RA_vs_control.csv") 
file.exists("Results/Tables/GO_Cellular_Component_RA_vs_control.csv")


# GO Biological Process dotplot maken ####
    # Dotplot van de top 8 GO Biological Process-termen op basis van padj.

go_bp_plot_df <- head(go_bp_sig, 8)

if (nrow(go_bp_plot_df) > 0) {
  
  # Lange termen afbreken zodat de figuur leesbaar blijft
  go_bp_plot_df$term_wrapped <- sapply(
    go_bp_plot_df$term,
    function(x) paste(strwrap(x, width = 45), collapse = "\n")
  )
  
  go_bp_plot <- ggplot(
    go_bp_plot_df,
    aes(
      x = GeneRatio,
      y = reorder(term_wrapped, GeneRatio),
      size = numDEInCat,
      color = padj
    )
  ) +
    geom_point() +
    labs(
      title = "GO Biological Process - RA vs controle",
      x = "GeneRatio",
      y = NULL,
      size = "Count",
      color = "padj"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Plot bekijken in RStudio
  go_bp_plot
  
  # Plot opslaan
  ggsave(
    filename = "Results/Figures/GO_BP_dotplot_RA_vs_control.png",
    plot = go_bp_plot,
    width = 14,
    height = 7,
    dpi = 300
  )
  
  # Controleren en openen
  file.exists("Results/Figures/GO_BP_dotplot_RA_vs_control.png")
  browseURL("Results/Figures/GO_BP_dotplot_RA_vs_control.png")
  
} else {
  message("Er zijn geen significante GO Biological Process-termen gevonden met padj < 0.05.")
}

# RStudio versie & Package versies ####
capture.output(
  {
    cat("RStudio versie:\n")
    print(RStudio.Version()$long_version)
    
    cat("\nR versie en sessie-informatie:\n")
    print(sessionInfo())
    
    cat("\nBelangrijkste package versies:\n")
    cat("BiocManager:", as.character(packageVersion("BiocManager")), "\n")
    cat("Rsubread:", as.character(packageVersion("Rsubread")), "\n")
    cat("Rsamtools:", as.character(packageVersion("Rsamtools")), "\n")
    cat("DESeq2:", as.character(packageVersion("DESeq2")), "\n")
    cat("KEGGREST:", as.character(packageVersion("KEGGREST")), "\n")
    cat("EnhancedVolcano:", as.character(packageVersion("EnhancedVolcano")), "\n")
    cat("pathview:", as.character(packageVersion("pathview")), "\n")
    cat("org.Hs.eg.db:", as.character(packageVersion("org.Hs.eg.db")), "\n")
    cat("AnnotationDbi:", as.character(packageVersion("AnnotationDbi")), "\n")
    cat("goseq:", as.character(packageVersion("goseq")), "\n")
    cat("rtracklayer:", as.character(packageVersion("rtracklayer")), "\n")
    cat("GenomicRanges:", as.character(packageVersion("GenomicRanges")), "\n")
    cat("GO.db:", as.character(packageVersion("GO.db")), "\n")
    cat("ggplot2:", as.character(packageVersion("ggplot2")), "\n")
  },
  file = "Results/sessionInfo_transcriptomics_RA.txt"
)

file.exists("Results/sessionInfo_transcriptomics_RA.txt")



# Kleine DESeq2-previewtabel maken voor GitHub ####
# Dit blok kan los worden uitgevoerd als de volledige DESeq2-resultaten al bestaan.

setwd("C:/Users/31613/OneDrive - NHL Stenden/00 jaar 2/periode 4/Transcriptomics/Github/J2P4-Transcriptomics-RA/")

# Volledige DESeq2-resultaten opnieuw inladen
deseq2_resultaten <- read.csv(
  "Results/Tables/DESeq2_resultaten_RA_vs_control.csv",
  header = TRUE,
  check.names = FALSE
)

# Genen zonder padj verwijderen en sorteren op laagste aangepaste p-waarde
deseq2_top25 <- deseq2_resultaten[
  !is.na(deseq2_resultaten$padj),
]

deseq2_top25 <- deseq2_top25[
  order(deseq2_top25$padj, decreasing = FALSE),
]

# Top 25 bewaren
deseq2_top25 <- head(deseq2_top25, 25)

# Opslaan als klein bestand dat GitHub wel makkelijk kan tonen
write.csv(
  deseq2_top25,
  file = "Results/Tables/DESeq2_top25_laagste_padj_RA_vs_control.csv",
  row.names = FALSE
)

# Controleren
file.exists("Results/Tables/DESeq2_top25_laagste_padj_RA_vs_control.csv")