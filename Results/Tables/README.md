# Resultaattabellen

Deze map bevat de tabellen die zijn gegenereerd tijdens de RNA-seq analyse van RA-synovium versus controlesynovium.

## DESeq2-resultaten

* `DESeq2_resultaten_RA_vs_control.csv`
  Volledige DESeq2-output voor RA versus controle.

* `DESeq2_resultaten_met_EntrezID.csv`
  DESeq2-resultaten aangevuld met Entrez ID’s voor functionele analyses.

* `DESeq2_top25_laagste_padj_RA_vs_control.csv`
  Compacte tabel met de 25 genen met de laagste aangepaste p-waarde. Deze tabel is toegevoegd omdat de volledige DESeq2-tabel te groot is om makkelijk in GitHub te bekijken.

* `DEG_aantallen_RA_vs_control.csv`
  Samenvatting van het aantal significant differentieel tot expressie komende genen.

* `significante_genen_RA_vs_control.csv`
  Alle significante genen op basis van `padj < 0.05` en `|log2FoldChange| > 1`.

* `significante_genen_RA_vs_control_met_EntrezID.csv`
  Significante genen aangevuld met Entrez ID’s.

## Toplijsten

* `hoogste_fold_change_RA.csv`
  Genen met de hoogste positieve log2 fold change, dus hogere expressie in RA.

* `laagste_fold_change_RA.csv`
  Genen met de meest negatieve log2 fold change, dus lagere expressie in RA.

* `laagste_p_waarde_RA.csv`
  Genen met de laagste aangepaste p-waarden (`padj`).

## Gene Ontology-resultaten

* `GO_Biological_Process_RA_vs_control.csv`
  Significante GO Biological Process-termen.

* `GO_Molecular_Function_RA_vs_control.csv`
  Significante GO Molecular Function-termen.

* `GO_Cellular_Component_RA_vs_control.csv`
  Significante GO Cellular Component-termen.

* `GO_Biological_Process_RA_vs_control_all_goseq.csv`
  Volledige goseq-output voor Biological Process.

* `GO_Molecular_Function_RA_vs_control_all_goseq.csv`
  Volledige goseq-output voor Molecular Function.

* `GO_Cellular_Component_RA_vs_control_all_goseq.csv`
  Volledige goseq-output voor Cellular Component.

## KEGG-resultaten

* `KEGG_pathways_met_overlap.csv`
  Overzicht van KEGG-pathways die overlap vertonen met de significant veranderde genen.

## Opmerking

De belangrijkste resultaten zijn samengevat in de hoofd-README. Deze map bevat de onderliggende tabellen waarmee de resultaten gecontroleerd kunnen worden.
