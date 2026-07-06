# GitHub beheer

Voor dit transcriptomics-project is GitHub gebruikt om scripts, metadata, resultaten en documentatie gestructureerd te beheren. De casus is niet als traditioneel verslag ingeleverd, maar als GitHub-repository. Hierdoor kan de beoordelaar de workflow, resultaten en data stewardship-documentatie op één plek terugvinden.

## Repository en structuur

De repository `J2P4-Transcriptomics-RA` is lokaal gekoppeld aan GitHub Desktop. Bestanden zijn lokaal toegevoegd en daarna met commits en `Push origin` online gezet.

De repository is opgebouwd met een vaste structuur:

```text
Data/
Mapped_reads/
Reference/
Reference_articles/
Results/
Scripts/
Data_stewardship/
```

Scripts staan in `Scripts/`, resultaten in `Results/`, metadata en verwerkte data in `Data/`, bronnen in `Reference_articles/` en uitleg over databeheer in `Data_stewardship/`.

## README-bestanden

De hoofd-README geeft een overzicht van het project, de methode, resultaten, conclusie, gebruikte bronnen en data stewardship. In meerdere submappen staan extra README-bestanden. Deze leggen uit waarvoor de map bedoeld is en welke bestanden wel of niet op GitHub staan.

Voorbeelden zijn:

```text
Data/raw_fastq/README.md
Mapped_reads/README.md
Reference/README.md
Data_stewardship/README.md
```

## Versiebeheer

Wijzigingen zijn vastgelegd met commits in GitHub Desktop. Bij elke commit is kort beschreven wat er is aangepast, bijvoorbeeld het toevoegen van resultaten, het verbeteren van de README of het bijwerken van bronverwijzingen. Hierdoor is terug te zien hoe het project is opgebouwd en welke bestanden zijn gewijzigd.

## .gitignore

Het `.gitignore`-bestand is gebruikt om grote of lokale bestanden uit de repository te houden, zoals FASTQ-, BAM-, BAI-, FASTA-, GTF-, index- en RStudio-bestanden. Daardoor blijft GitHub overzichtelijk en worden alleen relevante scripts, metadata, resultaten en documentatie gedeeld.

## Reproduceerbaarheid en terugvindbaarheid

De repository ondersteunt reproduceerbaarheid doordat de analysecode, metadata, resultaten, bronnen en `sessionInfo` gescheiden en duidelijk gelinkt zijn. Resultaattabellen staan in `Results/Tables/`, figuren in `Results/Figures/`, pathway-output in `Results/Pathways/` en de gebruikte R-code in `Scripts/`.

Door deze structuur is duidelijk welke data is gebruikt, welke stappen zijn uitgevoerd en waar de resultaten terug te vinden zijn.
