# GitHub beheer

## Inleiding

Voor dit transcriptomics-project is GitHub gebruikt om de projectbestanden, scripts, resultaten en documentatie gestructureerd te beheren. De casus is niet uitgewerkt als traditioneel verslag, maar als GitHub-repository. Hierdoor is de analyse overzichtelijk te bekijken en kunnen bestanden op één centrale plek worden teruggevonden.

GitHub is in dit project gebruikt voor versiebeheer, projectorganisatie en reproduceerbaarheid.

## Repository aanmaken en koppelen

Eerst is op GitHub een repository aangemaakt met de naam `J2P4-Transcriptomics-RA`. Deze repository is gekoppeld aan GitHub Desktop. Via GitHub Desktop is de online repository lokaal op de computer gezet. Daardoor konden bestanden in Windows Verkenner worden toegevoegd en daarna via GitHub Desktop naar de online repository worden gepusht.

De repository is eerst private gehouden. Hierdoor kon de projectstructuur worden opgebouwd voordat de pagina eventueel openbaar werd gemaakt of ingeleverd.

## Structuur van de repository

De repository is opgebouwd met een vaste mappenstructuur:

```text
Data/
Mapped_reads/
Reference/
Reference_articles/
Results/
Scripts/
Documentation/
```

Deze structuur zorgt ervoor dat de verschillende onderdelen van het project gescheiden blijven. Scripts staan in `Scripts/`, resultaten in `Results/`, metadata en verwerkte data in `Data/`, en uitleg over databeheer in `Documentation/`.

De mappen `Data/Raw_FASTQ/`, `Mapped_reads/` en `Reference/` bevatten README-bestanden waarin staat welke lokale bestanden daar horen. Grote bestanden zoals FASTQ-, BAM- en referentiebestanden zijn niet toegevoegd aan GitHub.

## Gebruik van README-bestanden

README-bestanden zijn gebruikt om de repository begrijpelijker te maken. De hoofd-README geeft een overzicht van het project, de onderzoeksvraag, de methode, resultaten en conclusies. Daarnaast zijn er README-bestanden in specifieke mappen geplaatst om uit te leggen waarvoor die mappen bedoeld zijn.

Voorbeelden hiervan zijn:

```text
Data/Raw_FASTQ/README.md
Mapped_reads/README.md
Reference/README.md
Documentation/README.md
```

Hierdoor kan iemand die de repository opent sneller begrijpen welke bestanden aanwezig zijn en welke bestanden bewust niet zijn geüpload.

## Versiebeheer

Met GitHub Desktop zijn wijzigingen vastgelegd via commits. Bij een commit is een korte omschrijving toegevoegd, bijvoorbeeld:

```text
Projectbestanden en resultaten toegevoegd
Databeheerstructuur en bronnen toegevoegd
Beheerdocumentatie toegevoegd
```

Door commit messages te gebruiken is zichtbaar welke aanpassingen op welk moment zijn gedaan. Dit maakt het project beter controleerbaar. Als er later iets verkeerd gaat, is het ook makkelijker om terug te kijken welke bestanden zijn gewijzigd.

Na het maken van een commit is steeds `Push origin` gebruikt om de lokale wijzigingen online zichtbaar te maken op GitHub.

## Gebruik van .gitignore

Het bestand `.gitignore` is gebruikt om te bepalen welke bestanden niet naar GitHub mogen worden geüpload. Dit is belangrijk omdat sommige bestanden lokaal nodig zijn voor de analyse, maar niet geschikt zijn voor de repository.

Voorbeelden van bestanden die zijn uitgesloten:

```text
*.fastq
*.fq
*.fastq.gz
*.bam
*.bai
*.sam
*.fasta
*.fna
*.gtf
.RData
.Rhistory
```

Door `.gitignore` blijven grote sequencingbestanden, alignmentbestanden en referentiebestanden lokaal. Tegelijk blijft de repository overzichtelijk en beter geschikt om te delen.

## Reproduceerbaarheid

De GitHub-pagina draagt bij aan reproduceerbaarheid doordat de analysecode, metadata, resultaten en documentatie op een gestructureerde manier beschikbaar zijn. Iemand anders kan zien welke stappen zijn uitgevoerd, welke packages zijn gebruikt en waar de outputbestanden zijn opgeslagen.

De R-scripts bevatten commentaarregels waarin de analysestappen worden uitgelegd. Ook worden outputbestanden automatisch in de juiste mappen opgeslagen. Hierdoor is de workflow makkelijker te volgen en opnieuw uit te voeren.

## Terugvindbaarheid

Door de vaste mappenstructuur zijn bestanden makkelijk terug te vinden. Resultaattabellen staan bijvoorbeeld in `Results/Tables/`, figuren in `Results/Figures/` en pathway-output in `Results/Pathways/`. Bronnen staan in `Reference_articles/` en documentatie over databeheer staat in `Documentation/`.

Deze structuur maakt de repository bruikbaar voor iemand die het project wil beoordelen of de analyse opnieuw wil uitvoeren.

## Conclusie

GitHub is gebruikt als hulpmiddel om het transcriptomics-project overzichtelijk en reproduceerbaar te beheren. Door duidelijke mappen, README-bestanden, `.gitignore`, commits en gestructureerde scripts is inzichtelijk gemaakt welke data is gebruikt, welke stappen zijn uitgevoerd en waar de resultaten te vinden zijn. Dit ondersteunt zowel de inhoudelijke analyse als de competentie beheren.
