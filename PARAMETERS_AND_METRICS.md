# Paramètres et Métriques de la Pipeline ChIP-seq AR

Ce document décrit en détail chaque paramètre utilisé, sa justification, et les métriques de qualité attendues à chaque étape de la pipeline.

---

## Étape 0 — Téléchargement des données

### Outils et versions

| Outil | Commande | Rôle |
|-------|----------|------|
| SRA Toolkit | `prefetch`, `fasterq-dump` | Téléchargement et conversion SRA → FASTQ |
| wget | `wget` | Téléchargement du génome hg38 |
| Bowtie2 | `bowtie2-build` | Construction de l'index du génome |
| Samtools | `samtools faidx` | Indexation FASTA et extraction des tailles chromosomiques |

### Paramètres

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `--threads` | 8 | Accélère la conversion SRA → FASTQ. Ajustable selon la machine. |
| `--split-files` (PE) | activé si paired-end | Sépare les reads R1 et R2 dans deux fichiers distincts |
| Génome | hg38 (GRCh38) | Version la plus récente et la plus utilisée du génome humain |
| Source génome | UCSC `hgdownload.soe.ucsc.edu` | Inclut les noms de chromosomes au format UCSC (chr1, chr2...) compatibles avec MACS2 et HOMER |

### Fichiers générés

| Fichier | Description |
|---------|-------------|
| `{sample}.fastq.gz` | Reads bruts compressés |
| `hg38.fa` | Séquence du génome de référence (~3 Go) |
| `hg38_bt2.*.bt2` | 6 fichiers d'index Bowtie2 (~4 Go total) |
| `hg38.chrom.sizes` | Tailles des chromosomes (2 colonnes) |

---

## Étape 1 — Contrôle qualité (FastQC)

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| FastQC | `fastqc` | Évaluation de la qualité des reads individuels |
| MultiQC | `multiqc` | Agrégation de tous les rapports FastQC en un seul rapport |

### Paramètres

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-t` | 8 | Nombre de threads pour paralléliser l'analyse |
| `-o` | `results/fastqc/` | Répertoire de sortie des rapports |

### Métriques évaluées

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Score Phred moyen par base | > 28 | Qualité moyenne de chaque position dans le read. Un score < 20 indique un taux d'erreur > 1%. |
| Qualité par séquence | Pic principal > Q30 | Distribution des scores de qualité moyens par read |
| Contenu en adaptateurs | < 5% | Pourcentage de reads contenant des séquences d'adaptateurs Illumina |
| Taux de duplication | Variable | Attendu élevé en ChIP-seq (enrichissement naturel). Ne pas filtrer ici. |
| Contenu en GC | Proche de 40-42% | Distribution du contenu GC. Un décalage indique un biais de PCR ou une contamination. |
| Distribution de la longueur | Pic unique | Tous les reads devraient avoir la même longueur avant trimming |
| Séquences surreprésentées | < 1% | Indicateur de contamination ou de problèmes de librairie |

---

## Étape 2 — Trimming (Trim Galore)

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| Trim Galore | `trim_galore` | Wrapper autour de Cutadapt avec détection automatique des adaptateurs |
| Cutadapt | appelé en interne | Suppression des adaptateurs et des bases de basse qualité |

### Paramètres

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `--quality` | 20 | Supprime les bases en 3' avec un score Phred < 20 (taux d'erreur > 1%). La valeur 20 est le standard pour ChIP-seq. |
| `--length` | 20 | Élimine les reads de moins de 20 bp après trimming. Les reads courts s'alignent mal et génèrent du bruit. |
| `--fastqc` | activé | Lance automatiquement FastQC après trimming pour vérifier la qualité post-trim |
| `--cores` | 4 | Parallélisation. Trim Galore est limité à 4 cores maximum. |
| `--paired` (PE) | activé si paired-end | Traite les deux fichiers R1/R2 simultanément et synchronise les paires |
| Adaptateur | auto-détecté | Trim Galore détecte automatiquement le type d'adaptateur (Illumina TruSeq `AGATCGGAAGAGC`, Nextera, etc.) |

### Métriques post-trimming

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Reads restants | > 90% | Pourcentage de reads conservés après trimming |
| Bases trimmées | < 10% | Pourcentage de bases supprimées en 3' |
| Longueur moyenne post-trim | ≥ 25 bp | Longueur moyenne des reads après suppression des adaptateurs |

---

## Étape 3 — Alignement (Bowtie2)

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| Bowtie2 | `bowtie2` | Alignement des reads sur le génome hg38 |
| Samtools | `samtools view` | Conversion SAM → BAM en streaming (pipe) |

### Paramètres

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-x` | `hg38_bt2` | Chemin vers l'index Bowtie2 du génome |
| `--very-sensitive` | activé | Équivalent à `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50`. Maximise la sensibilité au prix de la vitesse. Recommandé pour ChIP-seq où on veut détecter tous les sites de liaison. |
| `--no-mixed` (PE) | activé | Interdit les alignements où un seul read de la paire s'aligne |
| `--no-discordant` (PE) | activé | Interdit les alignements discordants (orientation ou distance anormale) |
| `--maxins` (PE) | 500 | Taille maximale de l'insert en paired-end. 500 bp est approprié pour les librairies ChIP-seq standard. |
| `--threads` | 8 | Parallélisation de l'alignement |
| `-bS` (samtools) | activé | `-b` : sortie BAM ; `-S` : entrée SAM (via pipe) |

### Pourquoi Bowtie2 plutôt que BWA

Bowtie2 est préféré pour le ChIP-seq car il est optimisé pour les reads courts (50-150 bp), utilise moins de mémoire que BWA-MEM, et le mode `--very-sensitive` est bien validé pour les analyses de facteurs de transcription.

### Métriques d'alignement

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Taux d'alignement global | > 80% | Pourcentage de reads alignés au moins une fois. Un taux < 70% suggère une contamination ou un problème de librairie. |
| Alignements uniques | > 60% | Reads alignés à une seule position. Les multi-mappés sont ambigus. |
| Taux d'alignement concordant (PE) | > 70% | Paires alignées avec la bonne orientation et distance |

---

## Étape 4 — Post-traitement de l'alignement

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| Samtools | `sort`, `view`, `index`, `flagstat` | Tri, filtrage, indexation, statistiques |
| Picard | `MarkDuplicates` | Identification et suppression des duplicats de PCR |

### Paramètres de filtrage (samtools view)

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-q` | 10 | Qualité de mapping minimale (MAPQ ≥ 10). MAPQ 10 signifie une probabilité ≤ 10% que le read soit mal aligné. Élimine les alignements ambigus. |
| `-F` | 1804 | Flag binaire excluant : reads non mappés (4), read mate non mappé (8), alignement non primaire (256), échec QC (512), duplicat (1024). Total = 4+8+256+512+1024 = 1804. |

### Paramètres Picard MarkDuplicates

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `REMOVE_DUPLICATES` | true | Supprime physiquement les duplicats de PCR au lieu de les marquer. En ChIP-seq, les duplicats de PCR faussent le signal d'enrichissement. |
| `VALIDATION_STRINGENCY` | LENIENT | Tolère les BAM avec des en-têtes légèrement non-conformes. Évite les erreurs avec certains aligner. |

### Filtrage des chromosomes

| Chromosomes exclus | Raison |
|-------------------|--------|
| chrM | L'ADN mitochondrial est souvent surreprésenté (~30% des reads) car il est circulaire et très abondant |
| chrUn_* | Contigs non placés — génèrent des alignements ambigus |
| *_random | Contigs aléatoires non assignés à un chromosome |
| *_alt | Haplotypes alternatifs — créent des multi-mappings |
| *_hap | Haplotypes non-référence |
| *_fix | Corrections de la séquence de référence |

### Métriques post-filtrage

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Reads finaux (ChIP) | > 10 millions | Nombre minimum recommandé par ENCODE pour les facteurs de transcription |
| Reads finaux (Input) | > 10 millions | L'Input doit avoir au moins autant de reads que le ChIP |
| Taux de duplication | 10-50% | Attendu pour ChIP-seq. > 80% indique une librairie trop amplifiée. |

---

## Étape 5 — Contrôle qualité ChIP-seq (deepTools)

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| deepTools | `plotFingerprint` | Évalue l'enrichissement du signal ChIP vs Input |
| deepTools | `bamCoverage` | Génère des fichiers bigWig pour visualisation (IGV) |
| deepTools | `multiBamSummary` | Calcule la corrélation entre réplicats |
| deepTools | `plotCorrelation` | Visualise la matrice de corrélation |

### Paramètres plotFingerprint

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-b` | ChIP.bam Input.bam | Fichiers BAM à comparer |
| `--labels` | ChIP Input | Étiquettes pour la légende |
| `--numberOfProcessors` | 8 | Parallélisation |

### Paramètres bamCoverage

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `--normalizeUsing` | RPKM | Normalisation Reads Per Kilobase per Million. Permet de comparer la couverture entre échantillons de profondeurs différentes. |
| `--binSize` | 10 | Résolution de 10 bp pour le fichier bigWig. Bon compromis entre résolution et taille de fichier. |

### Métriques de qualité ChIP-seq

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Fingerprint plot | Courbe ChIP bien au-dessus de l'Input | Le ChIP doit montrer un enrichissement clair : la courbe doit s'écarter de la diagonale. L'Input doit suivre la diagonale (pas d'enrichissement). |
| Corrélation entre réplicats (Pearson) | > 0.8 | Mesure la reproductibilité. < 0.8 suggère une variabilité technique excessive. |
| NSC (Normalized Strand Coefficient) | > 1.05 | Mesure la qualité de l'enrichissement par cross-corrélation des brins. Valeurs ENCODE : > 1.05 acceptable, > 1.1 bon. |
| RSC (Relative Strand Coefficient) | > 0.8 | Ratio entre le pic de cross-corrélation au fragment size et le pic fantôme. > 0.8 acceptable, > 1.0 bon. |

---

## Étape 6 — Peak Calling (MACS2)

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| MACS2 | `callpeak` | Identification des régions de liaison AR enrichies par rapport au contrôle |
| BEDTools | `merge` | Fusion des régions chevauchantes |

### Paramètres MACS2

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-t` | ChIP BAM | Fichier(s) BAM de traitement (ChIP AR) |
| `-c` | Input BAM | Fichier BAM contrôle (Input DNA) |
| `-f` | BAM | Format du fichier d'entrée |
| `-g` | hs | Taille effective du génome humain (2.7×10⁹ bp). MACS2 utilise cette valeur pour estimer le bruit de fond. `hs` = human, `mm` = mouse, `ce` = C. elegans, `dm` = Drosophila. |
| `-q` | 0.01 | Seuil FDR (False Discovery Rate) de 1%. Chaque pic avec une q-value ≤ 0.01 a moins de 1% de chance d'être un faux positif. Plus strict que le défaut (0.05). |
| `--keep-dup` | 1 | Conserve au maximum 1 read à chaque position. Puisque nous avons déjà supprimé les duplicats avec Picard, ce paramètre est une sécurité supplémentaire. |
| `--call-summits` | activé | Identifie le point le plus enrichi (sommet) dans chaque pic. Crucial pour la recherche de motifs car le site de liaison est centré sur le sommet. |
| `--bdg` | activé | Génère des fichiers bedGraph de la couverture normalisée. Utile pour la visualisation dans IGV. |
| `-n` | AR_ChIP | Préfixe des fichiers de sortie |

### Filtrage post-peak calling

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Fold Enrichment | > 4 | Colonne 7 du fichier narrowPeak. FE = signal ChIP / signal Input. Un FE > 4 élimine les pics faiblement enrichis qui sont souvent des artefacts. |
| Fenêtre sommet | ± 150 bp | Crée une région de 300 bp centrée sur chaque sommet de pic. Cette taille est optimale pour la découverte de motifs de facteurs de transcription : assez large pour capturer le motif et le contexte, assez étroite pour minimiser le bruit. |

### Format du fichier narrowPeak (10 colonnes)

| Colonne | Nom | Description |
|---------|-----|-------------|
| 1 | chrom | Chromosome |
| 2 | chromStart | Début du pic |
| 3 | chromEnd | Fin du pic |
| 4 | name | Nom du pic |
| 5 | score | Score (0-1000) |
| 6 | strand | Brin (. = les deux) |
| 7 | signalValue | Fold Enrichment |
| 8 | pValue | -log10(p-value) |
| 9 | qValue | -log10(q-value) |
| 10 | peak | Position du sommet relative au début du pic |

### Métriques de peak calling

| Métrique | Seuil acceptable | Description |
|----------|-----------------|-------------|
| Nombre de pics totaux | 10 000 – 200 000 | Pour AR ChIP-seq dans LNCaP, on attend typiquement 20 000-80 000 pics bruts |
| Nombre de pics filtrés (FE > 4) | 5 000 – 50 000 | Pics de haute confiance après filtrage |
| FRiP (Fraction of Reads in Peaks) | > 1% | Pourcentage des reads ChIP tombant dans les pics. Standard ENCODE : > 1%. Typiquement 5-30% pour un bon ChIP-seq AR. |
| Fold Enrichment médian | > 4 | Enrichissement médian des pics. Un FE médian < 2 indique un ChIP-seq de mauvaise qualité. |

---

## Étape 7 — Annotation des pics

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| HOMER | `annotatePeaks.pl` | Annotation fonctionnelle des pics (position génomique, gène le plus proche) |
| ChIPseeker (R) | `annotatePeak()` | Annotation détaillée avec visualisation et enrichissement GO |
| clusterProfiler (R) | `enrichGO()`, `enrichKEGG()` | Analyse d'enrichissement fonctionnel des gènes cibles |

### Paramètres HOMER annotatePeaks

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Entrée | narrowPeak filtré | Pics de haute confiance uniquement |
| Génome | hg38 | Utilise les annotations HOMER pour hg38 |
| `-annStats` | fichier de sortie | Génère les statistiques de distribution des annotations |

### Paramètres ChIPseeker (R)

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| TxDb | `TxDb.Hsapiens.UCSC.hg38.knownGene` | Base de données de gènes UCSC pour hg38 |
| annoDb | `org.Hs.eg.db` | Mapping des identifiants de gènes humains |
| level | "gene" | Annotation au niveau du gène (pas du transcrit) |
| Promoter window | -3000 / +3000 bp du TSS | Fenêtre utilisée pour le profil de liaison autour du TSS |

### Paramètres enrichissement GO

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| ont | "BP" | Biological Process — le plus informatif pour les facteurs de transcription |
| pAdjustMethod | "BH" | Correction de Benjamini-Hochberg pour les tests multiples |
| qvalueCutoff | 0.05 | Seuil FDR de 5% pour les termes GO significatifs |
| readable | TRUE | Convertit les identifiants de gènes en symboles lisibles |

### Distribution attendue des pics AR

| Catégorie génomique | Pourcentage attendu | Explication |
|--------------------|--------------------:|-------------|
| Intronique | 30-45% | Les enhancers introniques sont des sites majeurs de liaison AR |
| Intergénique | 25-40% | Enhancers distaux, typiques des récepteurs nucléaires |
| Promoteur (≤ 1 kb du TSS) | 5-15% | Liaison directe aux promoteurs |
| Promoteur (1-3 kb du TSS) | 3-8% | Région promotrice étendue |
| Exonique | 1-3% | Rare pour les facteurs de transcription |
| 5'UTR / 3'UTR | 1-3% | Régions non traduites |

---

## Étape 8 — Découverte de motifs

### Outils

| Outil | Commande | Rôle |
|-------|----------|------|
| BEDTools | `getfasta` | Extraction des séquences FASTA sous les pics |
| HOMER | `findMotifsGenome.pl` | Découverte de motifs de novo et comparaison aux motifs connus |
| MEME-ChIP | `meme-chip` | Suite complète de découverte de motifs |
| TOMTOM | `tomtom` | Comparaison des motifs découverts aux bases de données |
| FIMO | `fimo` | Scan du génome pour les occurrences des motifs découverts |

### Paramètres HOMER findMotifsGenome

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Entrée | BED des sommets | Fichier BED des régions autour des sommets |
| Génome | hg38 | Génome de référence pour l'extraction des séquences |
| `-size` | 200 | Taille de la fenêtre d'analyse en bp centrée sur le sommet. 200 bp est le standard pour les facteurs de transcription (motif de ~10-20 bp + contexte). |
| `-mask` | activé | Masque les régions répétées (repeat-masked). Évite de trouver des motifs dans les éléments transposables qui n'ont rien à voir avec AR. |
| `-p` | 8 | Nombre de processeurs pour la parallélisation |

### HOMER : ce qui se passe en interne

HOMER effectue automatiquement :

1. Extraction des séquences sous les pics
2. Génération d'un background (séquences aléatoires avec même composition GC)
3. Recherche de novo des motifs enrichis (algorithme basé sur les binômes)
4. Comparaison avec sa base de données interne de motifs connus
5. Calcul des p-values et du pourcentage de cibles vs background

### Paramètres MEME-ChIP

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `-oc` | répertoire de sortie | Output directory (overwrite si existe) |
| `-maxw` | 20 | Largeur maximale du motif recherché. 20 bp couvre l'ARE complet (15 bp) avec de la marge. |
| `-minw` | 6 | Largeur minimale. 6 bp capture les demi-sites hexamériques (AGAACA). |
| `-meme-nmotifs` | 10 | Nombre maximum de motifs à découvrir. 10 est suffisant pour identifier l'ARE + les co-facteurs principaux. |
| `-meme-mod` | zoops | Zero or One Occurrence Per Sequence. Modèle approprié pour les facteurs de transcription qui se lient 0 ou 1 fois par région. |
| `-centrimo-local` | activé | Analyse l'enrichissement des motifs dans toutes les régions, pas seulement au centre. Utile car AR peut se lier légèrement décalé du sommet MACS2. |
| `-db` | JASPAR (si disponible) | Base de données de motifs connus pour la comparaison. JASPAR est la référence pour les motifs de vertébrés. |

### MEME-ChIP : composants exécutés automatiquement

| Composant | Fonction |
|-----------|----------|
| MEME | Découverte de novo par Expectation-Maximization. Trouve les motifs les plus statistiquement enrichis. |
| STREME | Découverte de motifs courts enrichis (ancien DREME). Complémentaire à MEME pour les motifs courts. |
| CentriMo | Analyse de l'enrichissement central : vérifie si le motif est enrichi au centre des pics (attendu pour le motif primaire). |
| TOMTOM | Compare les motifs découverts aux bases de données connues (JASPAR, HOCOMOCO). |
| FIMO | Scanne chaque séquence pour les occurrences des motifs. |
| SpaMo | Analyse de l'espacement entre paires de motifs (co-facteurs). |

### Sous-échantillonnage

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Seuil | 10 000 séquences | Si plus de 10 000 pics, on ne garde que les 10 000 meilleurs (par score) pour la recherche de motifs. MEME est O(n²) en complexité, donc trop de séquences ralentit énormément sans améliorer les résultats. |
| Critère de tri | Score du pic (colonne 5) | Les pics les plus significatifs ont le signal le plus fort et contiennent plus probablement le motif ARE. |

### Paramètres TOMTOM

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Entrée | motifs découverts (format MEME) | Motifs de novo identifiés par MEME |
| Base de données | JASPAR vertebrates | Base de référence pour les motifs de facteurs de transcription de vertébrés |

### Paramètres FIMO

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| `--thresh` | 1e-4 | Seuil de p-value pour rapporter une occurrence. 1e-4 est un bon compromis entre sensibilité et spécificité. |
| Entrée | motifs MEME + génome hg38 | Scan génomique complet pour toutes les occurrences des motifs découverts |

---

## Étape 9 — Analyse et visualisation des motifs (R)

### Outils

| Outil | Fonction | Rôle |
|-------|----------|------|
| universalmotif | `read_meme()`, `compare_motifs()` | Lecture des motifs MEME et comparaison avec l'ARE de référence |
| ggseqlogo | `ggseqlogo()` | Génération des logos de séquence (visualisation des motifs) |
| ggplot2 | Graphiques | Visualisations générales |

### Paramètres de comparaison ARE

| Paramètre | Valeur | Justification |
|-----------|--------|---------------|
| Méthode de comparaison | PCC (Pearson Correlation Coefficient) | Mesure la similarité entre la matrice de poids (PWM) du motif découvert et l'ARE de référence. PCC va de -1 à 1, avec 1 = identique. |
| ARE de référence | AGAACAnnnTGTTCT | Matrice de poids construite à partir du consensus ARE canonique, avec 85% de probabilité pour la base consensus et 5% pour chaque autre base. Les 3 positions centrales (nnn) ont une probabilité uniforme de 25% par base. |

### Matrice de poids de l'ARE canonique

```
Position :  A     G     A     A     C     A     n     n     n     T     G     T     T     C     T
A        : 0.85  0.05  0.85  0.85  0.05  0.85  0.25  0.25  0.25  0.05  0.05  0.05  0.05  0.05  0.05
C        : 0.05  0.05  0.05  0.05  0.85  0.05  0.25  0.25  0.25  0.05  0.05  0.05  0.05  0.85  0.05
G        : 0.05  0.85  0.05  0.05  0.05  0.05  0.25  0.25  0.25  0.05  0.85  0.05  0.05  0.05  0.05
T        : 0.05  0.05  0.05  0.05  0.05  0.05  0.25  0.25  0.25  0.85  0.05  0.85  0.85  0.05  0.85
```

### Fichiers de sortie

| Fichier | Description |
|---------|-------------|
| `motif_summary.csv` | Tableau récapitulatif : consensus, largeur, E-value de chaque motif |
| `motif_N_logo.pdf` | Logo de séquence individuel pour chaque motif |
| `all_motifs.pdf` | Tous les logos sur une seule figure |
| `ARE_reference_vs_discovered.pdf` | Comparaison visuelle ARE canonique vs motif principal découvert |
| `homer_top_known.csv` | Top 10 motifs connus enrichis selon HOMER |

---

## Rapport HTML automatique

Le rapport `report.html` généré à la fin contient :

| Section | Contenu |
|---------|---------|
| Configuration | Tous les paramètres d'entrée utilisés |
| Peak Calling | Nombre de pics totaux et filtrés |
| Motif ARE attendu | Séquence consensus de référence |
| Motifs découverts | Tableau avec consensus et largeur |
| Timing | Durée de chaque étape |
| Fichiers clés | Liens vers tous les fichiers de résultats |
| ChIP-seq QC | Image du fingerprint plot |

---

## Résumé des seuils critiques

| Étape | Métrique | Seuil | Conséquence si non respecté |
|-------|----------|-------|-----------------------------|
| 1 | Phred moyen | > 28 | Augmenter la stringence du trimming |
| 2 | Reads conservés | > 90% | Vérifier la qualité de la librairie |
| 3 | Taux d'alignement | > 80% | Vérifier contamination ou mauvaise librairie |
| 4 | Reads finaux | > 10M | Augmenter la profondeur de séquençage |
| 4 | Taux de duplication | < 80% | Re-préparer la librairie (sur-amplification) |
| 5 | Fingerprint | Courbe décollée de la diagonale | Faible enrichissement ChIP |
| 5 | Corrélation réplicats | > 0.8 | Mauvaise reproductibilité |
| 6 | Nombre de pics | > 5 000 | Ajuster le seuil q-value ou vérifier le ChIP |
| 6 | FRiP | > 1% | Faible enrichissement |
| 8 | Motif ARE dans le top 3 | p-value < 1e-50 | Vérifier la qualité des données ou le design expérimental |
