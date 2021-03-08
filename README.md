# CrispritzWebApp

Versione corrente: crisprme_off.py

Seguire le istruzioni del paragrafo INSTALLATION (Phase 1) da [CRISPRitz](https://github.com/InfOmics/CRISPRitz) per installare conda. Seguire i primi 8 punti, fino a quando (base) appare nel terminale.

Modificare il file environment_for_python_3_8.yml, sostituendo name: base, con un nome da dare all'environment.

Modificare il file environment_for_python_3_8.yml, sostituendo prefix: /root/miniconda3 con la giusta directory dove conda è installato.

Eseguire

```
conda env create -f environment_for_python_3_8.yml
```

Per creare l'environment di conda dove poter eseguire CRISPRme.

Attivare il nuovo environment con

```
conda activate nome_environment
```

Scaricare CrispritzWebApp tramite questo github. Per poter funzionare, è necessario scaricare ulteriori files:

- Genoma: scaricare questo [zip](https://www.dropbox.com/s/01j6vg6dc75wkn0/genomes.zip?dl=0) e spostare le cartelle hg38_ref e hg38_ref+hg38_1000genomeproject all'interno della cartella Genomes
- Indici: scaricare questo [zip](https://www.dropbox.com/s/wd297qosnl82xto/genome_lib.zip?dl=0) e spostare le cartelle NGG_2_hg38_ref e NGG_2_hg38_ref+hg38_1000genomeproject all'interno della cartella genome_library
- Dizionari: scaricare questo [zip](https://www.dropbox.com/s/g2pe8tig7g6oj9c/dict.zip?dl=0), creare una cartella dal nome dictionaries , e spostare entrambi i file .json all'interno di questa cartella. NB la cartella dictionaries deve essere situata nello stesso livello della cartella CrispritzWebApp

Da terminale, spostarsi all'interno della cartella CrispritzWebApp, ed eseguire:

```
python crisprme_off.py
```

Per attivare CRISPRme

Usando un browser, la pagina principale si trova su:

```
127.0.0.1:8080
```

Esempio di pagina wait job (31/03/2020):

```
127.0.0.1:8080/load?job=Q47PXDTBC8
```

Esempio di result summary:

```
127.0.0.1:8080/result?job=Q47PXDTBC8
```
