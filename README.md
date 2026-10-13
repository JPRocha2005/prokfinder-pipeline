# Prokfinder

Pipeline para **busca, download e curadoria de genomas de procariotos** (Bacteria e Archaea) do NCBI, combinando
os metadados das montagens (NCBI), os dados fenotípicos das linhagens (BacDive) e a taxonomia filogenômica
(GTDB).

O objetivo é montar um banco de genomas **de alta qualidade, classificados taxonomicamente e sem redundância**,
selecionados por critérios biológicos. Por exemplo: todos os procariotos com temperatura ótima de crescimento
entre 45 e 80 °C isolados de fontes termais.

---

## Funcionalidades

| Etapa | Módulo | O que faz |
|---|---|---|
| 0 | `download-reference-db-ncbi.sh` | Baixa os metadados de todas as montagens de Bacteria/Archaea (NCBI Datasets, pacote *dehydrated*: só o `fetch.txt`, sem os FASTA). Une as linhas repetidas da mesma montagem, formata o cabeçalho (`_` no lugar de espaços), troca valores vazios ou inválidos por `NULL`, adiciona a linhagem NCBI (nome e TaxID de cada rank, do domínio à espécie, e o TaxID da cepa ou `NO_TAXID`) com o taxonkit e separa GenBank (`GCA_`) de RefSeq (`GCF_`) |
| 0 | `download-reference-db-bacdive.sh` | Baixa todas as linhagens do BacDive pela API e extrai 29 traços fenotípicos (temperatura, pH, salinidade, O₂, morfologia, metabolismo...) para um TSV |
| 0 | `download-reference-db-gtdb.sh` | Baixa a taxonomia e os metadados do GTDB (`ar53`/`bac120`, última release), confere o MD5 e separa a linhagem em uma coluna por rank |
| 1 | `search-genomes-in-reference-db.sh` | Busca os genomas que satisfazem os filtros do usuário (NCBI + BacDive + GTDB) |
| 2 | `download-selected-genomes-from-db.sh` | Baixa só os genomas selecionados (`datasets rehydrate`), apenas a sequência genômica (`*_genomic.fna`, nunca `cds_from_genomic`/`rna_from_genomic`) |
| 3 | `filter-low-quality-contigs-from-genomes.sh` | Remove contigs < 1000 pb, ordena por tamanho e renomeia (`Contig_1`, `Contig_2`, ...) com o seqkit |
| 4 | `evaluate-and-filter-genome-by-qc.sh` | seqkit (tamanho, nº de contigs, N50, bases N, GC) + CheckM2 (completude, contaminação) e filtragem por qualidade |
| 5 | `assign-taxonomy-gtdb.sh` | Taxonomia GTDB pela tabela do banco; os não encontrados vão para `NON-CLASSIFIED-GENOMES.txt` ou, com `--gtdbtk`, são classificados pelo GTDB-Tk |
| 6 | `build-taxonomic-dir-hierarchy.sh` | Hierarquia de diretórios `Domínio/Filo/Classe/Ordem/Família/Gênero/Espécie` |
| 7 | `dereplicate-genomes-skder.sh` | Desreplicação com o skDER dentro de cada táxon (diretório-folha) |
| 8 | `summarize-final-results.sh` | Lista final de genomas, tabelas de estatística e gráficos |

Critérios de qualidade padrão (etapa 4): **nº de contigs < 1000**, **N50 ≥ 5 kb**, **bases N < 100 000**,
**completude CheckM2 ≥ 75 %** e **contaminação CheckM2 ≤ 5 %**.

Cada etapa acrescenta suas colunas a uma **única tabela de metadados** (`OUTPUT/GENOME-METADATA/genome-metadata.tsv`,
uma linha por genoma). Ao final ela tem os metadados do NCBI, os traços do BacDive, o status do download, as
métricas de qualidade, a taxonomia GTDB, o diretório taxonômico e a indicação de genoma representativo.

---

## Estrutura do repositório

```
prokfinder-scripts/
├── main-prokfinder.sh          # controla a execução: opções, dependências, diretórios e ordem dos módulos
├── lib/
│   └── common.sh               # funções compartilhadas (log, tabelas, caminhos do banco de referência)
├── modules/                    # um script por etapa; cada um também pode ser executado sozinho (-h)
│   ├── download-reference-db-ncbi.sh
│   ├── download-reference-db-bacdive.sh
│   ├── download-reference-db-gtdb.sh
│   ├── search-genomes-in-reference-db.sh
│   ├── download-selected-genomes-from-db.sh
│   ├── filter-low-quality-contigs-from-genomes.sh
│   ├── evaluate-and-filter-genome-by-qc.sh
│   ├── assign-taxonomy-gtdb.sh
│   ├── build-taxonomic-dir-hierarchy.sh
│   ├── dereplicate-genomes-skder.sh
│   ├── summarize-final-results.sh
│   └── plot-statistics.R
└── README.md
```

---

## Dependências

| Programa | Usado em | Instalação (conda) |
|---|---|---|
| bash ≥ 4, awk, gawk, curl, unzip, gzip | todas as etapas | sistema |
| `datasets`, `dataformat` (NCBI Datasets) | etapas 0 e 2 | `conda install -c conda-forge ncbi-datasets-cli` |
| `taxonkit` | etapa 0 (NCBI) | `conda install -c bioconda taxonkit` |
| `jq` | etapa 0 (BacDive) | `conda install -c conda-forge jq` |
| `seqkit` | etapas 3 e 4 | `conda install -c bioconda seqkit` |
| `checkm2` + banco de dados | etapa 4 | `conda install -c bioconda checkm2` e `checkm2 database --download` |
| `skder` | etapa 7 | `conda install -c bioconda skder` |
| `gtdbtk` + pacote de referência (opcional) | etapa 5 com `--gtdbtk` | `conda install -c bioconda gtdbtk` |
| `Rscript` (opcional, só R base) | gráficos da etapa 8 | `conda install -c conda-forge r-base` |

Variáveis de ambiente:

| Variável | Obrigatória | Descrição |
|---|---|---|
| `CHECKM2DB` | sim (exceto com `-n`) | caminho do **arquivo** do banco do CheckM2 (`uniref100.KO.1.dmnd`) |
| `GTDBTK_DATA_PATH` | só com `--gtdbtk` | pasta do pacote de referência do GTDB-Tk |
| `NCBI_API_KEY` | não | chave da API do NCBI (aumenta o limite de 5 para 10 requisições/s) |

O `main-prokfinder.sh` confere todas as dependências necessárias para a execução pedida **antes** de começar e lista
todos os problemas de uma vez.

---

## Uso

```bash
export CHECKM2DB=/caminho/CheckM2_database/uniref100.KO.1.dmnd

# Termófilos (temperatura ótima 45-80 °C), em bastonete ou coco, apenas genomas de isolados
./main-prokfinder.sh -T 45:80 -O rod,coccus -g isolate_genome -p 32

# Anaeróbios Gram-negativos de uma lista de países
./main-prokfinder.sh -B gram_stain=negative -M anaerobe -C paises.txt

# Só selecionar (sem download): gênero Thermus no GTDB com completude CheckM2 ≥ 95 % segundo o GTDB
./main-prokfinder.sh --gtdb-taxon g__Thermus --gtdb-range checkm2_completeness=95: -n
```

Na primeira execução os bancos de referência (NCBI, BacDive e GTDB) são baixados para
`Prokfinder/INPUT/REFERENCE-METADATA`. Isso leva horas e ocupa alguns GB; as execuções seguintes reaproveitam esses
dados. A pasta de dados pode ser escolhida com `-d` (padrão: `./Prokfinder`) e o log completo fica em
`Prokfinder/prokfinder.log`.

### Filtros de seleção

Pelo menos um filtro é obrigatório. Valores separados por vírgula e as linhas de um arquivo de lista são combinados
com **OU**; filtros diferentes são combinados com **E**.

| Opção | Banco | Descrição |
|---|---|---|
| `-t, --taxid` / `-f, --taxid-file` | NCBI | TaxID do organismo **ou de qualquer rank da linhagem** (ex.: `1239` seleciona todos os Bacillota) |
| `-a, --accession` / `-A, --accession-file` | NCBI | código de acesso GCA_/GCF_ (versão opcional; GCA_x.1 também encontra o par GCF_x.1) |
| `--ncbi-taxon` | NCBI | nome de um táxon da linhagem NCBI, qualquer rank |
| `-i` / `-I` | NCBI | texto contido na fonte de isolamento (ex.: `"hot spring"`) |
| `-c` / `-C` | NCBI | país/local da localização geográfica |
| `-g` | NCBI | `isolate_genome`, `mag` ou `all` |
| `-G` | NCBI | faixa de GC (%) da montagem |
| `--ncbi-range COLUNA=min:max` / `--ncbi-match COLUNA=texto` | NCBI | qualquer outra coluna dos metadados do NCBI |
| `-s` | NCBI | tabela pesquisada: `genbank` (padrão) ou `refseq` |
| `-T` / `-P` / `-S` | BacDive | temperatura, pH e NaCl **ótimos** (`min:max`, um dos lados pode ficar vazio) |
| `--temp-growth` / `--ph-growth` / `--salt-growth` | BacDive | qualquer valor em que houve crescimento (mínimo, máximo, ótimo...) |
| `-O` / `-M` / `-B traço=classe` | BacDive | formato celular, tolerância a O₂ e qualquer traço dividido em classes (Gram, motilidade, esporulação, fonte de energia/carbono, doador de elétrons...) |
| `--bacdive-range coluna=min:max` | BacDive | qualquer traço numérico do BacDive |
| `--gtdb-taxon` | GTDB | nome de um táxon da linhagem GTDB, qualquer rank (`g__Thermus` ou `Thermus`) |
| `--gtdb-range` / `--gtdb-match` | GTDB | qualquer coluna dos metadados do GTDB (ex.: `checkm2_completeness=95:`) |

`./main-prokfinder.sh -h` mostra todas as opções e as classes disponíveis para cada traço.

**Como funcionam os filtros do BacDive:** uma linhagem pode ter vários valores para o mesmo traço, vindos de
referências diferentes. Na tabela do BacDive, `;` separa referências e `:` separa valores da mesma referência (ex.:
`temperature_opt = 37;45`). Uma linhagem é aceita quando **pelo menos um** dos seus valores satisfaz o filtro, e
linhagens sem valor para o traço filtrado são excluídas. O BacDive se liga ao NCBI pelo TaxID: uma montagem é
selecionada quando o TaxID do organismo, da cepa ou da espécie dela pertence a uma linhagem aceita.

### Arquivo de consulta

As opções de seleção são gravadas em `OUTPUT/genome-selection-query.txt`, que é a entrada do módulo de busca (uma
linha `CHAVE<TAB>VALOR` por filtro). Ele documenta exatamente o que foi pedido e permite repetir a busca com o
módulo isolado:

```
BACDIVE_RANGE	temperature_opt=45:80
BACDIVE_CLASS	cell_shape=rod,coccus
GENOME_TYPE	isolate_genome
```

---

## Estrutura de saída

```
Prokfinder/
├── prokfinder.log
├── INPUT/REFERENCE-METADATA/
│   ├── NCBI/
│   │   ├── BACTERIA-NCBI-GENOMES/fetch.txt            # só o manifesto, não os FASTA
│   │   ├── BACTERIA-NCBI-METADATA/{genbank,refseq}-metadata.tsv
│   │   ├── ARCHAEA-NCBI-GENOMES/  ARCHAEA-NCBI-METADATA/
│   │   └── TAXDUMP/                                   # taxonomia NCBI usada pelo taxonkit
│   ├── BACDIVE/
│   │   ├── BACDIVE-DOWNLOAD-BRUTE/                    # JSON da API (batch_*.jsonl) e lista de IDs
│   │   └── BACDIVE-FORMATTED-TSV/BACDIVE_DATA_FORMATTED.tsv
│   └── GTDB/
│       ├── {ARCHAEA,BACTERIA}-GTDB-TAXONOMY/          # ar53/bac120_taxonomy (bruto + formatado)
│       └── {ARCHAEA,BACTERIA}-GTDB-METADATA/          # ar53/bac120_metadata (bruto + formatado)
└── OUTPUT/
    ├── genome-selection-query.txt
    ├── GENOME-METADATA/genome-metadata.tsv            # tabela única, com as colunas de todas as etapas
    ├── GENOME-FILES/
    │   ├── raw_genomes/<ACESSO>.fna                   # FASTA completo
    │   └── clean_genomes/<ACESSO>.fna                 # sem contigs pequenos, ordenados e renomeados
    ├── GENOME-QUALITY/                                # seqkit, CheckM2, genome-quality.tsv
    ├── GENOME-TAXONOMY/                               # genome-taxonomy.tsv, NON-CLASSIFIED-GENOMES.txt
    ├── GENOME-DIR-HIERARCHY/
    │   └── Domínio/Filo/Classe/Ordem/Família/Gênero/Espécie/
    │       ├── genomes-files-symlink.txt              # acessos e caminhos dos genomas do táxon
    │       ├── GENOME-METADATA/genome-metadata.tsv    # metadados só desse táxon, com o representativo
    │       └── GENOME-DEREPLICATION/                  # saída do skDER
    └── FINAL-RESULTS/
        ├── high-quality-selected-genomes.txt          # acessos + caminho do FASTA dos genomas finais
        ├── high-quality-selected-genomes-metadata.txt
        └── STATISTICS/
            ├── TABLES/   # funil de seleção, QC, genomas por táxon, tabela de todos os procariotos
            └── GRAPHS/   # gráficos (R)
```

Os genomas limpos nunca saem de `clean_genomes/`: a hierarquia taxonômica e as listas finais guardam os caminhos.
O diretório-folha de um genoma é o menor rank conhecido (espécie ou, se ela não for conhecida, gênero, família...).
Genomas sem taxonomia GTDB não entram na hierarquia.

---

## Re-execução

O pipeline pode ser executado de novo sem refazer o que já foi concluído:

| Etapa | O que é reaproveitado |
|---|---|
| Bancos de referência | baixados só uma vez (checkpoint `COMPLETED-DOWNLOAD.txt`); downloads interrompidos são retomados |
| Busca | reaproveitada quando os filtros e os bancos de referência não mudaram |
| Download | genomas já baixados não são baixados de novo |
| Contigs, seqkit, CheckM2, GTDB-Tk | resultados guardados por genoma; só genomas novos ou alterados são processados. Os limites de QC são aplicados a cada execução, então mudá-los não exige rodar o CheckM2 de novo |
| skDER | reaproveitado por táxon quando a lista de genomas e os parâmetros não mudaram |

Para forçar uma etapa:

- `--force-refdb`: baixa de novo os bancos NCBI, BacDive e GTDB
- `--force-selection`: refaz a busca, mesmo com os mesmos filtros
- `--force-curation`: descarta os genomas limpos e todos os resultados de qualidade, taxonomia e desreplicação
  (os genomas baixados são mantidos)

---

## Configuração por variáveis de ambiente

| Variável | Padrão | Etapa |
|---|---|---|
| `QC_MAX_CONTIGS` / `QC_MIN_N50` / `QC_MAX_N_BASES` | 1000 / 5000 / 100000 | 4 |
| `QC_MIN_COMPLETENESS` / `QC_MAX_CONTAMINATION` | 75 / 5 | 4 |
| `CHECKM2_BATCH_SIZE` | 500 | 4 |
| `GTDBTK_EXTRA_ARGS` | vazio (ex.: `--skip_ani_screen`, conforme a versão do GTDB-Tk) | 5 |
| `SKDER_MODE` / `SKDER_PERCENT_IDENTITY` / `SKDER_ALIGNED_FRACTION` | greedy / 99.0 / 90.0 | 7 |
| `SKDER_THREADS_PER_JOB` / `SKDER_EXTRA_ARGS` | 4 / vazio | 7 |
| `NCBI_FIRST_YEAR` | 2017 (o download do NCBI é dividido por ano de publicação) | 0 |

## Executando um módulo isoladamente

Cada módulo é independente e mostra seu uso com `-h`. Exemplo: refazer só a formatação do BacDive a partir do JSON
já baixado.

```bash
modules/download-reference-db-bacdive.sh -o Prokfinder/INPUT/REFERENCE-METADATA/BACDIVE --skip-download
```

Sem o `main-prokfinder.sh`, cada módulo grava seu próprio arquivo `.log` na pasta de saída.

---

## Limitações conhecidas

- O GTDB-Tk precisa de muita RAM (> 200 GB com o pacote completo), por isso fica desativado por padrão.
- Nem todo táxon do GTDB tem correspondência com o NCBI. O filtro `--gtdb-taxon` busca por nome e o vínculo entre
  os dois bancos é feito pelo código de acesso.
- A extração do BacDive foi reescrita em Bash a partir do script #1 do BactoTraits (originalmente em R); um intervalo
  como "5-30 °C" é guardado apenas pelos extremos (5:30).

## Próximos passos

- Construção de um banco de k-mers dos genomas finais (ex.: `kounta`).
- Análises do banco em R.
