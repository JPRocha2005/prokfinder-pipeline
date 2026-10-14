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
│   ├── plot-statistics.R
│   └── submit-prokfinder-pbs.sh    # execução no cluster (PBS): flag --pbs
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

O Prokfinder usa duas pastas:

- **Bancos de referência** (NCBI, BacDive e GTDB): baixados **uma única vez** para `~/prokfinder-refdb` e
  compartilhados por todas as execuções e pastas de resultados. O download leva horas e ocupa alguns GB. Outro local
  pode ser escolhido com `--refdb-dir` ou com a variável `PROKFINDER_REFDB_DIR`. Se duas execuções precisarem baixar
  os bancos ao mesmo tempo, a segunda espera a primeira terminar e reaproveita o download.
- **Resultados de cada execução** (genomas selecionados, qualidade, taxonomia, resultados finais): `./prokfinder-results`,
  ou o nome escolhido com `-o`/`--results-dir`. Assim, cada consulta pode ter a sua pasta sem baixar os bancos de
  novo. O log completo fica em `prokfinder-results/prokfinder.log`.

```bash
./main-prokfinder.sh -T 45:80 -o termofilos            # resultados em ./termofilos
./main-prokfinder.sh -T 80: -o hipertermofilos         # outra consulta, mesmos bancos de referência
```

### Filtros de seleção

Regras:

1. É necessário **pelo menos uma flag de seleção** (NCBI, BacDive ou GTDB).
2. Para passar vários valores para a mesma flag, separe-os por **vírgula, sem espaço** (`-t 562,1280`). Um valor
   que contém espaço deve ficar **entre aspas duplas** (`-i "hot spring,geyser"`).
3. Um genoma é selecionado quando corresponde a **pelo menos um** dos valores de uma flag (ou das linhas do seu
   arquivo de lista) e, ao mesmo tempo, a **todas** as flags usadas. Exemplo: `-i "hot spring,geyser" -c Brazil`
   seleciona genomas isolados de fonte termal **ou** de gêiser, mas só os do Brasil.

Erros na linha de comando (flag desconhecida, flag sem valor, valor com espaço sem aspas) mostram só a mensagem de
erro, em vermelho; a ajuda completa aparece apenas com `-h` ou `--help`.

| Opção | Banco | Descrição |
|---|---|---|
| `-t, --taxid` / `-f, --taxid-file` | NCBI | TaxID do organismo **ou de qualquer rank da linhagem** (ex.: `1239` seleciona todos os Bacillota) |
| `-a, --accession` / `-A, --accession-file` | NCBI | código de acesso GCA_/GCF_ (versão opcional; GCA_x.1 também encontra o par GCF_x.1) |
| `--ncbi-taxon` | NCBI | nome de um táxon da linhagem NCBI, qualquer rank |
| `-i` / `-I` | NCBI | texto contido na fonte de isolamento (ex.: `"hot spring,hydrothermal vent"`) |
| `-c` / `-C` | NCBI | país/local da localização geográfica (ex.: `Brazil,Iceland`) |
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

Bancos de referência (compartilhados, `--refdb-dir`):

```
~/prokfinder-refdb/
├── NCBI/
│   ├── BACTERIA-NCBI-GENOMES/fetch.txt                # só o manifesto, não os FASTA
│   ├── BACTERIA-NCBI-METADATA/{genbank,refseq}-metadata.tsv
│   ├── ARCHAEA-NCBI-GENOMES/  ARCHAEA-NCBI-METADATA/
│   └── TAXDUMP/                                       # taxonomia NCBI usada pelo taxonkit
├── BACDIVE/
│   ├── BACDIVE-DOWNLOAD-BRUTE/                        # JSON da API (batch_*.jsonl) e lista de IDs
│   └── BACDIVE-FORMATTED-TSV/BACDIVE_DATA_FORMATTED.tsv
└── GTDB/
    ├── {ARCHAEA,BACTERIA}-GTDB-TAXONOMY/              # ar53/bac120_taxonomy (bruto + formatado)
    └── {ARCHAEA,BACTERIA}-GTDB-METADATA/              # ar53/bac120_metadata (bruto + formatado)
```

Resultados de uma execução (`-o`/`--results-dir`):

```
prokfinder-results/
├── prokfinder.log
├── PBS-JOBS/                                          # scripts .pbs e logs dos jobs (só com --pbs)
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
| Bancos de referência | baixados só uma vez em `~/prokfinder-refdb` e compartilhados por todas as pastas de resultados (checkpoint `COMPLETED-DOWNLOAD.txt`); downloads interrompidos são retomados |
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
| `GTDBTK_EXTRA_ARGS` | vazio (ex.: `--scratch_dir /tmp/gtdbtk_scratch`) | 5 |
| `SKDER_MODE` / `SKDER_PERCENT_IDENTITY` / `SKDER_ALIGNED_FRACTION` | greedy / 99.5 (padrão do skDER) / 90.0 | 7 |
| `SKDER_THREADS_PER_JOB` / `SKDER_EXTRA_ARGS` | 4 / vazio | 7 |
| `NCBI_FIRST_YEAR` | 2017 (o download do NCBI é dividido por ano de publicação) | 0 |
| `PROKFINDER_REFDB_DIR` | `~/prokfinder-refdb` (mesmo que `--refdb-dir`) | 0 |

Essas variáveis são lidas por cada módulo; quando não estão definidas, vale o padrão da tabela. Defina-as na mesma
linha do comando (vale só para aquela execução) ou com `export` (vale para a sessão do terminal):

```bash
# Exigir completude ≥ 90 % e contaminação ≤ 2 % só nesta execução
QC_MIN_COMPLETENESS=90 QC_MAX_CONTAMINATION=2 ./main-prokfinder.sh -T 45:80
```

Os limites de QC são reaplicados a cada execução sobre os resultados já guardados, então mudá-los não faz o CheckM2
rodar de novo.

### Opções extras dos programas (`SKDER_EXTRA_ARGS` e `GTDBTK_EXTRA_ARGS`)

O Prokfinder passa ao skDER e ao GTDB-Tk só as opções principais. Qualquer outra opção desses programas pode ser
acrescentada pelas variáveis `*_EXTRA_ARGS`: o texto é colado, sem alteração, no fim do comando do programa. Assim
dá para usar opções novas sem mudar o código.

**Exemplo:** ligar o filtro de elementos genéticos móveis (MGE) do skDER, que vem **desligado** por padrão. Com ele,
o skDER prevê e mascara profagos e plasmídeos antes de comparar os genomas.

```bash
# MGEs previstos com o PhiSpy (padrão do skDER; o PhiSpy precisa estar instalado)
SKDER_EXTRA_ARGS="--filter-mge" ./main-prokfinder.sh -T 45:80 -p 32

# MGEs previstos com o geNomad (no cluster da UFV, o banco indicado pelo módulo skder/1.3.6)
SKDER_EXTRA_ARGS="--filter-mge --genomad-database /data/db/skder/1.3.6/genomad/1.5" ./main-prokfinder.sh --pbs -T 45:80
```

As opções `--filter-mge` e `--genomad-database` foram conferidas no skDER 1.3.6 do cluster.

O comando executado fica, por exemplo:

```
skder --genomes ... --dereplication-mode greedy ... --threads 4 --symlink --filter-mge
```

Para o GTDB-Tk (no cluster, versão 2.7.2, que não precisa de opção extra), por exemplo para reduzir o uso de
memória do pplacer gravando em disco:

```bash
GTDBTK_EXTRA_ARGS="--scratch_dir /tmp/gtdbtk_scratch" ./main-prokfinder.sh -T 45:80 --gtdbtk
```

Cuidados:

- Os valores são separados por espaço, então um caminho com espaço não funciona dentro dessas variáveis.
- Mudar `SKDER_EXTRA_ARGS` faz o skDER rodar de novo nos táxons afetados (os parâmetros entram no cache). Mudar
  `GTDBTK_EXTRA_ARGS` **não** reclassifica genomas já classificados; para isso use `--force-curation`.

## Execução no cluster (PBS)

Em um cluster com o gerenciador de filas **PBS** (como o `cluster.dti.ufv.br` da UFV), o Prokfinder não deve rodar
no nó de login: ele precisa ser enviado como um *job* para um nó de cálculo. No cluster da UFV, qualquer comando
que passe de **10 minutos** no nó de login é interrompido, e o download dos bancos de referência leva horas. Por
isso, no cluster, rode sempre com `--pbs`, inclusive na primeira execução. A flag `--pbs` faz isso:

```bash
cd /pasta/de/trabalho                      # o job roda na pasta de onde foi submetido
./prokfinder-scripts/main-prokfinder.sh --pbs -T 45:80 -O rod,coccus
```

O que acontece:

1. As opções são conferidas na hora (um filtro inválido dá erro antes de ir para a fila).
2. O módulo [`modules/submit-prokfinder-pbs.sh`](modules/submit-prokfinder-pbs.sh) escreve um script de
   submissão em `prokfinder-results/PBS-JOBS/prokfinder_<data>_<hora>.pbs`. O script carrega os módulos do cluster
   (`module load seqkit checkm2 ncbi_datasets gtdbtk skder taxonkit`, além do `r` para os gráficos) e executa o
   `main-prokfinder.sh` com as mesmas opções, sem o `--pbs`.
3. O script é enviado à fila com `qsub`. Se o comando `qsub` não existir, o Prokfinder sai com erro.
4. A saída do job fica em `prokfinder-results/PBS-JOBS/prokfinder_<data>_<hora>.log`, e o log completo do pipeline em
   `prokfinder-results/prokfinder.log`. Acompanhe com `qstat -u $USER`.

Recursos padrão: **fila `roteamento`, 20 CPUs (`ppn`), 20 GB de RAM, walltime 2160:00:00**. O Prokfinder usa todas
as CPUs do job (`-p`), a não ser que você informe outro `-p`. Para mudar os recursos:

| Variável | Padrão | Descrição |
|---|---|---|
| `PROKFINDER_PBS_QUEUE` | `roteamento` | fila: `roteamento`, `qtime` ou `high_mem` |
| `PROKFINDER_PBS_PPN` | `20` | número de CPUs |
| `PROKFINDER_PBS_MEM` | `20gb` | memória RAM (`20gb`, `24000mb`...) |
| `PROKFINDER_PBS_WALLTIME` | `2160:00:00` (`72:00:00` na `qtime`) | tempo máximo do job |
| `PROKFINDER_PBS_JOB_NAME` | `prokfinder` | nome do job |
| `PROKFINDER_PBS_MODULES` | `seqkit checkm2 ncbi_datasets gtdbtk skder taxonkit` | módulos carregados no job |

```bash
# Mais CPUs e o máximo de RAM da fila roteamento
PROKFINDER_PBS_PPN=32 PROKFINDER_PBS_MEM=24000mb ./prokfinder-scripts/main-prokfinder.sh --pbs -T 45:80

# GTDB-Tk precisa de > 200 GB de RAM: fila high_mem
PROKFINDER_PBS_QUEUE=high_mem PROKFINDER_PBS_MEM=256gb PROKFINDER_PBS_PPN=32 \
    ./prokfinder-scripts/main-prokfinder.sh --pbs -T 45:80 --gtdbtk
```

Os limites de cada fila são conferidos antes da submissão: até 24024 MB por job nas filas `roteamento` e `qtime`,
no mínimo 25 GB na `high_mem` e walltime de no máximo 72 h na `qtime`. Os ajustes do Prokfinder definidos no
terminal (`QC_*`, `SKDER_*`, `GTDBTK_EXTRA_ARGS`...) são copiados para o job; a `NCBI_API_KEY` é repassada pelo
`qsub` e nunca é gravada no arquivo `.pbs`.

Outras opções do módulo (`modules/submit-prokfinder-pbs.sh -h` mostra todas):

```bash
# Job curto (fila qtime, 15 min) que só carrega os módulos e mostra a versão/ajuda de cada programa
./prokfinder-scripts/modules/submit-prokfinder-pbs.sh --test

# Só escrever o .pbs, sem submeter (para conferir ou submeter depois com qsub)
./prokfinder-scripts/modules/submit-prokfinder-pbs.sh --dry-run -- -T 45:80
```

## Executando um módulo isoladamente

Cada módulo é independente e mostra seu uso com `-h`. Exemplo: refazer só a formatação do BacDive a partir do JSON
já baixado.

```bash
modules/download-reference-db-bacdive.sh -o ~/prokfinder-refdb/BACDIVE --skip-download
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
