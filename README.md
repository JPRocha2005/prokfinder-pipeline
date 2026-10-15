# Prokfinder

Pipeline para **busca, download e curadoria de genomas de procariotos** (Bacteria e Archaea) do NCBI, combinando
os metadados das montagens (NCBI), os dados fenotípicos das linhagens (BacDive) e a taxonomia filogenômica
(GTDB).

O objetivo é montar um banco de genomas **de alta qualidade, classificados taxonomicamente e sem redundância**,
selecionados por critérios biológicos. Por exemplo: todos os procariotos com temperatura ótima de crescimento
entre 45 e 80 °C isolados de fontes termais.

---

## Sub-comandos

O Prokfinder é organizado em **sub-comandos**. Cada um tem as suas opções e a sua ajuda (`<sub-comando> -h`):

| Sub-comando | O que faz |
|---|---|
| `setup` | Instala os bancos de referência (NCBI, BacDive e GTDB). **Obrigatório antes de `search-genomes` e `download-wf`** (uma única vez); também atualiza bancos instalados por uma versão anterior do Prokfinder |
| `search-genomes` | Seleciona os genomas pelos filtros (NCBI, BacDive e GTDB), combinados com `AND`/`OR` |
| `download-wf` | Baixa os genomas selecionados, remove contigs pequenos, filtra por qualidade, atribui a taxonomia GTDB e desreplica cada táxon. Precisa de um `search-genomes` concluído |
| `run-on-cluster` | Executa os outros sub-comandos como um *job* em um nó de cálculo do cluster (PBS) |

Vários sub-comandos podem ser dados juntos; eles sempre rodam na ordem `setup` → `search-genomes` → `download-wf`:

```bash
./main-prokfinder.sh setup
./main-prokfinder.sh search-genomes -T 45:80 AND -O rod,coccus
./main-prokfinder.sh download-wf --min-completeness 90

# Tudo de uma vez
./main-prokfinder.sh setup search-genomes -T 45:80 AND -O rod,coccus download-wf --min-completeness 90
```

**Opções globais** (aceitas em qualquer posição, com qualquer sub-comando):

| Opção | Padrão | Descrição |
|---|---|---|
| `-o, --results-dir <dir>` | `./prokfinder-results` | pasta de resultados desta execução |
| `--refdb-dir <dir>` | `~/prokfinder-refdb` | pasta dos bancos de referência, compartilhada por todas as execuções |
| `-p, --threads <n>` | 4 (com `run-on-cluster`: as CPUs do job) | CPUs usadas por todos os programas (seqkit, CheckM2, GTDB-Tk, skDER) |
| `-h, --help` | | ajuda geral ou, depois de um sub-comando, a ajuda dele |

Erros na linha de comando mostram só a mensagem (em vermelho) e a dica de qual `-h` consultar; a ajuda completa
aparece apenas com `-h`/`--help`. Uma opção dada no sub-comando errado também é indicada (ex.: `--min-completeness`
pertence a `download-wf`).

---

## Etapas e módulos

| Etapa | Sub-comando | Módulo | O que faz |
|---|---|---|---|
| 0 | `setup` | `download-reference-db-ncbi.sh` | Baixa os metadados de todas as montagens de Bacteria/Archaea (NCBI Datasets, pacote *dehydrated*: só o `fetch.txt`, sem os FASTA). Une as linhas repetidas da mesma montagem, formata o cabeçalho, troca valores vazios ou inválidos por `NULL`, adiciona a linhagem NCBI (nome e TaxID de cada rank e o TaxID da cepa ou `NO_TAXID`) com o taxonkit e separa GenBank (`GCA_`) de RefSeq (`GCF_`). Grava também uma tabela **completa** por domínio (`complete-metadata.tsv`, todos os campos do NCBI) |
| 0 | `setup` | `download-reference-db-bacdive.sh` | Baixa todas as linhagens do BacDive pela API e extrai 29 traços fenotípicos, com o rank de cada TaxID. Une as referências de cada cepa e depois as cepas de cada TaxID, no estilo do BactoTraits (ver abaixo) |
| 0 | `setup` | `download-reference-db-gtdb.sh` | Baixa a taxonomia e os metadados do GTDB (`ar53`/`bac120`, última release), confere o MD5 e separa a linhagem em uma coluna por rank; a coluna `accession` do GTDB vira `gtdb_accession` |
| 1 | `search-genomes` | `search-genomes-in-reference-db.sh` | Busca os genomas que satisfazem a expressão de filtros |
| 2 | `download-wf` | `download-selected-genomes-from-db.sh` | Baixa só os genomas selecionados (`datasets rehydrate`), apenas a sequência genômica (`*_genomic.fna`) |
| 3 | `download-wf` | `filter-low-quality-contigs-from-genomes.sh` | Remove contigs pequenos (padrão < 1000 pb), ordena por tamanho e renomeia (`Contig_1`, ...) com o seqkit |
| 4 | `download-wf` | `evaluate-and-filter-genome-by-qc.sh` | seqkit (tamanho, nº de contigs, N50, bases N, GC) + CheckM2 (completude, contaminação) e filtragem por qualidade |
| 5 | `download-wf` | `assign-taxonomy-gtdb.sh` | Taxonomia GTDB pela tabela do banco; os não encontrados vão para `NON-CLASSIFIED-GENOMES-BY-GTDB.txt` ou, com `--gtdbtk`, são classificados pelo GTDB-Tk |
| 6 | `download-wf` | `build-taxonomic-dir-hierarchy.sh` | Hierarquia de diretórios `Domínio/Filo/Classe/Ordem/Família/Gênero/Espécie` |
| 7 | `download-wf` | `dereplicate-genomes-skder.sh` | Desreplicação com o skDER dentro de cada táxon (diretório-folha) |
| 8 | `download-wf` | `summarize-final-results.sh` | Lista final de genomas, tabelas de estatística e gráficos |
| | `run-on-cluster` | `submit-prokfinder-pbs.sh` | Escreve o script `.pbs` e o envia com `qsub` |

Cada etapa acrescenta suas colunas a uma **única tabela de metadados** (`OUTPUT/GENOME-METADATA/genome-metadata.tsv`,
uma linha por genoma): metadados do NCBI, colunas de seleção, traços do BacDive, status do download, métricas de
qualidade, taxonomia GTDB, diretório taxonômico e a indicação de genoma representativo.

---

## Dependências

| Programa | Usado em | Instalação (conda) |
|---|---|---|
| bash ≥ 4, awk, gawk, curl, unzip, gzip | todas as etapas | sistema |
| `datasets`, `dataformat` (NCBI Datasets) | `setup` e `download-wf` | `conda install -c conda-forge ncbi-datasets-cli` |
| `taxonkit` | `setup` (NCBI e BacDive) e filtro `-N` | `conda install -c bioconda taxonkit` |
| `jq` | `setup` (BacDive) | `conda install -c conda-forge jq` |
| `seqkit` | `download-wf` | `conda install -c bioconda seqkit` |
| `checkm2` + banco de dados | `download-wf` | `conda install -c bioconda checkm2` e `checkm2 database --download` |
| `skder` | `download-wf` | `conda install -c bioconda skder` |
| `gtdbtk` + pacote de referência (opcional) | `download-wf --gtdbtk` | `conda install -c bioconda gtdbtk` |
| `Rscript` (opcional, só R base) | gráficos do `download-wf` | `conda install -c conda-forge r-base` |

Os parâmetros do Prokfinder são dados **só por sub-comandos e flags**. As únicas variáveis de ambiente lidas são as
de configuração dos programas, que no cluster os `module load` já definem:

| Variável | Quando | Descrição |
|---|---|---|
| `CHECKM2DB` | `download-wf` | caminho do **arquivo** do banco do CheckM2 (`uniref100.KO.1.dmnd`) |
| `GTDBTK_DATA_PATH` | `download-wf --gtdbtk` | pasta do pacote de referência do GTDB-Tk |
| `NCBI_API_KEY` | opcional | chave da API do NCBI (aumenta o limite de 5 para 10 requisições/s) |

O `main-prokfinder.sh` confere todas as dependências necessárias **antes** de começar e lista todos os problemas de
uma vez.

---

## setup: bancos de referência

```bash
./main-prokfinder.sh setup                 # instala (ou atualiza) os bancos em ~/prokfinder-refdb
./main-prokfinder.sh setup --force-refdb   # baixa e instala tudo de novo
```

- Os bancos são baixados **uma única vez** e compartilhados por todas as execuções e pastas de resultados. O download
  leva horas; no cluster, use `run-on-cluster setup`.
- `search-genomes` e `download-wf` **não rodam sem o `setup`**: se os bancos não forem encontrados, o Prokfinder
  para com o erro `Reference databases not found in ~/prokfinder-refdb` e pede o `setup` (que pode vir no mesmo
  comando: `setup search-genomes ...`).
- Rodar `setup` de novo com os bancos já instalados não baixa nada. Para reinstalar é preciso `--force-refdb`.
- Bancos instalados por uma versão **anterior** do Prokfinder, cujas tabelas tinham outro formato, também precisam do
  `setup`: ele **atualiza** as tabelas a partir dos arquivos já baixados, **sem baixar nada de novo** (a versão do
  formato fica registrada em `COMPLETED-DOWNLOAD.txt`, linha `format=`).
- Se duas execuções precisarem instalar os bancos ao mesmo tempo, a segunda espera a primeira terminar.

**Tabelas do BacDive (estilo BactoTraits).** Uma cepa pode ter vários valores para o mesmo traço, vindos de
referências diferentes (em `BACDIVE_DATA_FORMATTED.tsv`, `;` separa referências e `:` separa valores da mesma
referência). Cada referência e depois cada cepa tem o mesmo peso:

| Tabela | Uma linha por | Números | `*_range` | Classes |
|---|---|---|---|---|
| `BACDIVE_CONCATENATED_REFERENCES_FOR_STRAIN.tsv` | cepa (`bacdive_id`) | média das referências | menor:maior valor | fração das referências que registraram cada classe (ex.: `coccus=0.25;rod=0.75`) |
| `BACDIVE_CONCATENATED_STRAINS_FOR_SPECIES.tsv` | TaxID do NCBI (`ncbi_tax_id`, `ncbi_tax_id_rank`, `bacdive_strains`) | média das cepas | menor:maior valor | média das frações das cepas |

---

## search-genomes: seleção dos genomas

```bash
./main-prokfinder.sh search-genomes -T 45:80 AND -O rod,coccus AND -g isolate_genome
./main-prokfinder.sh search-genomes -N taxons.txt OR "(" --temp-growth 70: AND -i "hot spring","hydrothermal vent" ")"
./main-prokfinder.sh -o termofilos search-genomes --gtdb-taxon g__Thermus AND -c Brazil,Iceland
```

### Como escrever os filtros

1. Cada flag de seleção é **um filtro**; é preciso pelo menos um.
2. Vários valores da **mesma flag** são separados por **vírgula, sem espaço**, e combinados com **OR** (o filtro
   passa se qualquer um dos valores corresponder): `-t 562,1280`. Cada valor com espaço fica entre as suas próprias
   aspas duplas: `-i "hot spring","hydrothermal vent",geyser`.
3. **Entre dois filtros é obrigatório escrever `AND` ou `OR`**; sem isso o Prokfinder dá erro:

   ```bash
   -i "hot spring",geyser AND -c Brazil    # fonte termal ou gêiser, E do Brasil
   -i "hot spring",geyser OR -c Brazil     # fonte termal ou gêiser, OU do Brasil
   ```

4. `AND` é avaliado antes de `OR` (`A AND B OR C` = `(A AND B) OR C`). **Parênteses** agrupam filtros; no terminal
   eles precisam estar entre aspas ou com `\`: `-T 45:80 AND "(" -O rod OR -M anaerobe ")"`.
5. Arquivos de lista (`-f`, `-A`, `-N`, `-I`, `-C`): um valor por linha. Os `\r` (fim de linha do Windows) e os
   espaços no começo/fim das linhas são removidos, e linhas vazias são ignoradas.

### Filtros

| Opção | Banco | Descrição |
|---|---|---|
| `-t, --taxid` / `-f, --taxid-file` | NCBI | TaxID do organismo **ou de qualquer rank da linhagem** (ex.: `1239` seleciona todos os Bacillota) |
| `-a, --accession` / `-A, --accession-file` | NCBI | código de acesso **com versão** (`GCA_000005845.2`). GCA_ e GCF_ podem vir misturados (GCF_x.1 seleciona o par GenBank GCA_x.1). Todas as linhas inválidas são listadas no erro |
| `--ncbi-taxon` | NCBI | nome de um táxon da linhagem NCBI, qualquer rank |
| `-N, --ncbi-taxon-file` | NCBI | lista de nomes de táxons (um por linha), convertidos em TaxID com `taxonkit name2taxid`; um nome que não é táxon de Bacteria/Archaea dá erro (com a lista dos nomes) |
| `-i` / `-I` | NCBI | texto contido na fonte de isolamento (ex.: `"hot spring","hydrothermal vent"`) |
| `-c` / `-C` | NCBI | país/local da localização geográfica (ex.: `Brazil,Iceland`) |
| `-g` | NCBI | `isolate_genome` ou `mag` |
| `-G` | NCBI | faixa de GC (%) da montagem (`min:max`) |
| `-T` / `-P` / `-S` | BacDive | temperatura, pH e NaCl **ótimos** (`min:max`, um dos lados pode ficar vazio) |
| `--temp-growth` / `--ph-growth` / `--salt-growth` | BacDive | valores em que houve crescimento |
| `-O` / `-M` / `-B traço=classe` | BacDive | formato celular, tolerância a O₂ e qualquer traço dividido em classes (Gram, motilidade, esporulação, fonte de energia/carbono...) |
| `--gtdb-taxon` | GTDB | nome de um táxon da linhagem GTDB, qualquer rank (`g__Thermus` ou `Thermus`) |

Outras opções: `-s genbank|refseq` (tabela do NCBI pesquisada; padrão `genbank`) e `--force-selection`.
`./main-prokfinder.sh search-genomes -h` mostra as classes disponíveis para cada traço.

**Filtros do BacDive:** a busca usa a tabela por TaxID. Um filtro numérico passa quando a **média** está no
intervalo; um filtro de classe, quando a classe foi registrada para o TaxID. A ligação com os genomas depende do
**rank** do TaxID: de **espécie** seleciona os genomas da espécie; **abaixo da espécie** (cepa, subespécie...)
seleciona os genomas daquele organismo/cepa; **acima da espécie** (gênero...) não é usado. O log mostra, por filtro,
quantos TaxIDs passaram por rank (ex.: `BacDive temperature_opt 45:80 -> 231 NCBI TaxIDs match (200 species, 31 strain)`).
Atenção: no BacDive a maioria das cepas registra só temperaturas de **cultivo**; poucas registram a temperatura
**ótima**, por isso `-T` encontra bem menos espécies que `--temp-growth`.

### Resultado

`OUTPUT/GENOME-METADATA/genome-metadata.tsv` tem uma linha por genoma selecionado e, para cada filtro, uma coluna
`TRUE`/`FALSE` (com `OR`, um genoma selecionado pode ter `FALSE` em algum filtro):

| Tipo de flag | Coluna | Exemplo |
|---|---|---|
| opção única | `SELECTED_BY_<FLAG>-<OPÇÃO>` | `SELECTED_BY_GENOME_TYPE-ISOLATE-GENOME` |
| numérica | `SELECTED_BY_<FLAG>-<FAIXA>` | `SELECTED_BY_OPT_TEMP-45:80`, `SELECTED_BY_GC-40:65` |
| vários valores | `SELECTED_BY_<FLAG>` + coluna com os valores que deram match (separados por `;`) | `SELECTED_BY_ISOLATION_SOURCE` + `ISOLATION_SOURCE_KEYWORD_MATCH` |

Se a mesma flag aparecer duas vezes, a segunda coluna ganha `_2`. As colunas do BacDive são `BACDIVE_TAXID`,
`BACDIVE_LINK_LEVEL` (rank), `BACDIVE_STRAINS` e os traços.

Ao terminar, a busca grava o checkpoint **`SEARCH-GENOMES-COMPLETED.txt`** na pasta de resultados (data, expressão,
tabela do NCBI usada e número de genomas). O `download-wf` só roda se esse arquivo existir. A busca é reaproveitada
quando os filtros e os bancos não mudaram; `--force-selection` remove o checkpoint e busca de novo.

A consulta fica gravada em `OUTPUT/genome-selection-query.txt` (um filtro por linha e a expressão), que é a entrada
do módulo de busca:

```
FILTER	F1	BACDIVE_RANGE	temperature_opt=45:80
FILTER	F2	BACDIVE_CLASS	cell_shape=rod,coccus
FILTER	F3	ISOLATION_SOURCE	hot spring
FILTER	F3	ISOLATION_SOURCE	geyser
EXPRESSION	F1 AND ( F2 OR F3 )
```

---

## download-wf: download, qualidade, taxonomia e desreplicação

```bash
./main-prokfinder.sh download-wf
./main-prokfinder.sh download-wf --min-completeness 90 --max-contamination 2 --derep-mode dynamic --filter-mge
```

Usa os genomas do último `search-genomes` da mesma pasta de resultados (`-o`). Opções:

| Etapa | Opção | Padrão | Descrição |
|---|---|---|---|
| Download | `-j, --download-threads <n>` | 10 | downloads simultâneos (1-30) |
| | `-z, --gzip` | | guarda os genomas comprimidos (`.fna.gz`) |
| Contigs | `--min-contig-length <pb>` | 1000 | contigs menores são removidos |
| Qualidade | `--max-contigs <n>` | 1000 | o genoma precisa ter menos contigs que isso |
| | `--max-n-bases <n>` | 100000 | o genoma precisa ter menos bases ambíguas (N) que isso |
| | `--min-n50 <pb>` | 5000 | N50 mínimo |
| | `--min-completeness <%>` | 75 | completude mínima (CheckM2) |
| | `--max-contamination <%>` | 5 | contaminação máxima (CheckM2) |
| Taxonomia | `--gtdbtk` | desligado | classifica com o GTDB-Tk os genomas que não estão nas tabelas do GTDB (**precisa de > 200 GB de RAM**); sem ela, eles só vão para `NON-CLASSIFIED-GENOMES-BY-GTDB.txt` |
| Desreplicação | `--derep-ani <%>` | 99,5 | identidade (ANI) acima da qual dois genomas do mesmo táxon são redundantes (padrão do skDER); valores menores deixam menos representantes |
| | `--derep-mode <modo>` | `greedy` | modo do skDER: `greedy` ou `dynamic` (fração alinhada fixa em 90 %, padrão do skDER) |
| | `--filter-mge` | desligado | o skDER mascara os elementos genéticos móveis (profagos e plasmídeos, previstos com o PhiSpy) antes de comparar os genomas |
| Re-execução | `--force-curation` | | descarta os genomas limpos e todos os resultados de qualidade, taxonomia e desreplicação (os genomas baixados são mantidos) |

Os limites de qualidade são reaplicados a cada execução sobre os resultados já guardados: mudá-los não faz o CheckM2
rodar de novo.

---

## run-on-cluster: execução no cluster (PBS)

Em um cluster com o gerenciador de filas **PBS** (como o `cluster.dti.ufv.br` da UFV), o Prokfinder não deve rodar
no nó de login: no cluster da UFV, qualquer comando que passe de **10 minutos** no nó de login é interrompido. Por
isso, no cluster, coloque `run-on-cluster` antes dos outros sub-comandos:

```bash
cd /pasta/de/trabalho                      # o job roda na pasta de onde foi submetido
./prokfinder-scripts/main-prokfinder.sh run-on-cluster setup
./prokfinder-scripts/main-prokfinder.sh run-on-cluster --queue qtime --walltime 24:00:00 search-genomes -T 45:80 AND -O rod
./prokfinder-scripts/main-prokfinder.sh run-on-cluster --queue high_mem --mem 256gb --cpus 32 download-wf --gtdbtk
```

O que acontece:

1. A linha de comando é conferida na hora (um filtro inválido dá erro antes de ir para a fila).
2. O módulo [`modules/submit-prokfinder-pbs.sh`](modules/submit-prokfinder-pbs.sh) escreve
   `<resultados>/PBS-JOBS/prokfinder_<data>_<hora>.pbs`, que carrega os módulos do cluster
   (`module load seqkit checkm2 ncbi_datasets gtdbtk skder taxonkit`, além do `r` para os gráficos) e executa o
   `main-prokfinder.sh` com a mesma linha de comando, sem o `run-on-cluster`.
3. O script é enviado à fila com `qsub` (se o `qsub` não existir, o Prokfinder sai com erro).
4. A saída do job fica em `<resultados>/PBS-JOBS/prokfinder_<data>_<hora>.log` e o log completo em
   `<resultados>/prokfinder.log`. Acompanhe com `qstat -u $USER`.

| Opção | Padrão | Descrição |
|---|---|---|
| `--queue <fila>` | `roteamento` | fila de submissão (ver abaixo) |
| `--mem <ram>` | `20gb` | memória RAM (`20gb`, `24000mb` ou só o número em GB) |
| `--cpus <n>` | 20 | CPUs do job; o Prokfinder usa todas, a não ser que `-p` seja dado |
| `--walltime <HH:MM:SS>` | `2160:00:00` (`72:00:00` na `qtime`) | tempo máximo do job |
| `--test` | | envia só um job curto (fila `qtime`, 15 min) que carrega os módulos e mostra a versão de cada programa |
| `--dry-run` | | só escreve o `.pbs`, sem submeter |

| Fila | Uso | Limites |
|---|---|---|
| `qtime` | testes e execuções curtas | walltime obrigatório, no máximo 72:00:00; até 24024 MB por job |
| `roteamento` | execuções comuns (padrão) | até 24024 MB por job; walltime até 2160:00:00 |
| `high_mem` | execuções com muita memória (ex.: `--gtdbtk`) | no mínimo 25 GB; até 512 GB e 64 CPUs por usuário; 2 jobs rodando |

Os limites de cada fila são conferidos antes da submissão. A `NCBI_API_KEY` é repassada pelo `qsub` e nunca é
gravada no arquivo `.pbs`.

---

## Estrutura de saída

Bancos de referência (compartilhados, `--refdb-dir`):

```
~/prokfinder-refdb/
├── NCBI/
│   ├── BACTERIA-NCBI-GENOMES/fetch.txt                # só o manifesto, não os FASTA
│   ├── BACTERIA-NCBI-METADATA/{genbank,refseq}-metadata.tsv   # tabelas usadas na busca (campos selecionados)
│   ├── BACTERIA-NCBI-METADATA/complete-metadata.tsv   # todos os campos do NCBI (GCA_ e GCF_)
│   ├── ARCHAEA-NCBI-GENOMES/  ARCHAEA-NCBI-METADATA/
│   └── TAXDUMP/                                       # taxonomia NCBI usada pelo taxonkit
├── BACDIVE/
│   ├── BACDIVE-DOWNLOAD-BRUTE/                        # JSON da API (batch_*.jsonl) e lista de IDs
│   └── BACDIVE-FORMATTED-TSV/
│       ├── BACDIVE_DATA_FORMATTED.tsv                 # uma linha por cepa, todas as referências
│       ├── BACDIVE_CONCATENATED_REFERENCES_FOR_STRAIN.tsv   # referências de cada cepa unidas
│       └── BACDIVE_CONCATENATED_STRAINS_FOR_SPECIES.tsv     # cepas de cada TaxID unidas (usada na busca)
└── GTDB/
    ├── {ARCHAEA,BACTERIA}-GTDB-TAXONOMY/              # ar53/bac120_taxonomy (bruto + formatado)
    └── {ARCHAEA,BACTERIA}-GTDB-METADATA/              # ar53/bac120_metadata (bruto + formatado)
```

Resultados de uma execução (`-o`/`--results-dir`):

```
prokfinder-results/
├── prokfinder.log
├── SEARCH-GENOMES-COMPLETED.txt                       # checkpoint da busca (necessário para o download-wf)
├── PBS-JOBS/                                          # scripts .pbs e logs dos jobs (só com run-on-cluster)
└── OUTPUT/
    ├── genome-selection-query.txt
    ├── GENOME-METADATA/genome-metadata.tsv            # tabela única, com as colunas de todas as etapas
    ├── GENOME-FILES/
    │   ├── raw_genomes/<ACESSO>.fna                   # FASTA completo
    │   └── clean_genomes/<ACESSO>.fna                 # sem contigs pequenos, ordenados e renomeados
    ├── GENOME-QUALITY/                                # seqkit, CheckM2, genome-quality.tsv
    ├── GENOME-TAXONOMY/                               # genome-taxonomy.tsv, NON-CLASSIFIED-GENOMES-BY-GTDB.txt
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

O Prokfinder pode ser executado de novo sem refazer o que já foi concluído:

| Etapa | O que é reaproveitado | Para forçar |
|---|---|---|
| `setup` | bancos baixados só uma vez (checkpoint `COMPLETED-DOWNLOAD.txt`); downloads interrompidos são retomados | `setup --force-refdb` |
| `search-genomes` | busca reaproveitada quando os filtros e os bancos não mudaram (checkpoint `SEARCH-GENOMES-COMPLETED.txt`) | `search-genomes --force-selection` |
| `download-wf` | genomas já baixados; contigs, seqkit, CheckM2 e GTDB-Tk por genoma; skDER por táxon (quando a lista de genomas e os parâmetros não mudaram) | `download-wf --force-curation` |

---

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
