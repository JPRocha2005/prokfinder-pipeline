# Prokfinder

Pipeline para **busca, download e curadoria de genomas de procariotos** (Bacteria e Archaea) do NCBI, combinando
os metadados das montagens (NCBI), os dados fenotípicos das linhagens (BacDive) e a taxonomia filogenômica
(GTDB).

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
| `--memory <ram>` | RAM do computador (com `run-on-cluster`: a RAM do job) | RAM que os programas podem usar (`16gb`, `16000mb` ou só o número em GB). Com menos de 16 GB o CheckM2 roda em modo de pouca memória (`--lowmem`); o skDER/CiDDER recebem o limite (`--max-memory`); com menos de 200 GB o posicionamento na árvore do GTDB-Tk grava em disco (`--scratch_dir`) |
| `--check` | | mostra quais dependências estão instaladas e sai (ver [Dependências](#dependências)) |
| `-h, --help` | | ajuda geral ou, depois de um sub-comando, a ajuda dele |

---

## Instalação

O Prokfinder é um conjunto de scripts em Bash: basta clonar o repositório (não há compilação).

```bash
git clone https://github.com/LBCM-UFV/prokfinder-pipeline.git
cd prokfinder-pipeline
chmod +x main-prokfinder.sh modules/*.sh modules/plot-statistics.R   # se o bit de execução tiver se perdido
./main-prokfinder.sh --check                                         # confere as dependências
```

Para atualizar uma cópia já clonada: `git pull`. Os bancos de referência ficam fora do repositório
(`~/prokfinder-refdb`, ou `--refdb-dir`) e são instalados depois com `./main-prokfinder.sh setup`.

### Instalação via Docker

A imagem Docker traz todos os programas (NCBI Datasets, taxonkit, seqkit, Pyrodigal, CheckM2, BUSCO, skDER/CiDDER,
GTDB-Tk e R), sem precisar instalar dependências na máquina. Os **bancos** não ficam na imagem: os do Prokfinder
(`~/prokfinder-refdb`), o do CheckM2 (`CHECKM2DB`) e o pacote de referência do GTDB-Tk (`GTDBTK_DATA_PATH`) ficam na
máquina e são montados no contêiner.

```bash
git clone https://github.com/LBCM-UFV/prokfinder-pipeline.git
cd prokfinder-pipeline
docker build -f docker/Dockerfile -t prokfinder .                       # imagem completa (~15-20 min; ~14 GB em disco)
docker build -f docker/Dockerfile --build-arg WITH_GTDBTK=0 -t prokfinder .   # sem GTDB-Tk (menor)

# Uso: a mesma linha de comando do main-prokfinder.sh
docker/prokfinder-docker.sh --check
docker/prokfinder-docker.sh setup
export CHECKM2DB=/caminho/CheckM2_database/uniref100.KO.1.dmnd           # banco do CheckM2 (checkm2 database --download)
docker/prokfinder-docker.sh search-genomes -T 45:80 AND -O rod download-wf
```

O `docker/prokfinder-docker.sh` roda o contêiner com o seu usuário (os resultados ficam com você como dono) e monta,
nos mesmos caminhos: a pasta atual, a pasta dos bancos (`--refdb-dir` ou `~/prokfinder-refdb`), a pasta de resultados
(`-o`) e as pastas de `CHECKM2DB` e `GTDBTK_DATA_PATH` (repassa também a `NCBI_API_KEY`). No cluster, onde o Docker
não pode ser usado diretamente, a mesma imagem roda com o Singularity, depois de publicada num registro de imagens.
Dentro do contêiner, `CheckM2`, `BUSCO` e `GTDB-Tk` ficam em ambientes conda separados (as dependências deles são
incompatíveis), chamados por pequenos *wrappers*.

---

## Dependências

Não é preciso instalar tudo: para **buscar genomas** (`setup` e `search-genomes`) bastam o NCBI Datasets e o
taxonkit (usado pelo `setup`). Os programas de análise de genomas só são necessários para o `download-wf`.

`./main-prokfinder.sh --check` mostra o que está instalado (com a versão), o que falta (com o comando de instalação)
e o estado dos bancos de referência, separando as dependências da busca e as do pipeline completo.

**Busca de genomas** (`setup` e `search-genomes`):

| Programa | Instalação (conda) |
|---|---|
| `datasets` e `dataformat` (NCBI Datasets) | `conda install -c conda-forge ncbi-datasets-cli` |
| `taxonkit` | `conda install -c bioconda taxonkit` |
| bash ≥ 4 e ferramentas Unix: `awk`, `gawk`, `jq`, `curl`, `unzip`, `tar`, `gzip`, `md5sum` | sistema (`apt install ...`) ou conda |

**Pipeline completo** (`download-wf`), além dos anteriores:

| Programa | Instalação (conda) |
|---|---|
| `seqkit` | `conda install -c bioconda seqkit` |
| Python 3 + `pyrodigal` (predição de genes) | `pip install pyrodigal` ou `conda install -c bioconda pyrodigal` (vem junto com o `skder`) |
| `checkm2` + banco de dados | `conda install -c bioconda checkm2` e `checkm2 database --download` |
| `busco` (dispensável com `--skip-busco`) | `conda install -c bioconda busco` |
| `skder` (inclui o `cidder`) | `conda install -c bioconda skder` |
| `gtdbtk` + pacote de referência (recomendado: ANI dos genomas que não estão nas tabelas do GTDB; obrigatório com `--use-gtdbtk-tree`) | `conda install -c bioconda gtdbtk` |
| `Rscript` (opcional, só R base: gráficos finais) | `conda install -c conda-forge r-base` |

Os parâmetros do Prokfinder são dados **só por sub-comandos e flags**. As únicas variáveis de ambiente lidas são as
de configuração dos programas, que no cluster os `module load` já definem:

| Variável | Quando | Descrição |
|---|---|---|
| `CHECKM2DB` | `download-wf` | caminho do **arquivo** do banco do CheckM2 (`uniref100.KO.1.dmnd`) |
| `GTDBTK_DATA_PATH` | `download-wf` | pasta do pacote de referência do GTDB-Tk |
| `NCBI_API_KEY` | opcional | chave da API do NCBI (aumenta o limite de 5 para 10 requisições/s) |

Antes de cada execução, o Prokfinder confere só as dependências dos sub-comandos pedidos e lista todos os problemas
de uma vez. No cluster, o job de `run-on-cluster` carrega só os módulos `ncbi_datasets taxonkit` quando não há
`download-wf`.

---

## setup: bancos de referência

```bash
./main-prokfinder.sh setup                 # instala (ou atualiza) os bancos em ~/prokfinder-refdb
./main-prokfinder.sh setup --force-refdb   # baixa e instala tudo de novo
./main-prokfinder.sh setup --ncbi-complete-metadata   # opcional: também converte todos os metadados do NCBI para TSV
```

- Os bancos são baixados **uma única vez** e compartilhados por todas as execuções e pastas de resultados. O download
  leva horas; no cluster, use `run-on-cluster setup`.
- `search-genomes` e `download-wf` **não rodam sem o `setup`**: se os bancos não forem encontrados, o Prokfinder
  para com o erro `Reference databases not found in ~/prokfinder-refdb` e pede o `setup` (que pode vir no mesmo
  comando: `setup search-genomes ...`).
- Rodar `setup` de novo com os bancos já instalados não baixa nada. Para reinstalar é preciso `--force-refdb`.
- Com `--force-refdb`, cada banco é baixado e formatado numa **pasta temporária** (`~/prokfinder-refdb/.NCBI.new`,
  `.BACDIVE.new`, `.GTDB.new`) e só substitui o banco atual quando termina sem erro. Se o download falhar, o banco
  atual continua lá e funcionando; o mesmo comando retoma a instalação nova de onde parou, sem baixar de novo os
  bancos que já foram trocados (eles ficam anotados em `~/prokfinder-refdb/.force-refdb.state`, apagado no fim).
- Bancos instalados por uma versão **anterior** do Prokfinder, cujas tabelas tinham outro formato, também precisam do
  `setup`: ele **atualiza** as tabelas a partir dos arquivos já baixados, **sem baixar nada de novo** (a versão do
  formato fica registrada em `COMPLETED-DOWNLOAD.txt`, linha `format=`).
- Se duas execuções precisarem instalar os bancos ao mesmo tempo, a segunda espera a primeira terminar.
- `--ncbi-complete-metadata` (opcional): converte **todos** os campos dos metadados do NCBI (o
  `assembly_data_report.jsonl` inteiro, ~180 colunas) para TSV em `NCBI/<DOMÍNIO>-NCBI-METADATA/complete-metadata.tsv`.
  É lento (mais de 15 h só para Bacteria) e a busca não usa essa tabela; por isso fica desligado por padrão. Com os
  bancos já instalados, só essa tabela é gravada, a partir do JSONL já baixado.

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
./main-prokfinder.sh search-genomes --gtdb-taxon-file taxons.txt OR "(" --temp-growth 70: AND -i "hot spring","hydrothermal vent" ")"
./main-prokfinder.sh -o termofilos search-genomes --gtdb-taxon g__Thermus AND -c Brazil,Iceland
```

### Como escrever os filtros

1. Cada flag de seleção é **um filtro**; é preciso pelo menos um.
2. Vários valores da **mesma flag** são separados por **vírgula, sem espaço**, e combinados com **OR** (o filtro
   passa se qualquer um dos valores corresponder): `-c Brazil,Iceland`. Cada valor com espaço fica entre as suas próprias
   aspas duplas: `-i "hot spring","hydrothermal vent",geyser`.
3. **Entre dois filtros é obrigatório escrever `AND` ou `OR`**; sem isso o Prokfinder dá erro:

   ```bash
   -i "hot spring",geyser AND -c Brazil    # fonte termal ou gêiser, E do Brasil
   -i "hot spring",geyser OR -c Brazil     # fonte termal ou gêiser, OU do Brasil
   ```

4. `AND` é avaliado antes de `OR` (`A AND B OR C` = `(A AND B) OR C`). **Parênteses** agrupam filtros; no terminal
   eles precisam estar entre aspas ou com `\`: `-T 45:80 AND "(" -O rod OR -M anaerobe ")"`.
5. Arquivos de lista (`--gtdb-taxon-file`, `-A`, `-I`, `-C`): um valor por linha. Os `\r` (fim de linha do Windows)
   e os espaços no começo/fim das linhas são removidos, e linhas vazias são ignoradas. Em uma **tabela (TSV)** só a
   **primeira coluna** é lida. A primeira linha é tratada como **cabeçalho** e ignorada (com um aviso) quando não é
   um valor válido: um nome que não está na taxonomia do GTDB (`--gtdb-taxon-file`) ou um código de acesso inválido
   (`-A`). `-I` e `-C` não podem ter cabeçalho.

A **taxonomia usada pelo Prokfinder é a do GTDB**, tanto na busca (`--gtdb-taxon`) quanto na classificação dos
genomas (`download-wf`); não há filtros pela taxonomia do NCBI.

### Filtros

| Opção | Banco | Descrição |
|---|---|---|
| `--gtdb-taxon` / `--gtdb-taxon-file` | GTDB | nome de um táxon da linhagem GTDB (sem diferenciar maiúsculas). Com o prefixo do rank (`g__Thermus`, `s__Thermus aquaticus`), só aquele rank; sem prefixo (`Thermus`), qualquer rank. Um nome que não está nas tabelas do GTDB (ex.: `g_Thermus`, com um `_` só) é um valor inválido |
| `-a, --accession` / `-A, --accession-file` | NCBI | código de acesso **com versão** (`GCA_000005845.2`). GCA_ e GCF_ podem vir misturados (GCF_x.1 seleciona o par GenBank GCA_x.1). Todas as linhas inválidas são listadas no erro |
| `-i` / `-I` | NCBI | texto **ou expressão regular** (ERE do POSIX, sem diferenciar maiúsculas) procurado na fonte de isolamento (ex.: `"hot spring","hydrothermal vent"`, `"hot ?springs?"`, `"^soil$"`). Um texto simples também funciona; para procurar literalmente `( ) [ ] . * + ? \| ^ $`, use `\` antes (ex.: `"soil \(forest\)"`). Uma expressão inválida dá erro antes da busca |
| `-c` / `-C` | NCBI | país/local da localização geográfica (ex.: `Brazil,Iceland`) |
| `-g` | NCBI | `isolate_genome` ou `mag` |
| `-G` | NCBI | faixa de GC (%) da montagem (`min:max`) |
| `-T` / `-P` / `-S` | BacDive | temperatura, pH e NaCl **ótimos** (`min:max`, um dos lados pode ficar vazio) |
| `--temp-growth` / `--ph-growth` / `--salt-growth` | BacDive | valores em que houve crescimento |
| `-O, --cell-shape` | BacDive | formato celular: `rod`, `coccus`, `ovoid`, `curved`, `spiral`, `filament`, `pleomorphic`, `other` |
| `-M, --oxygen` | BacDive | tolerância a O₂: `aerobe`, `anaerobe`, `facultative`, `microaerophile` |
| `--gram-stain` | BacDive | `positive` ou `negative` |
| `--motility` | BacDive | células móveis: `yes` ou `no` |
| `--spore-formation` | BacDive | forma esporos: `yes` ou `no` |
| `--energy-source` | BacDive | `chemotroph` ou `phototroph` |
| `--electron-donor` | BacDive | `lithotroph` ou `organotroph` |
| `--carbon-source` | BacDive | `autotroph` ou `heterotroph` |

Outras opções: `-s genbank|refseq` (tabela do NCBI pesquisada; padrão `genbank`), `--force-selection` e
`--force-search` (abaixo).

### Filtros inválidos e `--force-search`

Todo valor ou filtro inválido é listado em **`SEARCH-GENOMES-ERRORS.tsv`**, na pasta de resultados (colunas `FILTER`,
`OPTION`, `VALUE`, `ERROR`). Exemplos: nome de táxon que não está na taxonomia do GTDB (`--gtdb-taxon`),
código de acesso inválido, expressão regular inválida (`-i`/`-I`), faixa inválida, classe
inexistente do BacDive, arquivo de lista que não existe ou vazio. O arquivo é apagado quando a busca não tem erros.

- **Sem `--force-search`** (padrão): a busca para e mostra os primeiros erros.
- **Com `--force-search`**: os valores inválidos são ignorados (com aviso) e a busca continua com os válidos. Um
  filtro que fica sem nenhum valor válido é **removido**, também da expressão: ele sai junto com o `AND` mais
  próximo (ou, se não houver, o `OR` mais próximo). Ex.: `F1 AND F2 OR F3` sem `F2` vira `F1 OR F3`;
  `F1 AND ( F2 OR F3 )` sem `F2` e `F3` vira `F1`. O log mostra a expressão final.

```bash
./main-prokfinder.sh search-genomes --gtdb-taxon-file cultivo.tsv OR --gtdb-taxon-file literatura.tsv OR -I fontes.txt --force-search
```

Erros de sintaxe da linha de comando (falta de `AND`/`OR`, parênteses, opção desconhecida) sempre param o Prokfinder.

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

Usa os genomas do último `search-genomes` da mesma pasta de resultados (`-o`). Etapas:

1. **Download** dos genomas (`datasets rehydrate`, a partir dos `fetch.txt` de Bacteria e Archaea). O domínio de cada
   genoma segundo o NCBI (o `fetch.txt` que o contém) vai para `GENOME-TAXONOMY/genome-domain-ncbi.tsv` e para a
   coluna `GENOME_DOMAIN_NCBI`; ele escolhe a linhagem do BUSCO (e, no futuro, as ferramentas de anotação).
2. **Filtro de contigs** pequenos (seqkit).
3. **Avaliação estrutural** (seqkit): número de contigs, N50 e bases N.
4. **Predição de genes** (Pyrodigal), **uma única vez por genoma**; as proteínas (`GENOME-GENES/<ACESSO>.faa`) são usadas
   pelo CheckM2, BUSCO, GTDB-Tk e CiDDER. Modo *single* com a tabela de tradução escolhida como no GTDB-Tk: tabela 4
   quando a densidade codificante com ela é mais de 5 pontos maior que com a tabela 11 e passa de 70 % (ex.:
   *Mycoplasma*); senão, tabela 11. Só genomas com menos de 100 kpb (curtos demais para treinar o Prodigal) usam o
   modo *meta*; os MAGs usam o modo *single*, como os isolados (ver abaixo).
5. **Avaliação biológica** sobre as proteínas: CheckM2 (`predict --genes`: completude e contaminação, que filtram) e
   BUSCO (`--mode proteins`, só nos genomas aprovados pelo CheckM2: ortólogos completos, *single-copy*, duplicados,
   fragmentados e ausentes; registrado, não filtra). Linhagem do BUSCO pelo domínio do NCBI: `bacteria_odb12` ou
   `archaea_odb12` (`odb10` com BUSCO < 6), baixada uma vez para `<refdb>/BUSCO`.
6. **Taxonomia GTDB**, em até três passos (ver [Classificação taxonômica](#classificação-taxonômica)): tabelas do
   GTDB → ANI contra os representantes do GTDB (`gtdbtk ani_rep`, pouca RAM) → árvore do GTDB-Tk (só com
   `--use-gtdbtk-tree`).
7. **Hierarquia** de diretórios Domínio/.../Espécie.
8. **Desreplicação** dentro de cada táxon: skDER (ANI, padrão) ou CiDDER (agrupamento das proteínas preditas).
9. Tabelas e gráficos finais.

| Etapa | Opção | Padrão | Descrição |
|---|---|---|---|
| Download | `-j, --download-threads <n>` | 10 | downloads simultâneos (1-30) |
| | `-z, --gzip` | | guarda os genomas comprimidos (`.fna.gz`) |
| Contigs | `--min-contig-length <pb>` | 1000 | contigs menores são removidos |
| Estrutura | `--max-contigs <n>` | 1000 | o genoma precisa ter menos contigs que isso |
| | `--max-n-bases <n>` | 100000 | o genoma precisa ter menos bases ambíguas (N) que isso |
| | `--min-n50 <pb>` | 5000 | N50 mínimo |
| Qualidade | `--min-completeness <%>` | 75 | completude mínima (CheckM2) |
| | `--max-contamination <%>` | 5 | contaminação máxima (CheckM2) |
| | `--skip-busco` | desligado | não roda o BUSCO |
| Taxonomia | `--use-gtdbtk-tree` | desligado | posiciona na árvore do GTDB os genomas que não estão nas tabelas nem foram classificados por ANI (**precisa de > 200 GB de RAM**: fila `high_mem`) |
| Desreplicação | `--derep-method <m>` | `skder` | `skder` (ANI) ou `cidder` (agrupamento das proteínas com CD-HIT, parâmetros padrão do CiDDER) |
| | `--derep-ani <%>` | 99,5 | skDER: identidade (ANI) acima da qual dois genomas do mesmo táxon são redundantes |
| | `--derep-mode <modo>` | `greedy` | skDER: `greedy` ou `dynamic` (fração alinhada fixa em 90 %) |
| | `--filter-mge` | desligado | mascara os elementos genéticos móveis (profagos e plasmídeos, previstos com o PhiSpy) antes de comparar os genomas |
| Re-execução | `--force-curation` | | descarta os genomas limpos, os genes preditos e todos os resultados de qualidade, taxonomia e desreplicação (os genomas baixados são mantidos) |

Os limites de qualidade são reaplicados a cada execução sobre os resultados já guardados: mudá-los não faz o CheckM2
rodar de novo.

**Por que os MAGs usam o modo *single*:** um MAG é o genoma de um único organismo, e o modo *single* treina o Prodigal
nele (prevê genes melhor que os modelos genéricos do modo *meta*, feito para sequências curtas demais para treinar ou
para misturas de organismos, como contigs de um metagenoma não separados em *bins*). É também o que o GTDB-Tk faz.
Só genomas com menos de 100 kpb vão para o modo *meta*.

### Classificação taxonômica

Cada genoma de alta qualidade recebe a taxonomia do GTDB pelo primeiro destes métodos que o classificar (coluna
`GTDB_CLASSIFICATION_METHOD`):

| Método | Como |
|---|---|
| `GTDB_TABLE` | o genoma está nas tabelas do GTDB do banco de referência (pelo mesmo acesso, pelo par GCA/GCF ou por outra versão) |
| `GTDBTK_ANI` | `gtdbtk ani_rep` (skani, sobre o genoma): o representante do GTDB mais próximo cumpre o critério de espécie do GTDB (ANI ≥ o raio de ANI da espécie, em geral 95 %, e fração alinhada ≥ 0,5). Usa pouca RAM; roda sempre que o GTDB-Tk e o `GTDBTK_DATA_PATH` existirem (senão, um aviso e o passo é pulado) |
| `GTDBTK_TREE` | só com `--use-gtdbtk-tree`: `identify --genes` → `align` → `classify --genes` (posicionamento na árvore, a partir das proteínas do Pyrodigal). Precisa de muita RAM |
| `NOT_CLASSIFIED` | nenhum dos anteriores: vai para `NON-CLASSIFIED-GENOMES-BY-GTDB.txt` e fica fora da hierarquia e da desreplicação |

`GENOME-TAXONOMY/taxonomy-classification.tsv` mostra, para cada genoma, o método, a fonte (versão do GTDB, ou versão
do GTDB-Tk e pacote de referência), o acesso encontrado na tabela do GTDB, o representante mais próximo, o ANI, a
fração alinhada, **os limites usados** (`ANI_THRESHOLD`, `AF_THRESHOLD`), o método da árvore e uma nota explicando
por que um genoma não foi classificado.

Os resultados do ANI e da árvore ficam em cache (`GENOME-TAXONOMY/GTDBTK-RESULTS/`). Assim, depois de uma execução
sem a árvore, **o mesmo comando com `--use-gtdbtk-tree` e a mesma pasta (`-o`)** só posiciona os genomas que ainda
estão em `NON-CLASSIFIED-GENOMES-BY-GTDB.txt` (o resto é reaproveitado), atualiza as tabelas e a lista, e segue com
a hierarquia e a desreplicação:

```bash
./main-prokfinder.sh run-on-cluster -o termofilos search-genomes --gtdb-taxon g__Thermus download-wf
./main-prokfinder.sh run-on-cluster --queue high_mem --mem 256gb --cpus 32 -o termofilos \
    search-genomes --gtdb-taxon g__Thermus download-wf --use-gtdbtk-tree
```

**CiDDER 1.3.6** (a versão do módulo `skder` do cluster) não consegue ler só proteomas (falha ao calcular o N50 dos
genomas); o **1.3.8** (o da imagem Docker) já lê. Quando a falha acontece, o Prokfinder roda o CiDDER de novo com os
genomas, e o próprio CiDDER prediz os genes (com a tabela 11), com um aviso no log e `CIDDER_INPUT = genomes` no
resumo da desreplicação.

---

## run-on-cluster: execução no cluster (PBS)

Em um cluster com o gerenciador de filas **PBS** (como o `cluster.dti.ufv.br` da UFV), o Prokfinder não deve rodar
no nó de login: no cluster da UFV, qualquer comando que passe de **10 minutos** no nó de login é interrompido. Por
isso, no cluster, coloque `run-on-cluster` antes dos outros sub-comandos:

```bash
cd ~/prokfinder-scripts                    # o job roda na pasta de onde foi submetido
./main-prokfinder.sh run-on-cluster setup  # resultados em ./prokfinder-results (fora do git: .gitignore)
./main-prokfinder.sh run-on-cluster --queue qtime --walltime 24:00:00 search-genomes -T 45:80 AND -O rod
./main-prokfinder.sh run-on-cluster --queue high_mem --mem 256gb --cpus 32 download-wf --use-gtdbtk-tree
```

O que acontece:

1. A linha de comando é conferida na hora (um filtro inválido dá erro antes de ir para a fila).
2. O módulo [`modules/submit-prokfinder-pbs.sh`](modules/submit-prokfinder-pbs.sh) escreve
   `<resultados>/PBS-JOBS/prokfinder_<data>_<hora>.pbs`, que carrega os módulos do cluster
   (`module load seqkit checkm2 busco ncbi_datasets gtdbtk taxonkit skder`, além do `r` para os gráficos; o Python do
   `skder` tem o `pyrodigal`; sem
   `download-wf`, só `module load ncbi_datasets taxonkit`) e executa o
   `main-prokfinder.sh` com a mesma linha de comando, sem o `run-on-cluster`, com `-p` e `--memory` iguais às CPUs e à
   RAM do job. O `prokfinder.log` registra o job, o nó, a fila, as CPUs, a RAM e o walltime.
3. O script é enviado à fila com `qsub` (se o `qsub` não existir, o Prokfinder sai com erro).
4. O andamento da execução é escrito **enquanto o job roda** em `<resultados>/prokfinder.log` (acompanhe com
   `tail -f <resultados>/prokfinder.log` e `qstat -u $USER`). A saída do próprio job,
   `<resultados>/PBS-JOBS/prokfinder_<data>_<hora>.log`, só aparece **quando o job termina** (também em caso de erro).
5. O código do job fica registrado em `<resultados>/PBS-JOBS/submitted-jobs.tsv`.

**Uma execução por pasta de resultados.** Dois processos escrevendo na mesma pasta corromperiam os resultados.
Antes de submeter ou executar, o Prokfinder procura na pasta (`-o`):

- jobs registrados em `PBS-JOBS/submitted-jobs.tsv` que ainda estão na fila ou rodando (`qstat`);
- uma execução fora do PBS (por exemplo, no nó de login), marcada em `<resultados>/.prokfinder-running` (apagado
  quando a execução termina, também com erro ou `qdel`).

Se encontrar uma, mostra qual é e pergunta se deve cancelá-la (`qdel <job>`, ou `kill <PID>` para um processo da
mesma máquina) e seguir com a nova execução. Sem confirmação (ou fora de um terminal, como dentro de um job), o
Prokfinder para e sugere cancelar a outra execução ou usar outra pasta de resultados (`-o <dir>`). Uma execução em
outra máquina sem PBS não pode ser conferida: se ela não estiver mais rodando, apague `.prokfinder-running`.
`run-on-cluster --dry-run` e `--test` não são conferidos (só escrevem o `.pbs`).

| Opção | Padrão | Descrição |
|---|---|---|
| `--queue <fila>` | `roteamento` | fila de submissão (ver abaixo) |
| `--mem <ram>` | `20gb` | memória RAM (`20gb`, `24000mb` ou só o número em GB); **mínimo 4 GB, e 8 GB com `download-wf`** (o CheckM2 sozinho passou de 3 GB com 166 genomas) |
| `--cpus <n>` | 20 | CPUs do job (**mínimo 4**) |

O PBS mata o job que usa mais CPUs ou RAM do que pediu. Por isso os programas usam as CPUs (`-p`) e a RAM
(`--memory`) do job, e o Prokfinder recusa um `-p` ou `--memory` maior que `--cpus` ou `--mem`.
| `--walltime <HH:MM:SS>` | `2160:00:00` (`72:00:00` na `qtime`) | tempo máximo do job |
| `--test` | | envia só um job curto (fila `qtime`, 15 min) que carrega os módulos e mostra a versão de cada programa |
| `--dry-run` | | só escreve o `.pbs`, sem submeter |

| Fila | Uso | Limites |
|---|---|---|
| `qtime` | testes e execuções curtas | walltime obrigatório, no máximo 72:00:00; até 24024 MB por job |
| `roteamento` | execuções comuns (padrão) | até 24024 MB por job; walltime até 2160:00:00 |
| `high_mem` | execuções com muita memória (ex.: `--use-gtdbtk-tree`) | no mínimo 25 GB; até 512 GB e 64 CPUs por usuário; 2 jobs rodando |

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
│   ├── BACTERIA-NCBI-METADATA/complete-metadata.tsv   # opcional (--ncbi-complete-metadata): todos os campos do NCBI
│   ├── ARCHAEA-NCBI-GENOMES/  ARCHAEA-NCBI-METADATA/
│   └── TAXDUMP/                                       # taxonomia NCBI (rank dos TaxIDs do BacDive, no setup)
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
├── SEARCH-GENOMES-ERRORS.tsv                          # valores e filtros inválidos da busca (só quando houver)
├── .prokfinder-running                                # execução em andamento (host, PID, job); apagado no fim
├── PBS-JOBS/                                          # scripts .pbs, logs e submitted-jobs.tsv (só com run-on-cluster)
└── OUTPUT/
    ├── genome-selection-query.txt
    ├── GENOME-METADATA/genome-metadata.tsv            # tabela única, com as colunas de todas as etapas
    ├── GENOME-FILES/
    │   ├── raw_genomes/<ACESSO>.fna                   # FASTA completo
    │   └── clean_genomes/<ACESSO>.fna                 # sem contigs pequenos, ordenados e renomeados
    ├── GENOME-GENES/                                  # <ACESSO>.faa/.ffn/.gff (Pyrodigal) e gene-prediction.tsv
    ├── GENOME-QUALITY/                                # genome-structure.tsv (seqkit), genome-quality.tsv (CheckM2, BUSCO)
    ├── GENOME-TAXONOMY/                               # genome-taxonomy.tsv, taxonomy-classification.tsv (método de
    │                                                  # cada genoma), NON-CLASSIFIED-GENOMES-BY-GTDB.txt, GTDBTK-RESULTS/
    │                                                  # e genome-domain-ncbi.tsv (domínio de cada genoma pelo NCBI)
    ├── GENOME-DIR-HIERARCHY/
    │   └── Domínio/Filo/Classe/Ordem/Família/Gênero/Espécie/
    │       ├── genomes-files-symlink.txt              # acessos e caminhos dos genomas do táxon
    │       └── GENOME-METADATA/genome-metadata.tsv    # metadados só desse táxon, com o representativo
    ├── GENOME-DEREPLICATION/                          # todos os resultados da desreplicação, num só lugar
    │   ├── dereplication-summary.tsv                  # uma linha por linhagem (ver abaixo)
    │   ├── representative-genomes.tsv                 # acesso e linhagem de cada genoma representativo
    │   ├── genome-dereplication.tsv                   # uma linha por genoma (representativo: yes/no)
    │   └── RUNS/Domínio/.../Espécie/                  # saída bruta do skDER/CiDDER de cada linhagem (cache)
    └── FINAL-RESULTS/
        ├── high-quality-selected-genomes.txt          # acessos + caminho do FASTA dos genomas finais
        ├── high-quality-selected-genomes-metadata.txt
        └── STATISTICS/
            ├── TABLES/   # funil de seleção, QC, genomas por táxon, tabela de todos os procariotos
            └── GRAPHS/   # gráficos (R)
```

`GENOME-DEREPLICATION/dereplication-summary.tsv` resume a desreplicação, uma linha por linhagem (o diretório-folha da
hierarquia):

| Coluna | Conteúdo |
|---|---|
| `TAXONOMIC_LINEAGE` / `TAXONOMIC_LEAF_RANK` | linhagem GTDB (`Bacteria;Pseudomonadota;...;Escherichia_coli`) e o rank da folha |
| `GENOMES_BEFORE_DEREPLICATION` / `GENOMES_AFTER_DEREPLICATION` | número de genomas antes e depois (representativos) |
| `REPRESENTATIVE_GENOMES` | acessos dos representativos, separados por `;` |
| `DEREPLICATION_SOFTWARE` | `skder`, `cidder` ou `NONE` (linhagem com um só genoma) |
| `SKDER_MODE` / `SKDER_ANI_CUTOFF` | modo (`greedy`/`dynamic`) e ANI de corte do skDER |
| `CIDDER_INPUT` | `proteomes` (proteínas do Pyrodigal) ou `genomes` (quando o CiDDER não lê os proteomas) |
| `MGE_FILTERED` / `MGE_SOFTWARE` | se os elementos genéticos móveis foram mascarados (`--filter-mge`) e com qual programa (PhiSpy) |
| `DEREPLICATION_STATUS` | `DEREPLICATED`, `SINGLE_GENOME` ou `FAILED` |

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
| `download-wf` | genomas já baixados; contigs, seqkit, genes, CheckM2, BUSCO e GTDB-Tk (ANI e árvore) por genoma; skDER/CiDDER por táxon (quando a lista de genomas e os parâmetros não mudaram) | `download-wf --force-curation` |

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

- O posicionamento na árvore do GTDB-Tk precisa de muita RAM (> 200 GB com o pacote completo), por isso só roda com
  `--use-gtdbtk-tree`; a classificação por ANI (`ani_rep`) usa pouca RAM e roda por padrão.
- Nem todo táxon do GTDB tem correspondência com o NCBI. O filtro `--gtdb-taxon` busca por nome e o vínculo entre
  os dois bancos é feito pelo código de acesso.
- A extração do BacDive foi reescrita em Bash a partir do script #1 do BactoTraits (originalmente em R); um intervalo
  como "5-30 °C" é guardado apenas pelos extremos (5:30).

## Próximos passos

- Construção de um banco de k-mers dos genomas finais (ex.: `kounta`).
- Análises do banco em R.
