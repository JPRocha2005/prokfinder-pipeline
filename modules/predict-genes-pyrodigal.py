#!/usr/bin/env python3
"""
predict-genes-pyrodigal.py - Prokfinder: gene prediction (Pyrodigal, the Python port of Prodigal)

Predicts the genes of each genome once, so that the next steps (CheckM2, BUSCO, GTDB-Tk and CiDDER) all use the
same predicted proteins. Called by predict-genes-pyrodigal.sh, which decides which genomes must be (re)predicted.

Gene calling mode and translation table (the heuristic of GTDB-Tk):
  - single mode ("-p single"): Pyrodigal is trained on the genome with translation table 11 and with table 4.
    Table 4 is used only when its coding density is more than 5 percentage points higher than with table 11 and
    above 70%; otherwise table 11 is used (e.g., Mycoplasma uses table 4).
  - meta mode ("-p meta", pre-trained models): only genomes shorter than 100 kbp (too short to train Prodigal).
    The meta mode is meant for sequences of several organisms or too short to train on; a MAG is the genome of a
    single organism, so it is predicted in single mode too (as GTDB-Tk does).

Input: a TSV file, one genome per line: <accession> <TAB> <genome FASTA, .fna or .fna.gz>.
Output, in --out-dir, for each genome: <ACCESSION>.faa (proteins, without the stop "*"), <ACCESSION>.ffn (gene
nucleotide sequences) and <ACCESSION>.gff; and a TSV table (--table) with one row per genome:
  ASSEMBLY_ACCESSION GENE_CALLING_MODE TRANSLATION_TABLE CODING_DENSITY_TABLE11 CODING_DENSITY_TABLE4
  PREDICTED_GENES PROTEIN_FILE
"""

import argparse
import gzip
import multiprocessing
import os
import sys

try:
    import pyrodigal
except ImportError:
    sys.exit("ERROR: the Python package 'pyrodigal' is not installed (pip install pyrodigal, or "
             "conda install -c bioconda pyrodigal).")

# GeneFinder replaced OrfFinder in Pyrodigal 3
GeneFinder = getattr(pyrodigal, "GeneFinder", None) or getattr(pyrodigal, "OrfFinder")

MIN_TRAINING_LENGTH = 100000   # shorter genomes are predicted in meta mode
TABLE4_DENSITY_GAIN = 0.05     # table 4 must raise the coding density by more than this...
TABLE4_MIN_DENSITY = 0.70      # ...and reach more than this


def read_fasta(path):
    """Returns the sequences of a FASTA file (plain or gzip) as a list of (id, sequence)."""
    opener = gzip.open if path.endswith(".gz") else open
    records, name, chunks = [], None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def coding_density(predictions, total_length):
    """Fraction of the genome covered by the predicted genes."""
    coding = sum(abs(gene.end - gene.begin) + 1 for genes in predictions for gene in genes)
    return coding / total_length if total_length else 0.0


def predict_single(records, total_length, table):
    finder = GeneFinder(meta=False)
    finder.train(*[seq for _, seq in records], translation_table=table)
    predictions = [finder.find_genes(seq) for _, seq in records]
    return predictions, coding_density(predictions, total_length)


def write_outputs(accession, records, predictions, out_dir):
    faa = os.path.join(out_dir, accession + ".faa")
    tmp = {ext: os.path.join(out_dir, "." + accession + ext + ".tmp") for ext in (".faa", ".ffn", ".gff")}
    with open(tmp[".faa"], "w") as f_faa, open(tmp[".ffn"], "w") as f_ffn, open(tmp[".gff"], "w") as f_gff:
        for index, ((name, _), genes) in enumerate(zip(records, predictions)):
            genes.write_translations(f_faa, sequence_id=name, include_stop=False)
            genes.write_genes(f_ffn, sequence_id=name)
            try:
                genes.write_gff(f_gff, sequence_id=name, header=(index == 0))
            except TypeError:   # older Pyrodigal: no "header" argument
                genes.write_gff(f_gff, sequence_id=name)
    # The files appear only when complete (an interrupted run leaves no truncated .faa)
    for ext, path in tmp.items():
        os.replace(path, os.path.join(out_dir, accession + ext))
    return faa


def predict_genome(task):
    accession, path, out_dir = task
    try:
        records = read_fasta(path)
        total_length = sum(len(seq) for _, seq in records)
        if total_length == 0:
            raise ValueError("empty genome")
        density11 = density4 = "NULL"

        if total_length < MIN_TRAINING_LENGTH:
            mode = "meta"
            finder = GeneFinder(meta=True)
            predictions = [finder.find_genes(seq) for _, seq in records]
            tables = {gene.translation_table for genes in predictions for gene in genes}
            table = ",".join(str(t) for t in sorted(tables)) or "NULL"
        else:
            mode = "single"
            predictions11, d11 = predict_single(records, total_length, 11)
            predictions4, d4 = predict_single(records, total_length, 4)
            density11, density4 = "%.4f" % d11, "%.4f" % d4
            if d4 - d11 > TABLE4_DENSITY_GAIN and d4 > TABLE4_MIN_DENSITY:
                table, predictions = "4", predictions4
            else:
                table, predictions = "11", predictions11

        faa = write_outputs(accession, records, predictions, out_dir)
        n_genes = sum(len(genes) for genes in predictions)
        return "\t".join([accession, mode, table, density11, density4, str(n_genes), os.path.abspath(faa)]), None
    except Exception as error:   # one bad genome must not stop the others
        return None, "%s: %s" % (accession, error)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", required=True, help="TSV: accession, genome FASTA")
    parser.add_argument("--out-dir", required=True, help="folder of the .faa/.ffn/.gff files")
    parser.add_argument("--table", required=True, help="TSV written with one row per genome predicted")
    parser.add_argument("--threads", type=int, default=1, help="genomes predicted at the same time")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    tasks = []
    with open(args.genomes) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0]:
                tasks.append((fields[0], fields[1], args.out_dir))

    failures = 0
    with open(args.table, "w") as out:
        out.write("ASSEMBLY_ACCESSION\tGENE_CALLING_MODE\tTRANSLATION_TABLE\tCODING_DENSITY_TABLE11\t"
                  "CODING_DENSITY_TABLE4\tPREDICTED_GENES\tPROTEIN_FILE\n")
        with multiprocessing.Pool(max(1, args.threads)) as pool:
            for done, (row, error) in enumerate(pool.imap_unordered(predict_genome, tasks), start=1):
                if row:
                    out.write(row + "\n")
                    out.flush()
                else:
                    failures += 1
                    print("WARNING: gene prediction failed for " + error, file=sys.stderr)
                if done % 100 == 0 or done == len(tasks):
                    print("  genes predicted for %d of %d genome(s)" % (done, len(tasks)), file=sys.stderr)
    print("pyrodigal %s: %d genome(s) predicted, %d failed" % (pyrodigal.__version__, len(tasks) - failures, failures),
          file=sys.stderr)


if __name__ == "__main__":
    main()
