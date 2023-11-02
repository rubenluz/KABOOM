#!/usr/bin/env python
#
# KABOOM (proKAryote BUSCO phylogenOMics)
# Utility script to construct phylogenies using BUSCO genes.
#
# 2023 Rúben Luz
# https://github.com/rubenluz/KABOOM.git
#

import argparse
import logging
import os
import sys
from os import listdir, chdir, mkdir, system
from os.path import abspath, basename, isdir, join
from time import gmtime, strftime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="Perform phylogenomic reconstruction using BUSCO genes")

    parser.add_argument("-i", "--input", type=str, help="Input directory containing assembled genomes or completed "
                                                        "BUSCO runs", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output directory to store results", required=True)
    parser.add_argument("-t", "--threads", type=str, help="Number of threads to use", required=True)
    parser.add_argument("-l", "--lineage", type=str,
                        help="Choose witch lineage to use by BUSCO. Default: bacteria_odb10)", default="bacteria_odb10")
    parser.add_argument("-bp", "--busco_extra_parameters", type=str, default="",
                        help="Extra parameters for BUSCO. You can use this to restart or force a busco analysis.")
    parser.add_argument("-psc", "--percent_single_copy", type=float, action="store", dest="psc", default=100.0,
                        help="BUSCO presence cut-off. BUSCOs that are complete and single-copy in at least [-psc] "
                             "percent of species will be included in the concatenated alignment [default=100.0]")
    parser.add_argument("-ap", "--alignment_program", type=str,
                        help="Choose witch program to use to align", default="mafft")
    parser.add_argument("-ma", "--mafft_algorithm", type=str,
                        help="Choose witch algorithm to use with mafft: fast, genafpair (E-INS-i), localpair (L-INS-i), "
                             "globalpair (G-INS-i) [default=globalpair]", default="globalpair")
    parser.add_argument("-am", "--analysis_mode", default="AA",
                        help="Choose to perform phylogenetic reconstruction using nucleotide (DNA) sequences or "
                             "aminoacids sequences (AA) or both (it will perform both analysis one after the other). "
                             "Default: AA")
    parser.add_argument("--trimal_strategy", type=str, action="store", dest="trimal_strategy", default="automated1",
                        help="trimal trimming strategy (automated1, gappyout, strict, strictplus) [default=automated1]")
    parser.add_argument("--missing_character", type=str, action="store", dest="missing_character",
                        help="Character to represent missing data [default='?']", default="?")
    parser.add_argument("--phylogenetics_only", action="store_true", help="It flags to only run the phylogenetic "
                                                                          "analysis. It expects to encounter BUSCO "
                                                                          "files")
    parser.add_argument("-j", "--jackknife", type=str, default=1000,
                        help="Choose the number of jackknife replicates to by done by IQ-TREE. Default: 1000")
    parser.add_argument("-bb", "--ufbootstrap", type=str, default=1000,
                        help="Choose the number of bootstrap replicates to by done by IQ-TREE. Default: 1000")

    args = parser.parse_args()

    logging.basicConfig(filename="log.txt", level=logging.INFO)
    print_message("Starting prokaryotic BUSCO Phylogenomics Pipeline")
    print_message("User provided arguments:", sys.argv)
    print_message("Parsed arguments:", vars(args))

    input_directory = abspath(args.input)
    working_directory = abspath(args.output)
    threads = args.threads
    lineage = args.lineage
    analysis_mode = args.analysis_mode
    jackknife = args.jackknife
    ufbootstrap = args.ufbootstrap
    busco_extra_parameters = args.busco_extra_parameters
    mafft_algorithm = args.mafft_algorithm
    fast = False

    if busco_extra_parameters != "":
        busco_extra_parameters = "-" + busco_extra_parameters

    # Check if the input directory exists
    if not isdir(input_directory):
        print_message("ERROR. Input directory", input_directory, "not found")
        sys.exit()

    # Check if the output directory already exists
    if not isdir(working_directory):
        mkdir(working_directory)
    else:
        with os.scandir(working_directory) as it:
            if any(it):
                print_message()
                print_message('The working directory was not empty. Old runs might cause conflicts...')
                print_message()

    # Check trimal parameter
    trimal_strategy = args.trimal_strategy.lower()

    if trimal_strategy not in ["automated1", "gappyout", "strict", "strictplus"]:
        print_message(
            "ERROR. Didn't understand --trimal_strategy parameter. Accepted options = automated1, gappyout, strict, strictplus")
        sys.exit()

    trimal_strategy = "-" + trimal_strategy

    # Check alignment program
    alignment = args.alignment_program.lower()
    if alignment not in ["mafft", "muscle"]:
        print_message(
            "ERROR. Didn't understand --alignment parameter. Accepted options = mafft, muscle")
        sys.exit()

    if alignment == "mafft":
        mafft_algorithm = args.mafft_algorithm.lower()
        if mafft_algorithm not in ["fast", "genafpair", "localpair", "globalpair"]:
            print_message(
                "ERROR. Didn't understand --mafft_algorithm. Accepted options = fast, genafpair (E-INS-i), localpair (L-INS-i), globalpair (G-INS-i)")
            sys.exit()
        if mafft_algorithm == "fast":
            fast = True

    if not args.phylogenetics_only:
        print_message(
            "Preparing for BUSCO. The script will check the file names por possible prohibit characters '=', '/' among others")

        print_message("Running BUSCO")
        busco_working_directory = working_directory + "/busco_output"
        if not isdir(busco_working_directory):
            mkdir(busco_working_directory)

        run_busco(args.input,
                  lineage,
                  args.output + "/busco_output",
                  threads,
                  busco_extra_parameters)

        chdir(busco_working_directory)
        print_message("Looking for BUSCO runs in", busco_working_directory)

    else:
        chdir(input_directory)
        print_message("Looking for BUSCO runs in", input_directory)

    busco_samples = []
    busco_sample_names = []

    for i in listdir("."):
        i = abspath(i)

        if isdir(i):
            for j in listdir(i):
                j = join(i, j)
                if isdir(j) and "run_" in j:
                    busco_samples.append(j)
                    busco_sample_names.append(basename(i))

    if len(busco_samples) == 0:
        print_message("ERROR. Didn't find any BUSCO directories")
        sys.exit()

    print_message("Found", len(busco_samples), "BUSCO directories.")
    # print_message()
    # for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
    #    print(busco_sample_name, busco_sample)
    # print_message()

    print_message("Identifying complete and single-copy BUSCO genes")

    all_nucleotide_buscos = set()
    all_protein_buscos = set()

    buscos_nucleotide = {}
    buscos_nucleotide_per_species = {}

    buscos_protein = {}
    buscos_protein_per_species = {}

    if "DNA" in analysis_mode:
        # Identify single-copy BUSCOs per species
        for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
            buscos_nucleotide_per_species[busco_sample_name] = []

            chdir(join(busco_sample, "busco_sequences", "single_copy_busco_sequences"))

            for f in listdir("."):
                if f.endswith(".fna"):
                    busco_name = f.rstrip(".fna")
                    all_nucleotide_buscos.add(busco_name)
                    if busco_name not in buscos_nucleotide:
                        buscos_nucleotide[busco_name] = []

                    record = SeqIO.read(f, "fasta")
                    new_record = SeqRecord(Seq(str(record.seq)), id=busco_sample_name, description="")

                    buscos_nucleotide_per_species[busco_sample_name].append(busco_name)
                    buscos_nucleotide[busco_name].append(new_record)
        # print_message()
        # print("Name\tComplete and single-copy BUSCO genes:")
        # print_message()
        # for busco_sample_name in busco_sample_names:
        # print(busco_sample_name, len(buscos_nucleotide_per_species[busco_sample_name]))

    print_message()
    if "AA" in analysis_mode:
        # Identify single-copy BUSCOs per species
        for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
            buscos_protein_per_species[busco_sample_name] = []

            chdir(join(busco_sample, "busco_sequences", "single_copy_busco_sequences"))

            for f in listdir("."):
                if f.endswith(".faa"):
                    busco_name = f.rstrip(".faa")
                    all_protein_buscos.add(busco_name)
                    if busco_name not in buscos_protein:
                        buscos_protein[busco_name] = []

                    record = SeqIO.read(f, "fasta")
                    new_record = SeqRecord(Seq(str(record.seq)), id=busco_sample_name, description="")

                    buscos_protein_per_species[busco_sample_name].append(busco_name)
                    buscos_protein[busco_name].append(new_record)

        # if "DNA" not in data_input:
        # print_message()
        # print("Name\tComplete and single-copy BUSCO genes:")
        # print_message()
        # for busco_sample_name in busco_sample_names:
        # print(busco_sample_name, len(buscos_protein_per_species[busco_sample_name]))
        # print_message()

    if "AA" in analysis_mode:
        number_of_buscos = buscos_protein
    else:
        number_of_buscos = buscos_nucleotide

    print_message(len(number_of_buscos), "unique BUSCO genes considered")
    chdir(working_directory)

    if args.psc == 100.0:
        print_message("Identifying BUSCO genes that are complete and single-copy in all species")
    else:
        print_message("Identifying BUSCO genes that are complete and single-copy in at least", args.psc,
                      "percent of species")
    print_message()
    single_copy_buscos = []

    for busco in number_of_buscos:
        percent_present = (len(number_of_buscos[busco]) / len(busco_sample_names) * 100)

        if percent_present >= args.psc:
            single_copy_buscos.append(busco)

    if len(single_copy_buscos) == 0:
        if args.psc == 100.0:
            print_message(
                "ERROR. Didn't identify any BUSCO genes that are complete and single-copy in all species")
            print_message(
                "ERROR. Didn't identify any BUSCO genes that are complete and single-copy in at least",
                args.psc, "percent of species")
        print_message(
            "You may want to adjust the --percent_single_copy parameter to allow for greater amounts of missing data "
            "if your dataset is patchy")
        sys.exit()

    if args.psc == 100.0:
        print_message("Identified", len(single_copy_buscos),
                      "BUSCO genes that are complete and single-copy in all species:")
        print_message("Identified", len(single_copy_buscos),
                      "BUSCO genes that are complete and single-copy in at least", args.psc,
                      "percent of species:")
    # print(",".join(single_copy_buscos))

    mkdir("supermatrix")
    chdir("supermatrix")
    print_message()

    if "DNA" in analysis_mode:
        mkdir("nucleotides")
        print_message("Writing nucleotide sequences to", join(working_directory, "supermatrix", "nucleotides"))
        for busco in single_copy_buscos:
            busco_records = buscos_nucleotide[busco]
            SeqIO.write(busco_records, join(working_directory, "supermatrix", "nucleotides", busco + ".fna"),
                        "fasta")
    if "AA" in analysis_mode:
        mkdir("proteins")
        print_message("Writing aminoacid sequences to", join(working_directory, "supermatrix", "proteins"))
        for busco in single_copy_buscos:
            busco_records = buscos_protein[busco]
            SeqIO.write(busco_records, join(working_directory, "supermatrix", "proteins", busco + ".faa"), "fasta")

    print_message()

    i = 0
    if "DNA" in analysis_mode:
        mkdir("nucleotide_alignments")
        print_message("Aligning nucleotide sequences using " + alignment + " with", threads, "parallel jobs to:",
                      join(working_directory, "supermatrix", "nucleotide_alignments"))
        for busco in single_copy_buscos:
            align_input_directory = join(working_directory, "supermatrix", "nucleotides", busco + ".fna")
            align_output_directory = join(working_directory, "supermatrix", "nucleotide_alignments", busco + ".aln")
            i = i + 1
            print_message("Aligning " + str(i) + " sequence out of " + str(len(single_copy_buscos)))
            if alignment == "mafft":
                run_mafft(threads, align_input_directory, align_output_directory, mafft_algorithm, fast)
            else:
                run_muscle(threads, align_input_directory, align_output_directory)

    i = 0
    if "AA" in analysis_mode:
        mkdir("protein_alignments")
        print_message("Aligning aminoacid sequences using " + alignment + " with", threads, "parallel jobs to:",
                      join(working_directory, "supermatrix", "protein_alignments"))
        for busco in single_copy_buscos:
            align_input_directory = join(working_directory, "supermatrix", "proteins", busco + ".faa")
            align_output_directory = join(working_directory, "supermatrix", "protein_alignments", busco + ".aln")
            i = i + 1
            print_message("Aligning " + str(i) + " sequence out of " + str(len(single_copy_buscos)))
            if alignment == "mafft":
                run_mafft(threads, align_input_directory, align_output_directory, mafft_algorithm, fast)
            else:
                run_muscle(threads, align_input_directory, align_output_directory)

    print_message()
    if "DNA" in analysis_mode:
        print_message("Trimming alignments using trimAL (" + trimal_strategy + ") to:",
                      join(working_directory, "supermatrix", "nucleotide_trimmed_alignments"))
        mkdir("nucleotide_trimmed_alignments")
        for busco in single_copy_buscos:
            trimal_input_directory = join(working_directory, "supermatrix", "nucleotide_alignments", busco + ".aln")
            trimal_output_directory = join(working_directory, "supermatrix", "nucleotide_trimmed_alignments",
                                           busco + ".trimmed.aln")
            run_trimal(trimal_input_directory, trimal_output_directory)

    if "AA" in analysis_mode:
        print_message("Trimming alignments using trimAL (" + trimal_strategy + ") to:",
                      join(working_directory, "supermatrix", "protein_trimmed_alignments"))
        mkdir("protein_trimmed_alignments")
        for busco in single_copy_buscos:
            trimal_input_directory = join(working_directory, "supermatrix", "protein_alignments", busco + ".aln")
            trimal_output_directory = join(working_directory, "supermatrix", "protein_trimmed_alignments",
                                           busco + ".trimmed.aln")
            run_trimal(trimal_input_directory, trimal_output_directory)

    if "DNA" in analysis_mode:
        print_message()
        print_message("Concatenating all trimmed nucleotide alignments")

        chdir("nucleotide_trimmed_alignments")

        nucleotide_alignments = {}
        nucleotide_partitions = []
        for busco_sample_name in busco_sample_names:
            nucleotide_alignments[busco_sample_name] = ""

        start = 1  # tracking lengths for partitions file

        # If psc == 100; we can just concatenate alignments (no missing data)
        if args.psc == 100:
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    for record in SeqIO.parse(alignment, "fasta"):
                        nucleotide_alignments[str(record.id)] += str(record.seq)
                    nucleotide_partitions.append(
                        [alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)
        else:
            # we need to handle missing data here
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    check_samples = busco_sample_names[:]

                    for record in SeqIO.parse(alignment, "fasta"):
                        nucleotide_alignments[str(record.id)] += str(record.seq)
                        check_samples.remove(str(record.id))

                    nucleotide_partitions.append(
                        [alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)

                    if len(check_samples) > 0:
                        # This means some species were missing this busco, fill alignment with missing character ("?" is default)
                        seq_len = len(str(record.seq))
                        for sample in check_samples:
                            nucleotide_alignments[sample] += (args.missing_character * seq_len)

        chdir(working_directory)
        chdir("supermatrix")

        alignment_records = []
        for sample in nucleotide_alignments:
            record = SeqRecord(Seq(nucleotide_alignments[sample]), id=sample, description="")
            alignment_records.append(record)

        print_message("Writing concatenated supermatrix alignment to fasta format:",
                      join(working_directory, "supermatrix", "supermatrix_dna.fasta"))
        SeqIO.write(alignment_records, "supermatrix_dna.fasta", "fasta")

        print_message("Writing concatenated supermatrix alignment to phylip format:",
                      join(working_directory, "supermatrix", "supermatrix_dna.phylip"))
        SeqIO.write(alignment_records, "supermatrix_dna.phylip", "phylip-relaxed")

        print_message("Writing partitions NEXUS file to:",
                      join(working_directory, "supermatrix", "supermatrix_dna.partitions.nex"))
        fo = open(join(working_directory, "supermatrix", "supermatrix_dna.partitions.nex"), "w")

        fo.write("#nexus\n")
        fo.write("begin sets;\n")

        for p in nucleotide_partitions:
            fo.write("\tcharset " + p[0] + " = supermatrix_dna.phylip: " + str(p[1]) + "-" + str(p[2]) + ";\n")

        fo.write("end;\n")

        fo.close()
        print_message("Concatenate files are ready!")
        print_message("Supermatrix alignment is ", str(len(nucleotide_alignments[sample])),
                      " nucleotides in length from " + str(len(single_copy_buscos)) + " BUSCO genes")
        print_message()

    if "AA" in analysis_mode:
        print_message("Concatenating all trimmed AMINOACIDS alignments")

        chdir("protein_trimmed_alignments")

        protein_alignments = {}
        protein_partitions = []
        for busco_sample_name in busco_sample_names:
            protein_alignments[busco_sample_name] = ""

        start = 1  # tracking lengths for partitions file

        # If psc == 100; we can just concatenate alignments (no missing data)
        if args.psc == 100:
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    for record in SeqIO.parse(alignment, "fasta"):
                        protein_alignments[str(record.id)] += str(record.seq)
                    protein_partitions.append(
                        [alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)
        else:
            # we need to handle missing data here
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    check_samples = busco_sample_names[:]

                    for record in SeqIO.parse(alignment, "fasta"):
                        protein_alignments[str(record.id)] += str(record.seq)
                        check_samples.remove(str(record.id))

                    protein_partitions.append(
                        [alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)

                    if len(check_samples) > 0:
                        # This means some species were missing this busco, fill alignment with missing character ("?" is default)
                        seq_len = len(str(record.seq))
                        for sample in check_samples:
                            protein_alignments[sample] += (args.missing_character * seq_len)

        chdir(working_directory)
        chdir("supermatrix")

        alignment_records = []
        for sample in protein_alignments:
            record = SeqRecord(Seq(protein_alignments[sample]), id=sample, description="")
            alignment_records.append(record)

        print_message("Writing concatenated supermatrix alignment to fasta format:",
                      join(working_directory, "supermatrix", "supermatrix_aminoacids.fasta"))
        SeqIO.write(alignment_records, "supermatrix_aminoacids.fasta", "fasta")

        print_message("Writing concatenated supermatrix alignment to phylip format:",
                      join(working_directory, "supermatrix", "supermatrix_aminoacids.phylip"))
        SeqIO.write(alignment_records, "supermatrix_aminoacids.phylip", "phylip-relaxed")

        print_message("Writing partitions NEXUS file to:",
                      join(working_directory, "supermatrix", "supermatrix_aminoacids.partitions.nex"))
        fo = open(join(working_directory, "supermatrix", "supermatrix_aminoacids.partitions.nex"), "w")

        fo.write("#nexus\n")
        fo.write("begin sets;\n")

        for p in protein_partitions:
            fo.write("\tcharset " + p[0] + " = supermatrix_aminoacids.phylip: " + str(p[1]) + "-" + str(p[2]) + ";\n")

        fo.write("end;\n")

        fo.close()

        print_message("Concatenate files are ready!")
        print_message("Supermatrix alignment is ", str(len(protein_alignments[sample])),
                      " aminoacids in length from " + str(len(single_copy_buscos)) + " BUSCO genes")
        print_message()

    if "DNA" in analysis_mode:
        print_message("Starting DNA phylogenetic analysis")
        iqtree_input_file = join(working_directory, "supermatrix", "supermatrix_dna.fasta")
        iqtree_output_directory = join(working_directory, "iqtree_dna_output")
        if not isdir(iqtree_output_directory):
            mkdir(iqtree_output_directory)

        chdir(iqtree_output_directory)
        run_iqtree(iqtree_input_file, "dna_tree", str(ufbootstrap), str(jackknife))

    if "AA" in analysis_mode:
        print_message()
        print_message("Starting Aminoacids phylogenetic analysis")
        iqtree_input_file = join(working_directory, "supermatrix", "supermatrix_aminoacids.fasta")
        iqtree_output_directory = join(working_directory, "iqtree_aminoacids_output")
        if not isdir(iqtree_output_directory):
            mkdir(iqtree_output_directory)

        chdir(iqtree_output_directory)
        run_iqtree(iqtree_input_file, "aminoacids_tree", str(ufbootstrap), str(jackknife))

    print_message()
    print_message("KABOOM finished")


def run_busco(busco_input_directory, lineage, busco_working_directory, threads, busco_extra_parameters):
    system(
        "busco -i " + busco_input_directory + " -l " + lineage + " -o " + busco_working_directory + " -m genome -c " + threads + " " + busco_extra_parameters)


def run_muscle(threads, muscle_input_directory, muscle_output_directory):
    system(
        "muscle -threads " + threads + " -align " + muscle_input_directory + " -output " + muscle_output_directory + " > /dev/null 2>&1")


def run_mafft(threads, mafft_input_directory, mafft_output_directory, algorithm, fast):
    if fast:
        system(
            "mafft --thread " + threads + " --retree 2 --quiet " + mafft_input_directory + " > " + mafft_output_directory)
    else:
        system(
            "mafft --thread " + threads + " --maxiterate 1000 --" + algorithm + " --quiet " + mafft_input_directory + " > " + mafft_output_directory)


def run_trimal(trimal_input_directory, trimal_output_directory):
    system("trimal -in " + trimal_input_directory + " -out " + trimal_output_directory + " -automated1 ")


def run_iqtree(iqtree_input_file, prefix, ufbootstrap, jackknife):
    system(
        "iqtree2 -T AUTO -s " + iqtree_input_file + " --prefix " + prefix + " --ufboot " + ufbootstrap + " --ufjack " + jackknife)


def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))
    logging.info(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))


if __name__ == "__main__":
    main()
