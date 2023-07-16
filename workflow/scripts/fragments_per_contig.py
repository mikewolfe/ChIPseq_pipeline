import sys


if __name__ == "__main__":
    contigs_dict = {}
    for line in sys.stdin:
        this_contig = line.split("\t")[2]
        this_counter = contigs_dict.get(this_contig, 0)
        this_counter += 1
        contigs_dict[this_contig] = this_counter
    sys.stdout.write("contig\tfragments\n")
    for (key, val) in contigs_dict.items():
        sys.stdout.write("%s\t%d\n"%(key, val))


