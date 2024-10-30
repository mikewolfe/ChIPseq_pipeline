import sys
import numpy as np
import fasta as fa
import fastq as fq
import bed_utils as bed
import gzip

class QualityDistro(object):
    def __init__(self, max_length = 500):
        self.center_per_pos = np.zeros(max_length, dtype=int)
        self.spread_per_pos = np.zeros(max_length, dtype=float)

    def generate_quality(self, length, rng):
        quals = np.zeros(length, dtype=int)
        for loc, center, spread in zip(np.arange(length), self.center_per_pos, self.spread_per_pos):
            quals[loc] = int(rng.normal(center, spread))
        return "".join([chr(val + 33) for val in quals])

    def uniform_distro(self, center, spread):
        self.center_per_pos[:] = center
        self.spread_per_pos[:] = spread

class FragmentLengthDistro(object):
    def __init__(self, center, spread):
        self.center = center
        self.spread = spread

    def generate_length(self, rng):
        return int(rng.normal(self.center, self.spread))

class ReadLengthDistro(object):
    def __init__(self, center, spread):
        self.center = center
        self.spread = spread

    def generate_length(self, rng):
        return int(rng.normal(self.center, self.spread))

class ReadSampler(object):
    def __init__(self, length_distro, fragment_distro, quality_distro, fasta):
        self.length_distro = length_distro
        self.fragment_distro = fragment_distro
        self.quality_distro = quality_distro
        self.fasta = fasta

    def generate_read(self, name, chrm, loc, strand, rng, circ = False):
        name = "@%s"%(name)
        length = self.length_distro.generate_length(rng)
        this_chrm = self.fasta.pull_entry(chrm)
        chrm_length = len(this_chrm)
        if strand == "-":
            # have to shift register by one to match properly
            five_loc = max(0, loc-length+1)
            three_loc = loc+1
            seq = self.fasta.pull_entry(chrm).pull_seq(five_loc, three_loc, rc = True, circ = circ)
        else:
            five_loc = loc
            three_loc = min(loc + length, chrm_length-1)
            seq = self.fasta.pull_entry(chrm).pull_seq(five_loc, three_loc, rc = False, circ = circ)
        quality = self.quality_distro.generate_quality(len(seq), rng)

        return fq.FastqEntry(name = name, seq= seq, qual = quality)

    def generate_read_pair(self, chrm, loc, rng, frag_num, fragment_strand = ".", circ = False, chrm_to_num = None):

        fragment_length = self.fragment_distro.generate_length(rng)
        if fragment_strand == ".":
            fragment_strand = rng.choice(["-", "+"])
        five_loc = max(0, int(loc - fragment_length/2))
        three_loc = min(int(loc + fragment_length/2), len(self.fasta.pull_entry(chrm))-1)
        # for transcriptomes its sometimes better to convert to a numerical id
        if chrm_to_num:
            chrm_num = chrm_to_num[chrm]
        else:
            chrm_num = chrm
        shared_name = "SIM:%d:%s:0:0:%d:%d"%(frag_num, chrm_num, five_loc, three_loc)
        if fragment_strand == "-":
            R1 = self.generate_read(shared_name + " 1:N:0:ATCCGA", chrm, three_loc, fragment_strand, rng, circ = circ)
            R2 = self.generate_read(shared_name + " 2:N:0:ATCCGA", chrm, five_loc, "+", rng, circ = circ)
        else:
            R1 = self.generate_read(shared_name + " 1:N:0:ATCCGA", chrm, five_loc, "+", rng, circ = circ)
            R2 = self.generate_read(shared_name + " 2:N:0:ATCCGA", chrm, three_loc, "-", rng, circ = circ)
        return (R1, R2)
        
class LocSampler(object):

    def __init__(self, FastaFile = None, strands = ["+", "-"], initial_value = 0):
        self.weights = {}
        self.strands = strands
        if FastaFile is not None:
            for strand in self.strands:
                self.weights[strand] = {}
                for entry in FastaFile:
                    self.weights[strand][entry.chrm_name()] = np.zeros(len(entry), dtype=float) + initial_value
        self.probs = {}
        self.determine_probs()
        self.chrm_probs = {}
        self.determine_chrm_probs()
        self.strand_probs = {}
        self.determine_strand_probs()

    def determine_strand_probs(self):
        strand_sums = []
        for strand in self.strands:
            strand_sum = 0
            for chrm, weight_vec in self.weights[strand].items():
                strand_sum += np.sum(weight_vec)
            strand_sums.append(strand_sum)
        strand_probs = np.array(strand_sums)/np.sum(strand_sums)
        for strand, strand_prob in zip(self.strands, strand_probs):
            self.strand_probs[strand] = strand_prob

    def determine_chrm_probs(self):
        for strand in self.strands:
            self.chrm_probs[strand] = {}
            sum_weights = []
            chrms = []
            for chrm, weight_vec in self.weights[strand].items():
                chrms.append(chrm)
                sum_weights.append(np.sum(weight_vec))
            total = np.sum(sum_weights)
            chrm_probs = np.array(sum_weights) / total
            for chrm, chrm_prob in zip(chrms, chrm_probs):
                self.chrm_probs[strand][chrm] = chrm_prob

    def determine_probs(self):
        for strand in self.strands:
            self.probs[strand] = {}
            for chrm, chrm_weights in self.weights[strand].items():
                self.probs[strand][chrm] = chrm_weights/np.sum(chrm_weights)

    def add_enrichment_by_location(self, bedfile, footprint, max_enrichment = 16):
        strand = "."
        values = []
        for entry in bedfile:
            values.append(entry["score"])
        # scale values between 1 and max enrichment
        values = np.array(values)
        values = (max_enrichment - 1) *\
                ((values - np.min(values))/\
                (np.max(values)-np.min(values))) + 1

        half_footprint = int(footprint/2)
        for value, entry in zip(values, bedfile):
            center = int((entry["end"] - entry["start"])/2 + entry["start"])
            these_weights = self.weights[strand][entry["chrm"]]
            these_weights[(center-half_footprint):(center + half_footprint)] += value
            self.weights[strand][entry["chrm"]] = these_weights
        self.determine_probs()
        self.determine_chrm_probs()
        self.determine_strand_probs()

    def add_enrichment_by_region(self, bedfile):
        for entry in bedfile:
            self.weights[entry["strand"]][entry["chrm"]][entry["start"]:entry["end"]] += entry["score"]
        self.determine_probs()
        self.determine_chrm_probs()
        self.determine_strand_probs()

    def simulate(self, n, rng):
        locs = {}
        # which strands to pull from
        chosen_strands = rng.choice(self.strands, n, p = list(self.strand_probs.values()))
        strands, strand_n = np.unique(chosen_strands, return_counts = True)
        for this_strand, this_strand_n in zip(strands, strand_n):
            locs[this_strand] = {}
            # which chromosomes to pull from
            possible_chrms = list(self.chrm_probs[this_strand].keys())
            chrm_probs = list(self.chrm_probs[this_strand].values())
            chosen_chrms = rng.choice(possible_chrms, n, p = chrm_probs)
            chrm, chrm_n = np.unique(chosen_chrms, return_counts = True)
            # which locations to pull from for each chromosome
            for this_chrm, this_chrm_n in zip(chrm, chrm_n):
                locs[this_strand][this_chrm] = rng.choice(np.arange(len(self.probs[this_strand][this_chrm])), 
                        this_chrm_n, p = self.probs[this_strand][this_chrm])
        return locs

def ChIPseq_main(args):

    # input genome fasta
    infasta = args.infile
    # yaml file controlling the software
    out_prefix = args.out_pre
    n_fragments = args.num_fragments

    genome = fa.FastaFile()

    rng = np.random.default_rng(args.rng_seed)

    with open(infasta, mode = "r") as inf:
        genome.read_whole_file(inf)
    locations = LocSampler(genome, strands = ["."], initial_value = 1)

    # input locations for true peaks
    if args.peaks:
        inlocs = args.peaks

        inbed = bed.BedFile()
        inbed.from_bed_file(inlocs)

        enrichment = args.max_enrichment

        locations.add_enrichment_by_location(inbed, args.footprint, max_enrichment = enrichment)
    qual_distro = QualityDistro()
    qual_distro.uniform_distro(args.quality_mean, args.quality_std)
    fragment_distro = FragmentLengthDistro(args.fragment_length_mean, args.fragment_length_std)
    read_length_distro = ReadLengthDistro(args.read_length_mean, args.read_length_std)

    reads = ReadSampler(read_length_distro, fragment_distro, qual_distro, genome)

    out_R1 = gzip.open(out_prefix + "_R1.fastq.gz", mode ="wb")
    out_R2 = gzip.open(out_prefix + "_R2.fastq.gz", mode = "wb")
    
    all_locs = locations.simulate(n_fragments, rng)
    frag_num = 0
    for strand in all_locs.keys():
        for chrm in all_locs[strand].keys():
            for loc in all_locs[strand][chrm]:
                R1, R2 = reads.generate_read_pair(chrm, loc, rng, frag_num, strand) 
                frag_num += 1
                out_R1.write(str(R1).encode())
                out_R2.write(str(R2).encode())
    out_R1.close()
    out_R2.close()


def RNAseq_main(args):

    # input genome fasta
    infasta = args.infile
    # yaml file controlling the software
    out_prefix = args.out_pre
    n_fragments = args.num_fragments

    genome = fa.FastaFile()

    rng = np.random.default_rng(args.rng_seed)

    with open(infasta, mode = "r") as inf:
        genome.read_whole_file(inf)

    inbed = bed.BedFile()
    inbed.from_bed_file(args.scores)
    strands = np.unique([entry["strand"] for entry in inbed])
    locations = LocSampler(genome, strands = strands, initial_value = 0)
    locations.add_enrichment_by_region(inbed)

    qual_distro = QualityDistro()
    qual_distro.uniform_distro(args.quality_mean, args.quality_std)
    fragment_distro = FragmentLengthDistro(args.fragment_length_mean, args.fragment_length_std)
    read_length_distro = ReadLengthDistro(args.read_length_mean, args.read_length_std)

    reads = ReadSampler(read_length_distro, fragment_distro, qual_distro, genome)

    out_R1 = gzip.open(out_prefix + "_R1.fastq.gz", mode ="wb")
    out_R2 = gzip.open(out_prefix + "_R2.fastq.gz", mode = "wb")

    # get a chrm to number conversion
    chrm_to_num = {name: number for number, name in enumerate(reads.fasta.names)}
    
    all_locs = locations.simulate(n_fragments, rng)
    frag_num = 0
    for strand in all_locs.keys():
        for chrm in all_locs[strand].keys():
            for loc in all_locs[strand][chrm]:
                R1, R2 = reads.generate_read_pair(chrm, loc, rng, frag_num, strand, 
                        chrm_to_num = chrm_to_num) 
                frag_num += 1
                out_R1.write(str(R1).encode())
                out_R2.write(str(R2).encode())
    out_R1.close()
    out_R2.close()



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simulate fastq reads from a genome")
    subparsers = parser.add_subparsers(help = "Types of simulation")

    parent_parser = argparse.ArgumentParser(add_help = False)
    parent_parser.add_argument("infile", type=str, help="path to input fasta genome file")
    parent_parser.add_argument("out_pre", type =str, help = "output prefix for files")
    parent_parser.add_argument("--num_fragments", type = int, 
            help = "Number of fragments to simulate. Default = 1e6",
            default = 1000000)
    parent_parser.add_argument("--quality_mean", type = int,
            help = "Average quality score per base (Default = 36)", default = 36)
    parent_parser.add_argument("--quality_std", type = int,
            help = "Std quality score per base (Default = 1)", default = 1)
    parent_parser.add_argument("--fragment_length_mean", type = int,
            help = "Average fragment size (Default = 300)", default = 300)
    parent_parser.add_argument("--fragment_length_std", type = int,
            help = "Std fragment size (Default = 10)", default = 10)
    parent_parser.add_argument("--read_length_mean", type = int,
            help = "Average read length (Default = 150)", default = 150)
    parent_parser.add_argument("--read_length_std", type = int,
            help = "Std read length (Default = 0.0001)", default = 0.0001)
    parent_parser.add_argument("--rng_seed", type = int, help = "Seed for sampling")

    parser_RNAseq = subparsers.add_parser("Regions", parents = [parent_parser],
            help = "Simulate stranded or unstranded regions of enrichment in illumina data")
    parser_RNAseq.add_argument("--scores", type = str,
            help = "Bed file defining locations and values for weights")
    parser_RNAseq.set_defaults(func = RNAseq_main)

    parser_ChIPseq = subparsers.add_parser("ChIPseq", parents = [parent_parser],
            help = "Simulate illumina ChIPseq data")
    parser_ChIPseq.add_argument("--peaks", type = str, 
            help = "Bed file specifying areas of enrichment")
    parser_ChIPseq.add_argument("--max_enrichment", type = int, 
            help="Maximum fold enrichment above background", default = 16)
    parser_ChIPseq.add_argument("--footprint", type = int,
            help = "Estimated size of protein footprint for binding", default = 150)
    parser_ChIPseq.set_defaults(func = ChIPseq_main)
    

    args = parser.parse_args()
    args.func(args)
