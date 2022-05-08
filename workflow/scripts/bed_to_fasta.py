import fasta
import bed_utils
import argparse

def determine_start_end(feature, chrm, args):
    if args.centered:
        center = round((feature["start"] + feature["end"])/2)
        start = center
        end = center
    else:
        start = feature["start"]
        end = feature["end"]

    if feature["strand"] == "-":
        final_start = max(start - args.downstream, 0)
        final_end = min(end + args.upstream, len(chrm))
        final_rc = True
    else:
        final_start = max(start - args.upstream, 0)
        final_end = min(end + args.downstream, len(chrm))
        final_rc = False
    return(final_start, final_end, final_rc)

def determine_header(chrm, start, end, strand, name, args):
    if args.meme_header:
        return ">%s:%d-%d"%(chrm, start+1, end)
    else:
        return ">%s;%d;%d;%s;%s"%(chrm, start, end, strand,name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert bed to fasta")
    parser.add_argument("infile", type=str, help="path to input file")
    parser.add_argument("infasta", type=str, help="path to input fasta")
    parser.add_argument("outfile", type=str, help="path to output file")
    parser.add_argument("--centered", action = "store_true", help = "make everything relative to region center")
    parser.add_argument("--upstream", type = int, help = "# of bp to pad upstream. - numbers cut into region", default = 0)
    parser.add_argument("--downstream", type = int, help = "# of bp to pad downstream. - numbers cut into region", default = 0)
    parser.add_argument("--min_size", type = int, help = "minimum sequence length. Filters anything shorter than this after padding is applied")
    parser.add_argument("--meme_header", action = "store_true", help = "make the fasta entry header conform to memes expectations")
    parser.add_argument("--include_unstranded_rc", action="store_true", 
            help="include reverse complement of unstranded features?")

    # Read in bed file
    args = parser.parse_args()
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.infile)

    # read in fasta file
    infasta = fasta.FastaFile()
    with open(args.infasta) as inf:
        infasta.read_whole_file(inf)
    outfasta = fasta.FastaFile()


    # convert beds to fastas
    for entry in inbed:
        this_chrm = infasta.pull_entry(entry["chrm"])
        # if you don't know the strand put both the complement and 
        # reverse complement
        this_start, this_end, this_rc = determine_start_end(entry, this_chrm, args) 
        if (args.min_size is not None):
            if (this_end - this_start) < args.min_size:
                continue
        if entry["strand"] == ".":
            new_fasta_c = fasta.FastaEntry()
            new_fasta_c.set_header(determine_header(entry["chrm"], this_start, this_end, "+", entry["name"], args))
            new_fasta_c.set_seq(this_chrm.pull_seq(this_start, this_end, rc = False))
            outfasta.add_entry(new_fasta_c)
            if args.include_unstranded_rc:
                new_fasta_rc = fasta.FastaEntry() 
                new_fasta_rc.set_header(determine_header(entry["chrm"], this_start, this_end, "-", entry["name"], args))
                new_fasta.rc.set_seq(this_chrm.pull_seq(this_start, this_end, rc = True))
                outfasta.add_entry(new_fasta_rc)

        else:
            new_fasta = fasta.FastaEntry()
            new_fasta.set_header(determine_header(entry["chrm"], this_start, this_end, entry["strand"], entry["name"], args))
            new_fasta.set_seq(this_chrm.pull_seq(this_start, this_end,rc = this_rc))
            outfasta.add_entry(new_fasta)

    with open(args.outfile, mode="w") as outf:
        outfasta.write(outf)
