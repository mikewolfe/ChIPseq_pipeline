# great module from deeptools that helps deal with bigwigs. I am only 
# adding support functions here for the module
import pyBigWig
import numpy as np
import arraytools

def write_arrays_to_bigwig(outfilename, arrays, chrm_dict, res = 1, dropNaNsandInfs = False):
    """
    Convert a set of arrays for each contig to a bigwig file
    
    Args:
        outfilename - (str) output file name
        arrays - (dict) a dictionary of numpy arrays
        bw - (dict) a dictionary of chromosome lengths
        res - resolution that the array is in
    Returns:
        Writes a bigwig file
    """
    with pyBigWig.open(outfilename, "w") as outf:
        chrm_list = [(key, chrm_dict[key]) for key in chrm_dict.keys()]
        outf.addHeader(chrm_list)
        for key in chrm_dict.keys():
            this_array = arrays[key]
            # make sure nans don't get added to the bw file
            if dropNaNsandInfs:
                the_finite = np.isfinite(this_array)
            else:
                the_finite = np.ones(len(this_array), dtype=bool)
            starts = np.arange(0, chrm_dict[key], res, dtype=np.int64)
            ends = np.arange(res, chrm_dict[key], res, dtype=np.int64)
            if len(ends) < len(starts):
                ends = np.append(ends, chrm_dict[key])
            names = np.array([key]*len(starts))
            outf.addEntries(names[the_finite], starts[the_finite], 
                    ends = ends[the_finite], values=this_array[the_finite])

def bigwig_to_arrays(bw, res = None):
    """
    Convert a bigwig to a dictionary of numpy arrays, one entry per contig
    
    Args:
        bw - pyBigWig object
        res - resolution that you want the array in
    Returns:
        arrays (dict) - a dictionary of numpy arrays
    """
    arrays = {}
    for chrm in bw.chroms().keys():
        arrays[chrm] = contig_to_array(bw, chrm, res)
    return arrays


def contig_to_array(bw, chrm, res = None):
    """
    Convert single basepair bigwig information to a numpy array
    

    Args:
        bw - a pyBigWig object
        chrm - name of chromosome you want 
        res - resolution you want data at in bp. 

    Returns:
        outarray - numpyarray at specified resolution

    """
    chrm_length = bw.chroms(chrm)
    # makes an array at 1 bp resolution
    out_array = bw.values(chrm, 0, chrm_length, numpy=True)
    if res:
        out_array = out_array[::res]
    return out_array


# Normalization/Transformation functions

def RobustZ_transform(arrays):
    """
    Transform data by the median and MAD of all contigs
    

    Args:
        arrays (dict) - dictionary of numpy arrays
    Returns:
        outarrays - dictionary of numpy arrays

    """
    from scipy import stats
    one_array = np.hstack(list(arrays.values()))
    median = np.nanmedian(one_array)
    mad = stats.median_abs_deviation(one_array, nan_policy='omit', scale='normal')
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], median, mad)
    return arrays

def Median_norm(arrays, pseudocount = 0):
    """
    Normalize data by the median value at all contigs
    

    Args:
        arrays (dict) - dictionary of numpy arrays
    Returns:
        outarrays - dictionary of numpy arrays

    """
    one_array = np.hstack(list(arrays.values()))
    one_array += pseudocount
    median = np.nanmedian(one_array)
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], -pseudocount, median)
    return arrays

def fixed_subtract(arrays, fixed_regions = None, res = 1, summary_func = np.nanmean):
    import bed_utils
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(fixed_regions)
    fixed_vals = []
    for region in inbed:
        fixed_vals.extend(arrays[region["chrm"]][region["start"]//res:region["end"]//res])
    fixed_vals = np.array(fixed_vals)
    subtract_val = summary_func(fixed_vals[np.isfinite(fixed_vals)])
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], subtract_val, 1)
    return arrays

def fixed_scale(arrays, fixed_regions = None, res = 1, summary_func = np.nanmean):
    import bed_utils
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(fixed_regions)
    fixed_vals = []
    for region in inbed:
        fixed_vals.extend(arrays[region["chrm"]][region["start"]//res:region["end"]//res])
    fixed_vals = np.array(fixed_vals)
    scale_val = summary_func(fixed_vals[np.isfinite(fixed_vals)])
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, scale_val)
    return arrays

def get_values_per_region(arrays, query_regions, res =1):
    import bed_utils
    inbed = bed_utils.BedFile()
    inbed.from_bed_file(query_regions)
    region_averages = []
    region_names = []
    for region in inbed:
        region_vals = arrays[region["chrm"]][region["start"]//res:region["end"]//res]
        region_average = np.nanmean(region_vals[np.isfinite(region_vals)])
        region_averages.append(region_average)
        region_names.append(region["name"])
    region_averages = np.array(region_averages)
    region_names = np.array(region_names)
    finite_vals = np.isfinite(region_averages)
    return(region_names[finite_vals], region_averages[finite_vals])

def query_subtract(arrays, res = 1, query_regions = None, number_of_regions = None, summary_func = np.nanmean):
    region_names, region_averages = get_values_per_region(arrays, query_regions, res)
    region_averages = np.array(region_averages)
    sort_indices = np.argsort(region_averages)
    region_names = np.array(region_names)
    vals = region_averages[sort_indices][:number_of_regions]
    regions = region_names[sort_indices][:number_of_regions]
    print("Bottom %s background regions"%number_of_regions)
    for region,val in zip(regions, vals):
        print(region + "\t" + "%s"%val)
    background_val = summary_func(vals)
    print("Background val: %s"%background_val)
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], background_val, 1)
    return arrays


def scale_max(arrays, num_top, res = 1):
    one_array = np.hstack(list(arrays.values()))
    one_array = one_array[np.isfinite(one_array)]
    ind = np.argpartition(one_array, -num_top//res)[-num_top//res:]
    scale_factor = np.nanmedian(one_array[ind])

    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, scale_factor)
    return arrays

def scale_region_max(arrays, number_of_regions, query_regions = None, res =1, summary_func = np.nanmean):
    region_names, region_averages = get_values_per_region(arrays, query_regions, res)
    region_averages = np.array(region_averages)
    sort_indices = np.argsort(region_averages)
    region_names = np.array(region_names)
    vals = region_averages[sort_indices][-number_of_regions:]
    regions = region_names[sort_indices][-number_of_regions:]
    print("Top %s background regions"%number_of_regions)
    for region,val in zip(regions, vals):
        print(region + "\t" + "%s"%val)
    scale_factor = summary_func(vals)
    print("Max val: %s"%scale_factor)
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, scale_factor)
    return arrays



def manipulate_main(args):

    summary_func_dict = { "mean": np.nanmean,
            "median": np.nanmedian,
            "sum": np.nansum }

    operation_dict={"RobustZ": RobustZ_transform, 
            "Median_norm": lambda x: Median_norm(x, args.pseudocount),
            "fixed_subtract": lambda x: fixed_subtract(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "fixed_scale" : lambda x : fixed_scale(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "scale_max": lambda x: scale_max(x, 1000, args.res),
            "query_scale": lambda x: scale_region_max(x, args.number_of_regions, args.query_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "query_subtract": lambda x: query_subtract(x, args.res, args.query_regions, args.number_of_regions, summary_func = summary_func_dict[args.summary_func]),
            "spike_scale": lambda x: fixed_scale(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func])}

    # read in file 
    inf = pyBigWig.open(args.infile)
    
    # convert to a dictionary of arrays
    arrays = bigwig_to_arrays(inf, res = args.res)
    
    # perform operation on arrays
    arrays = operation_dict[args.operation](arrays)

    # write out file
    write_arrays_to_bigwig(args.outfile, arrays, inf.chroms(), res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
    inf.close()

def read_multiple_bws(bw_files, res = 1):
    all_bws = {}
    open_fhandles = [pyBigWig.open(fname) for fname in bw_files]
    for fname, handle in zip(bw_files, open_fhandles):
        print(fname)
        all_bws[fname] = bigwig_to_arrays(handle, res = res)
    return (all_bws, open_fhandles)

def query_summarize_identity(all_bws, samp_names, samp_to_fname, inbed, res):
    outvalues = {fname: [] for fname in args.infiles}
    region_names = []
    coordinates = []
    for region in inbed:
        for fname in args.infiles:
            these_arrays = all_bws[fname]
            array_len = len(these_arrays[region["chrm"]])
            if region["strand"] == "-":
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
    
            these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
            outvalues[fname].extend(these_values)
        these_coordinates = list(range((left_coord//res)*res,\
                (right_coord//res)*res, res))
        coordinates.extend(these_coordinates)
        region_names.extend([region["name"]]*len(these_coordinates))

    with open(args.outfile, mode = "w") as outf:
        header = "region\tcoord\t%s"%("\t".join(samp_names))
        outf.write(header + "\n")
        for i, (region, coord) in enumerate(zip(region_names, coordinates)):
            values = "\t".join([str(outvalues[samp_to_fname[samp]][i]) for samp in samp_names])
            outf.write("%s\t%s\t%s\n"%(region, coord, values))

def query_summarize_single(all_bws, samp_names, samp_to_fname, inbed, res, summary_func = np.nanmean):
    outvalues = {fname: [] for fname in args.infiles}
    region_names = []
    for region in inbed:
        for fname in args.infiles:
            these_arrays = all_bws[fname]
            array_len = len(these_arrays[region["chrm"]])
            if region["strand"] == "-":
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
    
            these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
            this_summary = summary_func(these_values)
            outvalues[fname].append(this_summary)
        region_names.append(region["name"])

    with open(args.outfile, mode = "w") as outf:
        header = "region\t%s"%("\t".join(samp_names))
        outf.write(header + "\n")
        for i, region in enumerate(region_names):
            values = "\t".join([str(outvalues[samp_to_fname[samp]][i]) for samp in samp_names])
            outf.write("%s\t%s\n"%(region, values))

        
def query_main(args):
    import bed_utils

    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.regions)
    res = args.res
    all_bws, open_fhandles = read_multiple_bws(args.infiles, res = res)
    
    if args.samp_names:
        samp_to_fname = {samp_name : fname for fname, samp_name in zip(args.infiles, args.samp_names)}
        samp_names = args.samp_names
    else:
        samp_to_fname = {fname : fname for fname in args.infiles}
        samp_names = args.infiles

    summary_funcs = {'mean' : np.nanmean,
            'median': np.nanmedian,
            'max' : np.nanmax,
            'min' : np.nanmin}
    try:
        summary_func = summary_funcs[args.summary_func]
    except KeyError:
        KeyError("%s is not a valid option for --summary_func"%(args.summary_func))

    overall_funcs = {'identity' : query_summarize_identity,
        'single' : lambda x, y, z, a, b: query_summarize_single(x, y, z, a,b, summary_func)}
    try:
        overall_func = overall_funcs[args.summarize]
    except KeyError:
        KeyError("%s is not a valid option for --summarize"%(args.summarize))

    overall_func(all_bws, samp_names, samp_to_fname, inbed, res)
    
    for fhandle in open_fhandles:
        fhandle.close()


def summarize_main(args):

    res = args.res

    inf1 = pyBigWig.open(args.infile)
    
    arrays = bigwig_to_arrays(inf1, res = args.res)
    with open(args.outfile, mode = "w") as outf:
        outf.write("contig\tsignal\n")
        for chrm in arrays.keys():
            this_array = arrays[chrm]
            the_finite = np.isfinite(this_array)
            outf.write("%s\t%s\n"%(chrm, np.sum(this_array[the_finite])))

def compare_log2ratio(arrays1, arrays2):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = np.log2(arrays1[chrm]) - np.log2(arrays2[chrm])
    return out_array

def compare_divide(arrays1, arrays2):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = arrays1[chrm] / arrays2[chrm]
    return out_array

def compare_subtract(arrays1, arrays2):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = arrays1[chrm] - arrays2[chrm]
    return out_array

def compare_add(arrays1, arrays2):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = arrays1[chrm] + arrays2[chrm]
    return out_array

def compare_recipratio(arrays1, arrays2):
    out_array = compare_divide(arrays1, arrays2)
    for chrm in arrays1.keys():
        # figure out small values
        recip_values = out_array[chrm] < 1
        # take negative reciprocal and add 1
        out_array[chrm][recip_values] = -1/out_array[chrm][recip_values] + 1
        # subtract 1 from postive values
        out_array[chrm][~recip_values] = out_array[chrm][~recip_values] - 1

    return out_array

def compare_main(args):
    operation_dict ={"log2ratio": compare_log2ratio,
            "add": compare_add,
            "subtract": compare_subtract,
            "divide": compare_divide,
            "recipratio": compare_recipratio}
    #read in files
    inf1 = pyBigWig.open(args.infile1)
    inf2 = pyBigWig.open(args.infile2)

    arrays1 = bigwig_to_arrays(inf1, res = args.res)
    arrays2 = bigwig_to_arrays(inf2, res = args.res)

    # perform operation
    arrays_out = operation_dict[args.operation](arrays1, arrays2)


    # write out file
    write_arrays_to_bigwig(args.outfile, arrays_out, inf1.chroms(), \
            res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
    inf1.close()
    inf2.close()

 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Tools to work with BigWig files")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    # manipulate verb
    parser_manipulate = subparsers.add_parser("manipulate", help = "manipulate BigWig files")
    parser_manipulate.add_argument('infile', type=str, 
            help="file to convert from")
    parser_manipulate.add_argument('outfile', type=str, help="file to convert to")
    parser_manipulate.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_manipulate.add_argument('--operation', type=str, default=None,
            help="operation to perform before writing out file. \
            All operations, neccesitate conversion to array internally \
            options {'RobustZ', 'Median_norm', 'background_subtract', 'scale_max', \
            'query_subtract', 'query_scale'}")
    parser_manipulate.add_argument('--fixed_regions', type=str, default=None,
            help="bed file containing a fixed set of regions on which to average over.")
    parser_manipulate.add_argument('--query_regions', type=str, default=None,
            help="bed file containing regions to consider.")
    parser_manipulate.add_argument('--number_of_regions', type = int, default = 20,
            help = "number of regions to include in background or max")
    parser_manipulate.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_manipulate.add_argument('--pseudocount', default = 0, type = int,
            help = "Add value to all unmasked regions before normalization. Only\
                    applicable to the median normalization. Default = 0 ")
    parser_manipulate.add_argument('--summary_func', default = "mean", type = str,
            help = "For methods that use regions. What summary function do you want to use over \
                    said regions. Options: mean, median, sum \
                    Default = mean ")
    parser_manipulate.set_defaults(func=manipulate_main)


    # summarize verb
    parser_summarize = subparsers.add_parser("summarize", help = "Get total signal per contig")
    parser_summarize.add_argument('infile', type=str, 
            help="file to convert from")
    parser_summarize.add_argument('outfile', type=str, help="file to convert to")
    parser_summarize.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_summarize.set_defaults(func=summarize_main)


    # query verb
    parser_query = subparsers.add_parser("query", help = "Lookup values in BigWig files")
    parser_query.add_argument("outfile", type=str, help = "file to output to")
    parser_query.add_argument("infiles", nargs = "+", type=str, help = "input files")
    parser_query.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_query.add_argument('--regions', type=str, help = "regions to grab data from")
    parser_query.add_argument('--upstream', type=int, help = "bp upstream to add")
    parser_query.add_argument('--downstream', type = int, help = "bp downstream to add")
    parser_query.add_argument('--samp_names', type = str, nargs = "+", help = "sample names for each file")
    parser_query.add_argument('--summarize', type = str, help = "How to summarize data. 'identity' reports \
            each data point. 'single' gives a single number summary. Default = 'identity'", default = 'identity')
    parser_query.add_argument('--summary_func', type = str, help = "What function to use to summarize data when not using \
            'identity' summary. mean, median, max, min supported. Default = 'mean'", default = "mean")
    parser_query.set_defaults(func=query_main)

    # compare verb
    parser_compare = subparsers.add_parser("compare", help = "Calculate operations between two BigWig files")
    parser_compare.add_argument("infile1", type=str, help = "input file 1")
    parser_compare.add_argument("infile2", type=str, help = "input file 2")
    parser_compare.add_argument("outfile", type=str, help = "file to output to")
    parser_compare.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input files")
    parser_compare.add_argument('--operation', type=str, default="log2ratio",
            help="Default is log2ratio i.e. log2(infile1) - log2(infile2). Other \
                    options include: add, subtract, divide, recipratio (invert ratios less than 1. Center at 0)")
    parser_compare.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_compare.set_defaults(func=compare_main)
    
    args = parser.parse_args()
    args.func(args) 
