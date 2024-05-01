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


def write_arrays_to_bedgraph(outfilename, arrays, chrm_dict, res = 1, dropNaNsandInfs = False, header=False):
    """
    Convert a set of arrays for each contig to a bedgraph file
    
    Args:
        outfilename - (str) output file name
        arrays - (dict) a dictionary of numpy arrays
        chrm_dict - (dict) a dictionary of chromosome lengths
        res - resolution that the array is in
    Returns:
        Writes a bedgraph file
    """
    with open(outfilename, "w") as outf:
        if header:
            outf.write("chrom\tstart\tend\tvalue\n")
        chrm_list = [(key, chrm_dict[key]) for key in chrm_dict.keys()]
        for key in chrm_dict.keys():
            this_array = arrays[key]
            # make sure nans don't get added to the bedgraph file
            if dropNaNsandInfs:
                the_finite = np.isfinite(this_array)
            else:
                the_finite = np.ones(len(this_array), dtype=bool)
            starts = np.arange(0, chrm_dict[key], res, dtype=np.int64)
            ends = np.arange(res, chrm_dict[key], res, dtype=np.int64)
            if len(ends) < len(starts):
                ends = np.append(ends, chrm_dict[key])
            names = np.array([key]*len(starts))
            for name, start, end, value in zip(names[the_finite], starts[the_finite], ends[the_finite], this_array[the_finite]):
                outf.write("%s\t%s\t%s\t%s\n"%(name, start, end, value))

def write_arrays_to_wig(outfilename, arrays, chrm_dict, res = 1, dropNaNsandInfs = False, header=False):
    """
    Convert a set of arrays for each contig to a bedgraph file
    
    Args:
        outfilename - (str) output file name
        arrays - (dict) a dictionary of numpy arrays
        chrm_dict - (dict) a dictionary of chromosome lengths
        res - resolution that the array is in
    Returns:
        Writes a variableStep .wig file
    """
    with open(outfilename, "w") as outf:

        outf.write("track type=wiggle_0\n")
        chrm_list = [(key, chrm_dict[key]) for key in chrm_dict.keys()]
        for key in chrm_dict.keys():
            if key in arrays:
                this_array = arrays[key]
            else:
                continue
            # make sure nans don't get added to the bedgraph file
            if dropNaNsandInfs:
                the_finite = np.isfinite(this_array)
            else:
                the_finite = np.ones(len(this_array), dtype=bool)
            starts = np.arange(0, chrm_dict[key], res, dtype=np.int64)
            ends = np.arange(res, chrm_dict[key], res, dtype=np.int64)
            if len(ends) < len(starts):
                ends = np.append(ends, chrm_dict[key])
            names = np.array([key]*len(starts))
            # write header for chromosome
            outf.write("variableStep chrom=%s\n"%(key))
            # Simplest wig, one row per datum. Put datum on start of bin.
            for name, start, end, value in zip(names[the_finite], starts[the_finite], ends[the_finite], this_array[the_finite]):
                outf.write("%s %s\n"%(start+1, value))

def change_array_resolution(arrays, chrm_dict, summary_func = np.nanmean, res_from = 1, res_to = 1):
    """
    Change the resolution of an array by applying a summary function over bins
    """
    if res_to < res_from:
        raise ValueError("Can't go from lower res to higher res (%s to %s)"%(res_from, res_to))
    chrm_list = [(key, chrm_dict[key]) for key in chrm_dict.keys()]
    for key in chrm_dict.keys():
        this_array = arraytools.to_lower_resolution(arrays[key], res_to//res_from, summary_func)
        arrays[key] = this_array
    return arrays

def bigwig_to_arrays(bw, res = None, nan_to_zero = False):
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
        arrays[chrm] = contig_to_array(bw, chrm, res, nan_to_zero = nan_to_zero)
    return arrays


def contig_to_array(bw, chrm, res = None, nan_to_zero = False):
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
    if nan_to_zero:
        out_array[np.isnan(out_array)] = 0
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
    if median == 0:
        raise ValueError("Median was zero. Consider adding a pseudocount for median calculation")
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, median)
    return arrays

def smooth(arrays, wsize, kernel_type, edge, sigma = None, scaler = None):
    """
    Smooth data using a convolution with a kernel (flat or gaussian).

    Args:
        arrays (dict) - dictionary of numpy arrays
        wsize (int) - size of half the window in res (i.e if res is 5 and wsize
                      is 5 then the wsize in bp is 5*5)
        kernel_type (str) - "gaussian" or "flat" for the kernel. sd is set to span the window
        edge (str) - "mirror" or "wrap" for dealing with the edges
    Returns:
        outarrays - dictionary of numpy arrays
    """
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.smooth_1D(arrays[chrm], wsize, kernel_type, edge, sigma)
        if scaler is not None:
            arrays[chrm] = scaler(arrays[chrm])
    return arrays


def smooth_spline(arrays, knots = None, scaler = None):
    """
    Smooth data using a cubic Bspline

    Args:
        arrays (dict) - dictionary of numpy arrays
        knots (list) - locations in array coordinates of knots
    Returns:
        outarrays - dictionary of numpy arrays
    """
    for chrm in arrays.keys():
        these_knots = knots.get(chrm, None)
        arrays[chrm] = arraytools.Bspline_1D(arrays[chrm], these_knots)
        if scaler is not None:
            arrays[chrm] = scaler(arrays[chrm])
    return arrays

def savgol(arrays, wsize, polyorder, edge, scaler = None):
    """
    Smooth data using a Savtizky-Golay filter

    Args:
        arrays (dict) - dictionary of numpy arrays
        wsize (int) - size of half the window in res (i.e if res is 5 and wsize
                      is 5 then the wsize in bp is 5*5)
        polyorder (int) - order of polynomial to fit within each window. Can't
                          be larger than the total window size. 
                          Higher order smooths less.
        edge (str) - "mirror" or "wrap" for dealing with the edges
    Returns:
        outarrays - dictionary of numpy arrays
    """
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.savgol_1D(arrays[chrm], wsize, polyorder=polyorder, edge=edge)
        if scaler is not None:
            arrays[chrm] = scaler(arrays[chrm])
    return arrays


def smooth_GLMGam(arrays, max_knots=10, scaler = None):
    """
    Smooth data using a GLMGam.
    Assumes you are working with count data since it fits
    a Poisson model

    Args:
        arrays (dict) - dictionary of numpy arrays
        max_knots (int) - maximum number of knots to start with
    Returns:
        outarrays - dictionary of numpy arrays
    """
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.GLMGam_1D(arrays[chrm], max_knots)
        if scaler is not None:
            arrays[chrm] = scaler(arrays[chrm])
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
    if scale_val == 0:
        raise ValueError("Scale factor value was zero. Consider using different regions")
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
    if scale_factor == 0:
        raise ValueError("Scale factor value was zero. Consider using different regions")

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
    if scale_factor == 0:
        raise ValueError("Scale factor value was zero. Consider using different regions")
    for chrm in arrays.keys():
        arrays[chrm] = arraytools.normalize_1D(arrays[chrm], 0, scale_factor)
    return arrays

def read_knot_file(knot_file, res = 1):
    if knot_file is not None:
        with open(knot_file, mode = "r") as inf:
            knot_dict = {}
            for line in inf:
                chrm, knot = line.split("\t")
                these_knots = knot_dict.get(chrm, [])
                these_knots.append(int(knot)//res)
                knot_dict[chrm] = these_knots
    else:
        knot_dict = None

    return knot_dict


def manipulate_main(args):

    summary_func_dict = { "mean": np.nanmean,
            "median": np.nanmedian,
            "sum": np.nansum }
    scaler_dict = {"min_max": lambda x: arraytools.min_max_1D(x),
                    "fractional": lambda x: arraytools.fractional_1D(x)}
    if args.smooth_scaler is not None:
        smooth_scaler = scaler_dict[args.smooth_scaler]
    else:
        smooth_scaler = None

    operation_dict={"RobustZ": RobustZ_transform, 
            "Median_norm": lambda x: Median_norm(x, args.pseudocount),
            "fixed_subtract": lambda x: fixed_subtract(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "fixed_scale" : lambda x : fixed_scale(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "scale_max": lambda x: scale_max(x, 1000, args.res),
            "query_scale": lambda x: scale_region_max(x, args.number_of_regions, args.query_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "query_subtract": lambda x: query_subtract(x, args.res, args.query_regions, args.number_of_regions, summary_func = summary_func_dict[args.summary_func]),
            "spike_scale": lambda x: fixed_scale(x, args.fixed_regions, args.res, summary_func = summary_func_dict[args.summary_func]),
            "gauss_smooth": lambda x: smooth(x, args.wsize, kernel_type = "gaussian", edge = args.edge, sigma = args.gauss_sigma, min_max = args.smooth_minmax),
            "flat_smooth": lambda x: smooth(x, args.wsize, kernel_type = "flat", edge = args.edge, scaler = smooth_scaler),
            "savgol_smooth": lambda x: savgol(x, args.wsize, polyorder = args.savgol_poly, edge = args.edge, scaler = smooth_scaler),
            "Bspline_smooth": lambda x: smooth_spline(x, read_knot_file(args.smooth_knots, args.res), scaler = smooth_scaler),
            "GLMGam_smooth": lambda x: smooth_GLMGam(x, args.gam_max_knots, scaler = smooth_scaler),
            "scale_byfactor": lambda x: scale_byfactor(x, args.scalefactor_table, args.scalefactor_id, args.pseudocount)}

    # Extra logic if trying to downsample everything which requires both strands
    # for stranded data
    # if stranded also open minus strand file
    if args.operation == "downsample":

        # read in file 
        inf = pyBigWig.open(args.infile)
    
        # convert to a dictionary of arrays
        arrays = bigwig_to_arrays(inf, res = args.res, nan_to_zero = True)
        if args.minus_strand is not None:
            inf_minus = pyBigWig.open(args.minus_strand)
            arrays_minus = bigwig_to_arrays(inf_minus, res = args.res, nan_to_zero = True)
            dwn_plus, dwn_minus = downsample(arrays, args.scalefactor_table, args.scalefactor_id, arrays_minus)

            write_arrays_to_bigwig(args.outfile, dwn_plus, inf.chroms(), res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
            write_arrays_to_bigwig(args.minus_strand_out, dwn_minus, inf_minus.chroms(), res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
            inf.close()
            inf_minus.close()
        else:
            dwn_plus, dwn_minus = downsample(arrays, args.scalefactor_table, args.scalefactor_id)
            write_arrays_to_bigwig(args.outfile, dwn_plus, inf.chroms(), res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
            inf.close()
    else:

        # read in file 
        inf = pyBigWig.open(args.infile)
    
        # convert to a dictionary of arrays
        arrays = bigwig_to_arrays(inf, res = args.res)

        # perform operation on arrays
        arrays = operation_dict[args.operation](arrays)
        
        # write out file
        write_arrays_to_bigwig(args.outfile, arrays, inf.chroms(), res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
        inf.close()

def read_multiple_bws(bw_files, res = 1, res_to = None):
    all_bws = {}
    open_fhandles = [pyBigWig.open(fname) for fname in bw_files]
    for fname, handle in zip(bw_files, open_fhandles):
        print(fname)
        all_bws[fname] = bigwig_to_arrays(handle, res = res)
        if res_to is not None:
            all_bws[fname] = change_array_resolution(all_bws[fname], handle.chroms(), res_from = res, res_to =res_to)
    return (all_bws, open_fhandles)


def query_summarize_identity(all_bws, samp_names, samp_to_fname, inbed, res, gzip = False, coord = "absolute", minus_bws = None, samp_to_fname_minus = None, antisense = False):
    if minus_bws is not None:
        stranded = True
    else:
        stranded = False
    outvalues = {samp: [] for samp in samp_names}
    region_names = []
    coordinates = []
    chrm_names = []
    for region in inbed:
        for samp in samp_names:
            fname = samp_to_fname[samp]
            these_arrays = all_bws[fname]
            array_len = len(these_arrays[region["chrm"]])
            if region["strand"] == "-":
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
            if stranded:
                if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                    minus_fname = samp_to_fname_minus[samp]
                    these_minus_arrays = minus_bws[minus_fname]
                    these_values = these_minus_arrays[region["chrm"]][left_coord//res:right_coord//res]
                else:
                    these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
            else:
                these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]

            outvalues[samp].extend(these_values)
        these_coordinates = np.arange((left_coord//res)*res,\
                (right_coord//res)*res, res)
        if coord == "relative_start":
            if region["strand"] == "-":
                these_coordinates = ((region["end"])//res)*res - these_coordinates - res
            else:
                these_coordinates = these_coordinates - (region["start"]//res)*res
        elif coord == "relative_end":
            if region["strand"] == "-":
                these_coordinates = (region["start"]//res)*res - these_coordinates 
            else:
                these_coordinates = these_coordinates - ((region["end"])//res)*res + res
        elif coord == "relative_center":
            center_coord = (int((region["start"] + region["end"])/ 2)//res)*res
            if region["strand"] == "-":
                these_coordinates = center_coord - these_coordinates
            else:
                these_coordinates = these_coordinates - center_coord
        elif coord == "absolute":
            these_coordinates = these_coordinates
        else:
            raise ValueError("Coord must be one of 'relative_start', 'relative_center', 'relative_start', or 'absolute'. Not %s"%(coord))
        coordinates.extend(these_coordinates)
        region_names.extend([region["name"]]*len(these_coordinates))
        chrm_names.extend([region["chrm"]]*len(these_coordinates))

    header = "chrm\tregion\tcoord\t%s"%("\t".join(samp_names)) + "\n"
    values_func = lambda i: "\t".join([str(outvalues[samp][i]) for samp in samp_names])
    if gzip:
        with gzip.open(args.outfile, mode = "wb") as outf:
            outf.write(header.encode())
            for i, (chrm, region, coord) in enumerate(zip(chrm_names, region_names, coordinates)):
                line = "%s\t%s\t%s\t%s\n"%(chrm, region, coord, values_func(i))
                outf.write(line.encode())
    else:
        with open(args.outfile, mode = "w") as outf:
            outf.write(header)
            for i, (chrm, region, coord) in enumerate(zip(chrm_names, region_names, coordinates)):
                line = "%s\t%s\t%s\t%s\n"%(chrm, region, coord, values_func(i))
                outf.write(line)


def relative_polymerase_progression(array):
    return arraytools.weighted_center(array, only_finite = True, normalize = True) 

def traveling_ratio(array, res, wsize, peak_loc = None, upstream = 0, out = "ratio"):
    if peak_loc is not None:
        peak_loc = peak_loc + upstream
        if peak_loc < wsize:
            raise ValueError("Relative location (%s) must be > wsize (%s) for fixed traveling ratio."%(relative_location, wsize))
        return arraytools.traveling_ratio(array, wsize = wsize//res, peak = peak_loc//res, out = out)
    else:
        return arraytools.traveling_ratio(array, wsize = wsize//res, out = out)

def summit_loc(array, res, wsize, upstream):
    loc = arraytools.relative_summit_loc(array, wsize = wsize//res)
    return loc*res - upstream

def query_summarize_single(all_bws, samp_names, samp_to_fname, inbed, res, summary_func = np.nanmean, frac_na = 0.25, gzip = False, minus_bws = None, samp_to_fname_minus = None, antisense = False):
    if minus_bws is not None:
        stranded = True
    else:
        stranded = False
    outvalues = {samp: [] for samp in samp_names}
    region_names = []
    chrm_names = []
    for region in inbed:
        for samp in samp_names:
            fname = samp_to_fname[samp]
            these_arrays = all_bws[fname]
            array_len = len(these_arrays[region["chrm"]])
            if region["strand"] == "-":
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)

            if stranded:
                if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                    minus_fname = samp_to_fname_minus[samp]
                    these_minus_arrays = minus_bws[minus_fname]
                    these_values = these_minus_arrays[region["chrm"]][left_coord//res:right_coord//res]
                else:
                    these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
            else:
                these_values = these_arrays[region["chrm"]][left_coord//res:right_coord//res]
    
            # add a filter for regions that have high amounts of nans
            if (np.sum(np.isnan(these_values)) / len(these_values)) < frac_na:
                if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                    this_summary = summary_func(these_values[::-1])
                else:
                    this_summary = summary_func(these_values)
            else:
                this_summary = np.nan
            outvalues[samp].append(this_summary)
        region_names.append(region["name"])
        chrm_names.append(region["chrm"])

    header = "chrm\tregion\t%s"%("\t".join(samp_names)) + "\n"
    values_func = lambda i: "\t".join([str(outvalues[samp][i]) for samp in samp_names])
    if gzip:
        with gzip.open(args.outfile, mode = "wb") as outf:
            outf.write(header.encode())
            for i, (chrm, region) in enumerate(zip(chrm_names, region_names)):
                line = "%s\t%s\t%s\n"%(chrm, region, values_func(i))
                outf.write(line.encode())
    else:

        with open(args.outfile, mode = "w") as outf:
            outf.write(header)
            for i, (chrm, region) in enumerate(zip(chrm_names, region_names)):
                line = "%s\t%s\t%s\n"%(chrm, region, values_func(i))
                outf.write(line)


def gini_coefficient(data, unbiased = False):
    '''
    Calculate the gini coefficient for a 1D set of data
    Arguments:
        data - 1D numpy array of data
        unbiased - boolean to determine if unbiased estimator should be used
    References:
        https://en.wikipedia.org/wiki/Gini_coefficient
        https://www.statsdirect.com/help/nonparametric_methods/gini_coefficient.htm
        https://stackoverflow.com/questions/48999542/more-efficient-weighted-gini-coefficient-in-python
    Formula:
        Standard:
        $G = \frac{1}{n} (n + 1 - 2 \frac{\sum_{i=1}^n (n + 1 - i)y_i}{\sum_{i=1}^n y_i})$
        Unbiased:                                                                
        $G = \frac{1}{n-1} (n + 1 - 2 \frac{\sum_{i=1}^n (n + 1 - i)y_i}{\sum_{i=1}^n y_i})$
    '''
    data_sorted = np.sort(data)
    n = len(data)
    # the sum of the cumulative sum vector gives the numerator when the data
    # is sorted in ascending order. The float here prevents integer overflows
    cumulative = np.cumsum(data_sorted, dtype = float)
    # the last value in a cumulative vector is the total sum
    stat = (n + 1 - 2 * np.sum(cumulative)/ cumulative[-1])
    if unbiased:
        stat = stat/(n-1)
    else:
        stat = stat/n
    return stat


def query_main(args):
    import bed_utils

    inbed = bed_utils.BedFile()
    inbed.from_bed_file(args.regions)
    res = args.res
    all_bws, open_fhandles = read_multiple_bws(args.infiles, res = args.res, res_to = args.res_to)
    if args.gzip:
        import gzip
    
    if args.samp_names:
        samp_to_fname = {samp_name : fname for fname, samp_name in zip(args.infiles, args.samp_names)}
        samp_names = args.samp_names
    else:
        samp_to_fname = {fname : fname for fname in args.infiles}
        samp_names = args.infiles

    if args.minus_strand:
        if len(args.minus_strand) != len(args.infiles):
            raise ValueError("Amount of files for plus and minus strand do not match. Need one file for each strand")
        if args.samp_names is None or len(args.samp_names) != len(args.infiles):
            raise ValueError("Need to give one sample name for each file in stranded analyses")
        all_minus_bws, open_minus_fhandles = read_multiple_bws(args.minus_strand, res = res)
        minus_stfn = {samp_name : fname for fname, samp_name in zip(args.minus_strand, args.samp_names)}
    else:
        all_minus_bws = None
        minus_stfn = None


    if args.res_to is not None:
        res = args.res_to
    else:
        res = args.res

    summary_funcs = {'mean' : np.nanmean,
            'median': np.nanmedian,
            'max' : np.nanmax,
            'min' : np.nanmin,
            'RPP' : relative_polymerase_progression,
            'TR' : lambda array: traveling_ratio(array, res, args.wsize, args.TR_A_center, args.upstream, out = "ratio"),
            'TR_A': lambda array: traveling_ratio(array, res, args.wsize, args.TR_A_center, args.upstream, out = "A") ,
            'TR_B': lambda array: traveling_ratio(array, res, args.wsize, args.TR_A_center, args.upstream, out = "B"),
            # kept for compatibility
            'TR_fixed' : lambda array: traveling_ratio(array, res, args.wsize, args.upstream, args.TR_A_center, out = "ratio"),
            'summit_loc': lambda array: summit_loc(array, res, args.wsize, args.upstream),
            'Gini': lambda array: gini_coefficient(array)}
    try:
        summary_func = summary_funcs[args.summary_func]
    except KeyError:
        KeyError("%s is not a valid option for --summary_func"%(args.summary_func))

    overall_funcs = {'identity' : lambda bws, names, smtofn, bed, res, gzip: \
            query_summarize_identity(bws, names, smtofn, bed, res, gzip, \
            args.coords, all_minus_bws, minus_stfn, args.antisense),
        'single' : lambda bws, names, smtofn, bed, res, gzip: \
                query_summarize_single(bws, names, smtofn, bed, res, \
                summary_func, args.frac_na, gzip, all_minus_bws, minus_stfn, args.antisense)}
    try:
        overall_func = overall_funcs[args.summarize]
    except KeyError:
        KeyError("%s is not a valid option for --summarize"%(args.summarize))

    overall_func(all_bws, samp_names, samp_to_fname, inbed, res, gzip)
    
    for fhandle in open_fhandles:
        fhandle.close()
    if args.minus_strand:
        for fhandle in open_minus_fhandles:
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

def compare_log2ratio(arrays1, arrays2, pseudocount = 0):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = np.log2(arrays1[chrm] + pseudocount) - np.log2(arrays2[chrm]+pseudocount)
    return out_array

def compare_divide(arrays1, arrays2, pseudocount = 0):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = (arrays1[chrm]+pseudocount) / (arrays2[chrm]+pseudocount)
    return out_array

def compare_subtract(arrays1, arrays2, pseudocount=0):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = arrays1[chrm] - arrays2[chrm]
    return out_array

def compare_add(arrays1, arrays2, pseudocount = 0):
    out_array = {}
    for chrm in arrays1.keys():
        out_array[chrm] = arrays1[chrm] + arrays2[chrm]
    return out_array

def compare_recipratio(arrays1, arrays2, pseudocount = 0):
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
    arrays_out = operation_dict[args.operation](arrays1, arrays2, args.pseudocount)


    # write out file
    write_arrays_to_bigwig(args.outfile, arrays_out, inf1.chroms(), \
            res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
    inf1.close()
    inf2.close()

def check_same_chromosomes(array_list):
    chrms = array_list[0].keys()
    out = True
    for array in array_list:
        if set(array) != set(chrms):
            out = False
    return out

def open_multiple_bigwigs(file_list):
    out_handles = []
    for inf in file_list:
        out_handles.append(pyBigWig.open(inf))
    return out_handles

def convert_bigwigs_to_arrays(handle_list, res = 5):
    out_arrays = []
    for handle in handle_list:
        out_arrays.append(bigwig_to_arrays(handle, res = res))
    return out_arrays
        
def close_multiple_bigwigs(handle_list):
    for handle in handle_list:
        handle.close()

def multicompare_within_group(array_list, function = np.nanmean):
    # check to make sure the arrays list is ok
    if not check_same_chromosomes(array_list):
        raise ValueError("Can't merge bigwigs with different chromosomes")
    out_array = {}
    for chrm in array_list[0].keys():
        chrm_list = [array[chrm] for array in array_list]
        out_array[chrm] = function(chrm_list, axis = 0)
    return out_array


def multicompare_main(args):
    within_operation_dict = {"mean": np.nanmean,
            "median": np.nanmedian,
            "max": np.nanmax,
            "min": np.nanmin}
    between_operation_dict ={"log2ratio": compare_log2ratio,
            "add": compare_add,
            "subtract": compare_subtract,
            "divide": compare_divide,
            "recipratio": compare_recipratio}
    #read in files
    groupA_handles = open_multiple_bigwigs(args.groupA)
    groupA_arrays = convert_bigwigs_to_arrays(groupA_handles, res = args.res)

    # convert to one set of arrays
    groupA_array = multicompare_within_group(groupA_arrays, function = within_operation_dict[args.within_operation])

    if args.groupB is not None:
        groupB_handles = open_multiple_bigwigs(args.groupB)
        groupB_arrays = convert_bigwigs_to_arrays(groupB_handles, res = args.res)
        groupB_array = multicompare_within_group(groupB_arrays, function = within_operation_dict[args.within_operation])


        # perform between operation
        arrays_out = between_operation_dict[args.between_operation](groupA_array, groupB_array)
        close_multiple_bigwigs(groupB_handles)
    else:
        arrays_out = groupA_array

    # write out file
    write_arrays_to_bigwig(args.outfile, arrays_out, groupA_handles[0].chroms(), \
            res = args.res, dropNaNsandInfs = args.dropNaNsandInfs)
    close_multiple_bigwigs(groupA_handles)

def single_num_summary(arrays, function, contigs = None, expected_locs = None, res = None):
    import bed_utils
    if contigs is None:
        contigs = arrays.keys()
    if expected_locs is not None:
        out_array = []
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(expected_locs)
        for region in inbed:
            this_chrm = region["chrm"]
            array_len = len(arrays[this_chrm])
            left_coord = max(region["start"] - args.upstream, 0)
            right_coord = min(region["end"] + args.downstream, array_len * res)
            out_array.append(arrays[this_chrm][left_coord//res:right_coord//res])
        out_array = np.concatenate(out_array) 
    else:
        out_array = np.zeros(np.sum([arrays[chrm].size for chrm in contigs]), float)
        i = 0
        for chrm in contigs:
            arrsize = arrays[chrm].size
            out_array[i:i+arrsize] = arrays[chrm][:]
            i = i + arrsize
    return function(out_array)


def stranded_single_num_summary(arrays_plus, arrays_minus, function, contigs = None, expected_locs = None, res = None, antisense = False):
    import bed_utils
    if contigs is None:
        contigs = arrays_plus.keys()

    if expected_locs is not None:
        out_array = []
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(expected_locs)
        for region in inbed:
            this_chrm = region["chrm"]
            array_len = len(arrays_plus[this_chrm])
            if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                left_coord = max(region["start"] - args.downstream, 0)
                right_coord = min(region["end"] + args.upstream, array_len * res)
                out_array.append(arrays_minus[this_chrm][left_coord//res:right_coord//res])
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
                out_array.append(arrays_plus[this_chrm][left_coord//res:right_coord//res])
        out_array = np.concatenate(out_array)
    else:
        out_array = np.zeros(np.sum([arrays_plus[chrm].size for chrm in contigs])*2, float)
        i = 0
        for chrm in contigs:
            arrsize = arrays_plus[chrm].size
            out_array[i:i+arrsize] = arrays_plus[chrm][:]
            i = i + arrsize
        for chrm in contigs:
            arrsize = arrays_minus[chrm].size
            out_array[i:i+arrsize] = arrays_minus[chrm][:]
            i = i + arrsize

    return function(out_array)

def get_TMM_factor(sample_array, reference_array, total_sa, total_ra, spike_contigs, minus_sa = None, minus_ra = None, expected_locs = None, res = None, antisense = False, M_trim = 0.30, A_trim = 0.05, pseudocount = 0, A_cutoff = -1e10):
    from scipy import stats
    import bed_utils
    # turn everything into one big array
    sa = []
    ra = []
    if expected_locs is not None:
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(expected_locs)
        for region in inbed:
            this_chrm = region["chrm"]
            array_len = len(sample_array[this_chrm])
            # stranded locations
            if minus_sa is not None:
                if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                    left_coord = max(region["start"] - args.downstream, 0)
                    right_coord = min(region["end"] + args.upstream, array_len * res)
                    sa.append(minus_sa[this_chrm][left_coord//res:right_coord//res])
                    ra.append(minus_ra[this_chrm][left_coord//res:right_coord//res])
                else:
                    left_coord = max(region["start"] - args.upstream, 0)
                    right_coord = min(region["end"] + args.downstream, array_len * res)
                    sa.append(sample_array[this_chrm][left_coord//res:right_coord//res])
                    ra.append(reference_array[this_chrm][left_coord//res:right_coord//res])
            # unstranded locations
            else:
                left_coord = max(region["start"] - args.upstream, 0)
                right_coord = min(region["end"] + args.downstream, array_len * res)
                sa.append(sample_array[this_chrm][left_coord//res:right_coord//res])
                ra.append(reference_array[this_chrm][left_coord//res:right_coord//res])

    else:
        for contig in spike_contigs:
            sa.append(sample_array[contig])
            ra.append(reference_array[contig])
        if minus_sa is not None:
            sa.append(minus_sa[contig])
            ra.append(minus_ra[contig])

    sa = np.concatenate(sa)
    ra = np.concatenate(ra)
    # Add a pseudo count
    sa = sa + pseudocount
    ra = ra + pseudocount
    

    # adapted from EdgeR code here: https://rdrr.io/bioc/edgeR/src/R/calcNormFactors.R

    # calculate M values
    Mvals = np.log2(sa/total_sa) - np.log2(ra/total_ra)
    # Calculate A values
    Avals = (np.log2(sa/total_sa) + np.log2(ra/total_ra))/2
    # calculate estimated variance
    v = (total_sa - sa)/total_sa/sa + (total_ra - ra)/total_ra/ra

    # get finite values
    mask = np.logical_and(np.logical_and(np.isfinite(Mvals),np.isfinite(Avals)), Avals > A_cutoff)

    Mvals = Mvals[mask]
    Avals = Avals[mask]
    v = v[mask]

    if(max(abs(Mvals)) < 1e-6):
        return 1

    # Determiming trimming based on A and M
    # Following EdgeRs copy of the original mean function
    total_M = len(Mvals)
    lo_M = np.floor(total_M*M_trim) + 1
    hi_M = total_M + 1 - lo_M
    lo_A = np.floor(total_M*A_trim) + 1
    hi_A = total_M + 1 - lo_A
    M_trim_mask = np.logical_and(stats.rankdata(Mvals) >= lo_M, stats.rankdata(Mvals)<= hi_M)
    A_trim_mask = np.logical_and(stats.rankdata(Avals) >= lo_A, stats.rankdata(Avals)<= hi_A)
    Mvals = Mvals[np.logical_and(M_trim_mask, A_trim_mask)]
    Avals = Avals[np.logical_and(M_trim_mask, A_trim_mask)]
    v = v[np.logical_and(M_trim_mask, A_trim_mask)]

    # Calculate TMM values for final TMM calc
    TMM = np.power(2, np.sum(Mvals/v)/np.sum(1/v))
    if not np.isfinite(TMM):
        TMM = 1
    return TMM

def geometric_mean(arrays, axis = 0, pseudocount = 1.0):
    arrays = [np.log(array + pseudocount) for array in arrays]
    return np.exp(np.nanmean(arrays, axis=axis))

def add_pseudocount(array, pseudocount):
    out_array = {}
    for chrm in array.keys():
        out_array[chrm] = array[chrm] + pseudocount
    return out_array

def scale_array(array, scale_factor, pseudocount):
    if scale_factor == 0 or scale_factor is None:
        raise ValueError("Scale factor value was zero or undefined.")
    out_array = {}
    for chrm in array.keys():
        out_array[chrm] = (array[chrm] + pseudocount)/scale_factor
    return out_array

def scale_byfactor(array, scalefactor_table, scalefactor_id, pseudocount):
    import pandas as pd

    sf_table = pd.read_csv(scalefactor_table, sep = "\t")
    scale_factor = sf_table.loc[sf_table["sample_name"] == scalefactor_id[0], scalefactor_id[1]].values[0]
    return scale_array(array, scale_factor, pseudocount)

def to_one_array(array, minus_array = None):
    indices = {}
    indices["+"] = {}
    start = 0
    out_array = []
    for chrm in array.keys():
        indices["+"][chrm] = [start, start + array[chrm].size]
        start += array[chrm].size
        out_array.append(array[chrm])

    if minus_array is not None:
        indices["-"] = {}
        for chrm in minus_array.keys():
            indices["-"][chrm] = [start, start + minus_array[chrm].size]
            start += minus_array[chrm].size
            out_array.append(minus_array[chrm])
    out_array = np.concatenate(out_array)

    return (indices, out_array)


def from_one_array(array, indices):
    plus_out_array = {}
    minus_out_array = {}
    for key in indices["+"]:
        chrm_idx = indices["+"][key]
        plus_out_array[key] = array[chrm_idx[0]:chrm_idx[1]]
    if "-" in indices.keys():
        for key in indices["-"]:
            chrm_idx = indices["-"][key]
            minus_out_array[key] = array[chrm_idx[0]:chrm_idx[1]]

    return (plus_out_array, minus_out_array)


def downsample_array(array, total, minus_array = None, seed = 42):
    indices, one_array = to_one_array(array, minus_array)

    total_count = int(np.sum(one_array))

    count_to_remove = total_count - int(total)

    if count_to_remove <= 0:
        raise ValueError("Trying to downsample more (%s) than total count (%s)"%(total, total_count))

    overall_idx = np.arange(0, one_array.size, dtype = int)

    rng = np.random.default_rng(seed = seed)

    dwnsample_idx = rng.choice(overall_idx, size = count_to_remove, p = one_array/total_count, replace = True)
    for idx in dwnsample_idx:
        one_array[idx] = one_array[idx] - 1

    return from_one_array(one_array, indices)


def downsample(array, scalefactor_table, scalefactor_id, minus_array = None):
    import pandas as pd

    sf_table = pd.read_csv(scalefactor_table, sep = "\t")
    if scalefactor_id[1] == "nonspike_frags" or scalefactor_id[1] == "spike_frags":
        total = sf_table.loc[sf_table["sample_name"] == scalefactor_id[0], "total_frags"].values[0]
        spike = sf_table.loc[sf_table["sample_name"] == scalefactor_id[0], "spike_frags"].values[0]
        if scalefactor_id[1] == "nonspike_frags":
            scale_factor = total/(total-spike)
        else:
            scale_factor = total/spike
    else:
        scale_factor = 1
    dwn_total = int(np.min(sf_table.loc[:, scalefactor_id[1]])*scale_factor)

    return downsample_array(array, dwn_total, minus_array)

def get_regression_estimates(ext_array, input_array, spike_contigs, expected_locs= None, res = None, upstream = 0, downstream = 0):
    import bed_utils
    from sklearn.linear_model import LinearRegression
    if expected_locs is not None:
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(expected_locs)
    if spike_contigs is None:
        spike_contigs = ext_array.keys()
    mask_array = {}
    for contig in spike_contigs:
        mask_array[contig] = np.ones(len(ext_array[contig]), bool)
    if expected_locs is not None:
        for region in inbed:
            if region["chrm"] in spike_contigs: 
                array_len = len(mask_array[contig])
                if region["strand"] == "-":
                    left_coord = max(region["start"] - args.downstream, 0)
                    right_coord = min(region["end"] + args.upstream, array_len * res)
                else:
                    left_coord = max(region["start"] - args.upstream, 0)
                    right_coord = min(region["end"] + args.downstream, array_len * res)
                mask_array[region["chrm"]][left_coord//res:right_coord//res] = 0
    

    # concatenate all locations where enrichment is not expected into a single array
    X = []
    y = []
    for contig in spike_contigs:
        y.append(ext_array[contig][mask_array[contig]])
        X.append(input_array[contig][mask_array[contig]])

    y = np.concatenate(y)
    X = np.concatenate(X)

    # fit input vs. extracted
    # make sure there are no NaNs from masked areas
    remove_nansandinfs = np.logical_and(np.isfinite(X), np.isfinite(y))
    y = y[remove_nansandinfs].reshape(-1,1)
    X = X[remove_nansandinfs].reshape(-1,1)
    reg = LinearRegression(fit_intercept = False).fit(X,y)
    regress_slope = reg.coef_[0,0]

    # predict expected bound regions using input 
    if expected_locs is not None:
        new_mask_array = {contig: ~mask_array[contig] for contig in mask_array}
    else:
        new_mask_array = mask_array
    full_X = []
    full_y = []
    for contig in spike_contigs:
        full_y.append(ext_array[contig][new_mask_array[contig]])
        full_X.append(input_array[contig][new_mask_array[contig]])

    full_y = np.concatenate(full_y)
    full_X = np.concatenate(full_X)
    # make sure there are no NaNs from masked areas
    remove_nansandinfs = np.logical_and(np.isfinite(full_X), np.isfinite(full_y))
    # get residuals
    resids = full_y[remove_nansandinfs] - full_X[remove_nansandinfs]*regress_slope
    
    # get sum of positive residuals and return
    if expected_locs is not None:
        resids[resids < 0] = 0
    return np.mean(resids)


def get_stranded_regression_estimates(ext_array_plus, ext_array_minus, inp_array_plus, inp_array_minus,  spike_contigs, expected_locs, res, antisense = False):
    import bed_utils
    from sklearn.linear_model import LinearRegression
    if expected_locs is not None:
        inbed = bed_utils.BedFile()
        inbed.from_bed_file(expected_locs)
    mask_array = {"-" : {}, "+" : {}}
    for contig in spike_contigs:
        mask_array["-"][contig] = np.ones(len(ext_array_plus[contig]), bool)
        mask_array["+"][contig] = np.ones(len(ext_array_minus[contig]), bool)

    if expected_locs is not None:
        for region in inbed:
            if region["chrm"] in spike_contigs:
                array_len = len(mask_array["+"][contig])
                if (region["strand"] == "-" and not antisense) or (region["strand"] == "+" and antisense):
                    left_coord = max(region["start"] - args.downstream, 0)
                    right_coord = min(region["end"] + args.upstream, array_len * res)
                    mask_array["-"][region["chrm"]][left_coord//res:right_coord//res] = 0
                else:
                    left_coord = max(region["start"] - args.upstream, 0)
                    right_coord = min(region["end"] + args.downstream, array_len * res)
                    mask_array["+"][region["chrm"]][left_coord//res:right_coord//res] = 0

    

    # concatenate all locations where enrichment is not expected into a single array
    X = []
    y = []
    for contig in spike_contigs:
        y.append(ext_array_plus[contig][mask_array["+"][contig]])
        y.append(ext_array_minus[contig][mask_array["-"][contig]])
        X.append(inp_array_plus[contig][mask_array["+"][contig]])
        X.append(inp_array_minus[contig][mask_array["-"][contig]])

    y = np.concatenate(y)
    X = np.concatenate(X)

    # fit input vs. extracted
    # make sure there are no NaNs from masked areas
    remove_nansandinfs = np.logical_and(np.isfinite(X), np.isfinite(y))
    y = y[remove_nansandinfs].reshape(-1,1)
    X = X[remove_nansandinfs].reshape(-1,1)
    reg = LinearRegression(fit_intercept = False).fit(X,y)
    regress_slope = reg.coef_[0,0]
    
    # predict expected bound regions using input 
    full_X = []
    full_y = []

    if expected_locs is not None:
        new_mask_array = {"-" : {}, "+" : {}}
        for strand in ["+", "-"]:
            new_mask_array[strand] = {contig: ~mask_array[strand][contig] for contig in mask_array[strand]}
    else:
        new_mask_array = mask_array


    for contig in spike_contigs:
        full_y.append(ext_array_plus[contig][new_mask_array["+"][contig]])
        full_y.append(ext_array_minus[contig][new_mask_array["-"][contig]])
        full_X.append(inp_array_plus[contig][new_mask_array["+"][contig]])
        full_X.append(inp_array_minus[contig][new_mask_array["-"][contig]])

    full_y = np.concatenate(full_y)
    full_X = np.concatenate(full_X)

    # make sure there are no NaNs from masked areas
    remove_nansandinfs = np.logical_and(np.isfinite(full_X), np.isfinite(full_y))
    # get residuals
    resids = full_y[remove_nansandinfs] - full_X[remove_nansandinfs]*regress_slope
    

    # get sum of positive residuals and return
    resids[resids < 0] = 0
    return np.mean(resids)

def normfactor_main(args):
    import pandas as pd 

    frag_table = pd.read_csv(args.fragCountTable, sep = "\t")
    md = pd.read_csv(args.metaDataTable)

    # get the total number of frags per each sample
    overall = frag_table.assign(total_frag = \
            lambda x: x.groupby(["sample_name"])['fragments']\
            .transform(lambda x: x.sum()))\
    .drop(["contig","fragments"], axis=1)\
    .drop_duplicates()
    total_frag_sf_column = pd.DataFrame(data = {'total_frag_sfs': overall["total_frag"]/np.mean(overall["total_frag"]), 'sample_name': overall["sample_name"]})
    overall = overall.merge(total_frag_sf_column, on = 'sample_name', how = 'left')


    if args.spikecontigs is not None:
        # spike fragments
        spike_frags = frag_table[frag_table["contig"].isin(args.spikecontigs)].assign(spike_frag = \
                lambda x: x.groupby(["sample_name"])['fragments']\
                .transform(lambda x: x.sum()))\
        .drop(["contig","fragments"], axis=1)\
        .drop_duplicates()
    
        overall = overall.merge(spike_frags, on='sample_name', how = 'left')
        overall = overall.assign(nonspike_frag = overall['total_frag'] - overall['spike_frag'])

        spike_frag_sf_column = pd.DataFrame(data = {'spike_frag_sfs': overall["spike_frag"]/np.mean(overall["spike_frag"]), 'sample_name': overall["sample_name"]})
        nonspike_frag_sf_column = pd.DataFrame(data = {'nonspike_frag_sfs': overall["nonspike_frag"]/np.mean(overall["nonspike_frag"]), 'sample_name': overall["sample_name"]})

        overall = overall.merge(spike_frag_sf_column, on = 'sample_name', how = 'left')
        overall = overall.merge(nonspike_frag_sf_column, on = 'sample_name', how = 'left')

        # get the change in the fraction of spike-in fragments from ext to inp
        diff_in_spike_values = []
        for sample in args.samples:
            #print(overall.loc[overall["sample_name"] == sample, "spike_frag"], sample)
            diff_in_spike_values.append(np.power(2, np.log2(\
                    overall.loc[overall["sample_name"] == sample, "spike_frag"].values[0]/ \
                    overall.loc[overall["sample_name"] == sample, "total_frag"].values[0])- \
                    np.log2( \
                    overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0], "spike_frag"].values[0]/\
                    overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0], "total_frag"].values[0])))

        diff_in_spike = pd.DataFrame(data = {'diff_spike_sfs': diff_in_spike_values, 'sample_name': args.samples})
        overall = overall.merge(diff_in_spike, on = 'sample_name', how = 'left')
        diff_in_spike_rpm = pd.DataFrame(data = {'diff_spike_rpm_sfs': overall["diff_spike_sfs"]*overall["total_frag_sfs"],\
                'sample_name': overall['sample_name']})

        overall = overall.merge(diff_in_spike_rpm, on = 'sample_name', how = 'left')

    # DEseq2 size factors
    # stranded version
    if args.ext_bws_minus is not None:
        bw_plus_handles = open_multiple_bigwigs(args.ext_bws)
        bw_plus_arrays = convert_bigwigs_to_arrays(bw_plus_handles, res = args.res)

        bw_minus_handles = open_multiple_bigwigs(args.ext_bws_minus)
        bw_minus_arrays = convert_bigwigs_to_arrays(bw_minus_handles, res = args.res)

        geom_means_plus = multicompare_within_group(bw_plus_arrays, lambda x, axis: geometric_mean(x, axis, pseudocount = args.pseudocount))
        geom_means_minus = multicompare_within_group(bw_minus_arrays, lambda x, axis: geometric_mean(x, axis, pseudocount = args.pseudocount))
        deseq2_sfs = []
        for array_plus, array_minus, sample in zip(bw_plus_arrays, bw_minus_arrays, args.samples):
            compared_minus = compare_divide(add_pseudocount(array_minus, args.pseudocount), geom_means_minus)
            compared_plus = compare_divide(add_pseudocount(array_plus, args.pseudocount), geom_means_plus)
            deseq2_sfs.append(stranded_single_num_summary(compared_plus, compared_minus, np.nanmedian))

        deseq2_column = pd.DataFrame(data = {'deseq2_sfs': deseq2_sfs, 'sample_name': args.samples})
        overall = overall.merge(deseq2_column, on = 'sample_name', how = 'left')
    # unstranded version
    else:
        bw_handles = open_multiple_bigwigs(args.ext_bws)
        bw_arrays = convert_bigwigs_to_arrays(bw_handles, res = args.res)
        geom_means = multicompare_within_group(bw_arrays, lambda x, axis: geometric_mean(x, axis, pseudocount = args.pseudocount))
        deseq2_sfs = []
        for array, sample in zip(bw_arrays, args.samples):
            deseq2_sfs.append(single_num_summary(compare_divide(add_pseudocount(array, args.pseudocount), geom_means), np.nanmedian))
        deseq2_column = pd.DataFrame(data = {'deseq2_sfs': deseq2_sfs, 'sample_name': args.samples})
        overall = overall.merge(deseq2_column, on = 'sample_name', how = 'left')

    # TMMs
    # find sample with highest frag count
    all_frags = []
    for sample_name in args.samples:
        all_frags.append(overall.loc[overall["sample_name"] == sample_name, "total_frag"].values[0])
    high_frag_index = np.argmax(all_frags)
    ref_sample = args.samples[high_frag_index]
    # stranded version
    if args.ext_bws_minus is not None:

        tmm_sfs = []
        lib_size = []

        for array_plus, array_minus, sample in zip(bw_plus_arrays, bw_minus_arrays, args.samples):
            lib_size.append(overall.loc[overall["sample_name"] == sample, "total_frag"].values[0])
            if sample == ref_sample:
                tmm_sfs.append(1)
            else:
                tmm_sfs.append(get_TMM_factor(array_plus, bw_plus_arrays[high_frag_index], \
                        overall.loc[overall["sample_name"] == sample, "total_frag"].values[0],
                        overall.loc[overall["sample_name"] == ref_sample, "total_frag"].values[0],
                        array_plus.keys(),
                        minus_sa = array_minus,
                        minus_ra = bw_minus_arrays[high_frag_index],
                        pseudocount = args.pseudocount))

    # unstranded version
    else:

        tmm_sfs = []
        lib_size = []

        for array, sample in zip(bw_arrays, args.samples):
            lib_size.append(overall.loc[overall["sample_name"] == sample, "total_frag"].values[0])
            if sample == ref_sample:
                tmm_sfs.append(1)
            else:
                tmm_sfs.append(get_TMM_factor(array, bw_arrays[high_frag_index], \
                        overall.loc[overall["sample_name"] == sample, "total_frag"].values[0],
                        overall.loc[overall["sample_name"] == ref_sample, "total_frag"].values[0],
                        array.keys(),
                        pseudocount = args.pseudocount))



    tmm_column = pd.DataFrame(data = {'tmm_sfs': tmm_sfs, 'sample_name': args.samples})
    lib_scaled = np.multiply(tmm_sfs, lib_size)
    tmm2_column = pd.DataFrame(data = {'tmm_rpm_sfs': lib_scaled/(np.mean(lib_scaled)), 'sample_name': args.samples})
    overall = overall.merge(tmm_column, on = 'sample_name', how = 'left')
    overall = overall.merge(tmm2_column, on = 'sample_name', how = 'left')

    # DEseq2 size factors spike-in only
    if args.spikecontigs is not None or args.expected_regions is not None:
    
        # stranded version
        if args.ext_bws_minus is not None:

            deseq2_sfs = []
            for array_plus, array_minus, sample in zip(bw_plus_arrays, bw_minus_arrays, args.samples):
                compared_minus = compare_divide(add_pseudocount(array_minus, args.pseudocount), geom_means_minus)
                compared_plus = compare_divide(add_pseudocount(array_plus, args.pseudocount), geom_means_plus)
                deseq2_sfs.append(stranded_single_num_summary(compared_plus, compared_minus, np.nanmedian, contigs = args.spikecontigs,
                    expected_locs = args.expected_regions, res = args.res, antisense = args.antisense))

            deseq2_column = pd.DataFrame(data = {'deseq2_spike_sfs': deseq2_sfs, 'sample_name': args.samples})
            overall = overall.merge(deseq2_column, on = 'sample_name', how = 'left')
        # unstranded version
        else:
            deseq2_sfs = []
            for array, sample in zip(bw_arrays, args.samples):
                deseq2_sfs.append(single_num_summary(compare_divide(add_pseudocount(array, args.pseudocount), geom_means), np.nanmedian, 
                    contigs = args.spikecontigs, expected_locs = args.expected_regions, res = args.res))
            deseq2_column = pd.DataFrame(data = {'deseq2_spike_sfs': deseq2_sfs, 'sample_name': args.samples})
            overall = overall.merge(deseq2_column, on = 'sample_name', how = 'left')


        # TMMs
        # stranded version
        if args.ext_bws_minus is not None:
    
            tmm_sfs = []
            lib_size = []
    
            for array_plus, array_minus, sample in zip(bw_plus_arrays, bw_minus_arrays, args.samples):
                lib_size.append(overall.loc[overall["sample_name"] == sample, "total_frag"].values[0])
                if sample == ref_sample:
                    tmm_sfs.append(1)
                else:
                    tmm_sfs.append(get_TMM_factor(array_plus, bw_plus_arrays[high_frag_index], \
                            overall.loc[overall["sample_name"] == sample, "total_frag"].values[0],
                            overall.loc[overall["sample_name"] == ref_sample, "total_frag"].values[0],
                            args.spikecontigs,
                            minus_sa = array_minus,
                            minus_ra = bw_minus_arrays[high_frag_index],
                            expected_locs = args.expected_regions,
                            res = args.res,
                            antisense = args.antisense,
                            pseudocount = args.pseudocount))
    
        # unstranded version
        else:
    
            tmm_sfs = []
            lib_size = []
    
            for array, sample in zip(bw_arrays, args.samples):
                lib_size.append(overall.loc[overall["sample_name"] == sample, "total_frag"].values[0])
                if sample == ref_sample:
                    tmm_sfs.append(1)
                else:
                    tmm_sfs.append(get_TMM_factor(array, bw_arrays[high_frag_index], \
                            overall.loc[overall["sample_name"] == sample, "total_frag"].values[0],
                            overall.loc[overall["sample_name"] == ref_sample, "total_frag"].values[0],
                            args.spikecontigs,
                            expected_locs = args.expected_regions,
                            res = args.res,
                            pseudocount = args.pseudocount))
    

        tmm_column = pd.DataFrame(data = {'tmm_spike_sfs': tmm_sfs, 'sample_name': args.samples})
        lib_scaled = np.multiply(tmm_sfs, lib_size)
        tmm2_column = pd.DataFrame(data = {'tmm_spike_rpm_sfs': lib_scaled/(np.mean(lib_scaled)), 'sample_name': args.samples})
        overall = overall.merge(tmm_column, on = 'sample_name', how = 'left')
        overall = overall.merge(tmm2_column, on = 'sample_name', how = 'left')


    # enrichment based
    if args.spikecontigs is not None:
        # stranded version
        if args.ext_bws_minus is not None:
            inp_plus_bw_handles = open_multiple_bigwigs(args.inp_bws)
            inp_plus_bw_arrays = convert_bigwigs_to_arrays(inp_plus_bw_handles, res = args.res)

            inp_minus_bw_handles = open_multiple_bigwigs(args.inp_bws_minus)
            inp_minus_bw_arrays = convert_bigwigs_to_arrays(inp_minus_bw_handles, res = args.res)
            regress_resids = []
            regress_rpm = []

            regress_total = []
            regress_total_rpm = []

            for ext_array_plus, ext_array_minus, inp_array_plus, inp_array_minus, sample in zip(bw_plus_arrays, bw_minus_arrays, inp_plus_bw_arrays, inp_minus_bw_arrays, args.samples):
                regress_resids.append(get_stranded_regression_estimates(\
                        scale_array(ext_array_plus, overall.loc[overall["sample_name"] == sample, "spike_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(ext_array_minus, overall.loc[overall["sample_name"] == sample, "spike_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array_plus, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"spike_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array_minus, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"spike_frag"].values[0]/1e6, args.pseudocount),\
                        args.spikecontigs,
                        args.expected_regions,
                        args.res,
                        antisense = args.antisense))
                regress_rpm.append(overall.loc[overall["sample_name"] == sample, "nonspike_frag_sfs"].values[0])


                regress_total.append(get_stranded_regression_estimates(\
                        scale_array(ext_array_plus, overall.loc[overall["sample_name"] == sample, "total_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(ext_array_minus, overall.loc[overall["sample_name"] == sample, "total_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array_plus, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"total_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array_minus, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"total_frag"].values[0]/1e6, args.pseudocount),\
                        args.spikecontigs,
                        args.expected_regions,
                        args.res,
                        antisense = args.antisense))
                regress_rpm.append(overall.loc[overall["sample_name"] == sample, "total_frag_sfs"].values[0])

            regress_column = pd.DataFrame(data = {'regress_resids': regress_resids, 'sample_name':args.samples}) 
            regress_sfs = pd.DataFrame(data = {'regress_sfs': regress_resids/np.mean(regress_resids), 'sample_name': args.samples})   
            lib_scaled = (regress_resids/np.mean(regress_resids))*regress_rpm
            regress_rpm_sfs = pd.DataFrame(data = {'regress_rpm_sfs': lib_scaled, 'sample_name': args.samples})

            overall = overall.merge(regress_column, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_sfs, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_rpm_sfs, on = 'sample_name', how = 'left')


            regress_total_column = pd.DataFrame(data = {'regress_total_resids': regress_total, 'sample_name':args.samples}) 
            regress_total_sfs = pd.DataFrame(data = {'regress_total_sfs': regress_total/np.mean(regress_total), 'sample_name': args.samples})   
            lib_total_scaled = (regress_total/np.mean(regress_total))*regress_total_rpm
            regress_total_rpm_sfs = pd.DataFrame(data = {'regress_total_rpm_sfs': lib_total_scaled, 'sample_name': args.samples})

            overall = overall.merge(regress_total_column, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_total_sfs, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_total_rpm_sfs, on = 'sample_name', how = 'left')
    
        # unstranded version
        else:
            inp_bw_handles = open_multiple_bigwigs(args.inp_bws)
            inp_bw_arrays = convert_bigwigs_to_arrays(inp_bw_handles, res = args.res)

            regress_resids = []
            regress_rpm = []
            regress_spike_rpm = []

            regress_total = []
            regress_total_rpm = []

            for ext_array, inp_array, sample in zip(bw_arrays, inp_bw_arrays, args.samples):
                regress_resids.append(get_regression_estimates(scale_array(ext_array, overall.loc[overall["sample_name"] == sample, "spike_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"spike_frag"].values[0]/1e6, args.pseudocount),\
                        args.spikecontigs,
                        args.expected_regions,
                        args.res,
                        args.upstream,
                        args.downstream))
                regress_rpm.append(overall.loc[overall["sample_name"] == sample, "nonspike_frag_sfs"].values[0])
                regress_spike_rpm.append(overall.loc[overall["sample_name"] == sample, "spike_frag_sfs"].values[0])

                regress_total.append(get_regression_estimates(scale_array(ext_array, overall.loc[overall["sample_name"] == sample, "total_frag"].values[0]/1e6, args.pseudocount),\
                        scale_array(inp_array, overall.loc[overall["sample_name"] == md.loc[md["sample_name"] == sample, "input_sample"].values[0],"total_frag"].values[0]/1e6, args.pseudocount),\
                        args.spikecontigs,
                        args.expected_regions,
                        args.res,
                        args.upstream,
                        args.downstream))
                regress_total_rpm.append(overall.loc[overall["sample_name"] == sample, "total_frag_sfs"].values[0])

            regress_column = pd.DataFrame(data = {'regress_resids': regress_resids, 'sample_name':args.samples})
            regress_sfs = pd.DataFrame(data = {'regress_sfs': regress_resids/np.mean(regress_resids), 'sample_name': args.samples})    
            lib_scaled = (regress_resids/np.mean(regress_resids))*regress_rpm
            regress_rpm_sfs = pd.DataFrame(data = {'regress_rpm_sfs': lib_scaled, 'sample_name': args.samples})
            #spike_lib_scaled = (regress_resids/np.mean(regress_resids))*regress_spike_rpm
            #regress_spike_rpm_sfs = pd.DataFrame(data = {'regress_spike_rpm_sfs': spike_lib_scaled/np.mean(spike_lib_scaled), 'sample_name': args.samples})

            regress_total_column = pd.DataFrame(data = {'regress_total_resids': regress_total, 'sample_name':args.samples})
            regress_total_sfs = pd.DataFrame(data = {'regress_total_sfs': regress_total/np.mean(regress_total), 'sample_name':args.samples})
            total_lib_scaled = (regress_total/np.mean(regress_total))*regress_total_rpm
            regress_total_rpm_sfs = pd.DataFrame(data = {'regress_total_rpm_sfs': total_lib_scaled, 'sample_name': args.samples})
            
            overall = overall.merge(regress_column, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_sfs, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_rpm_sfs, on = 'sample_name', how = 'left')
            #overall = overall.merge(regress_spike_rpm_sfs, on = 'sample_name', how = 'left')

            overall = overall.merge(regress_total_column, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_total_sfs, on = 'sample_name', how = 'left')
            overall = overall.merge(regress_total_rpm_sfs, on = 'sample_name', how = 'left')

    overall.to_csv(args.outfile, index = False, sep = "\t", float_format='%.4f')
    
    if args.ext_bws_minus is not None: 
        close_multiple_bigwigs(bw_plus_handles)
        close_multiple_bigwigs(bw_minus_handles)
        close_multiple_bigwigs(inp_minus_bw_handles)
        close_multiple_bigwigs(inp_plus_bw_handles)
    else:
        close_multiple_bigwigs(bw_handles)
        close_multiple_bigwigs(inp_bw_handles)

def convert_main(args):
    summary_operation_dict = {"mean": np.nanmean,
            "sum": np.nansum,
            "median": np.nanmedian,
            "max": np.nanmax,
            "min": np.nanmin}
    # read in files
    inf = pyBigWig.open(args.infile)

    arrays = bigwig_to_arrays(inf, res = args.inres)
    # convert resolution
    if args.inres != args.outres:
        arrays = change_array_resolution(arrays, inf.chroms(), summary_operation_dict[args.summary_func], \
                                     res_from = args.inres, res_to = args.outres)
    # filter chrms you want
    if args.chrms is not None:
        out_arrays = {}
        for chrm in args.chrms:
            out_arrays[chrm] = arrays[chrm]
        arrays = out_arrays
    #write file out
    if args.outfmt == ".bedgraph":
        write_arrays_to_bedgraph(args.outfile, arrays, inf.chroms(), res = args.outres, dropNaNsandInfs = args.dropNaNsandInfs, header = args.header)
    elif args.outfmt == ".wig":
        write_arrays_to_wig(args.outfile, arrays, inf.chroms(), res = args.outres, dropNaNsandInfs = args.dropNaNsandInfs, header = args.header)
    elif args.outfmt == ".bw":
        write_arrays_to_bigwig(args.outfile, arrays, inf.chroms(), res = args.outres, dropNaNsandInfs = args.dropNaNsandInfs)
    else:
        raise ValueError("Only support --outfmt of '.bedgraph', '.wig', and '.bw'. Not %s"%(args.outfmt))
    inf.close()
 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Tools to work with BigWig files")
    subparsers = parser.add_subparsers(help = "Supported Commands")

    # manipulate verb
    parser_manipulate = subparsers.add_parser("manipulate", help = "manipulate BigWig files")
    parser_manipulate.add_argument('infile', type=str, 
            help="file to convert from")
    parser_manipulate.add_argument('outfile', type=str, help="file to convert to")
    parser_manipulate.add_argument("--minus_strand", type = str, help = "minus strand file for stranded data")
    parser_manipulate.add_argument("--minus_strand_out", type = str, help = "minus strand out file for stranded data")
    parser_manipulate.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_manipulate.add_argument('--operation', type=str, default=None,
            help="operation to perform before writing out file. \
            All operations, neccesitate conversion to array internally \
            options {'RobustZ', 'Median_norm', 'background_subtract', 'scale_max', \
            'query_subtract', 'query_scale', 'gauss_smooth', 'flat_smooth', 'savgol_smooth', 'Bspline_smooth',\
            'GLMGam_smooth','scale_byfactor', 'downsample'}")
    parser_manipulate.add_argument('--fixed_regions', type=str, default=None,
            help="bed file containing a fixed set of regions on which to average over.")
    parser_manipulate.add_argument('--query_regions', type=str, default=None,
            help="bed file containing regions to consider.")
    parser_manipulate.add_argument('--number_of_regions', type = int, default = 20,
            help = "number of regions to include in background or max")
    parser_manipulate.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_manipulate.add_argument('--pseudocount', default = 0, type = float,
            help = "Add value to all unmasked regions before normalization. Only\
                    applicable to the median and scale_byfactor normalization. Default = 0 ")
    parser_manipulate.add_argument('--summary_func', default = "mean", type = str,
            help = "For methods that use regions. What summary function do you want to use over \
                    said regions. Options: mean, median, sum \
                    Default = mean ")
    parser_manipulate.add_argument('--wsize', type = int, help = "For smoothing methods \
            specifies half the window size in units of res. Thus, a wsize of 5 with \
            resolution 5 is a 5*5 = 25 bp window. wsize is half the window so \
            the total window size is wsize*2 + 1")
    parser_manipulate.add_argument('--edge', default = "mirror", type = str,
            help = "For smoothing methods how to deal with the edge? 'mirror' takes \
            half the window on each end and mirrors it. 'wrap' wraps the entire array \
            around which is ideal for circular chromosomes. Default = 'mirror'"),
    parser_manipulate.add_argument('--savgol_poly', type = int,
            help = "For Savtizky-Golay smoothing what order polynomial? Must not \
            be larger than the full window size.")
    parser_manipulate.add_argument('--gauss_sigma', type = int,
            help = "For gaussian smoothing what is sigma? Must not \
            be larger than the full window size. Default is wsize*2/6")
    parser_manipulate.add_argument('--smooth_scaler', type = str,
            help = "Scale the smooth signal? min_max = min max scaler. fractional = array/sum(array)*len(array)",
            default = None)
    parser_manipulate.add_argument('--smooth_knots', type = str,
            help = "Text file with one row of chrm\\tknot per knot; Otherwise let smoother find optimal knots", default = None)
    parser_manipulate.add_argument('--gam_max_knots', type = int,
            help = "Maximum total knots for gam smoother. Default = 10", default = 10)
    parser_manipulate.add_argument('--scalefactor_table', type = str, default = None,
            help = "A table of scale factors for scale_byfactor")
    parser_manipulate.add_argument('--scalefactor_id', type = str, nargs = 2, default = None,
            help = "The sample and scale factor to choose from the table")
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
    parser_query.add_argument("--minus_strand", nargs = "+", type = str, help = "minus strand files for stranded data")
    parser_query.add_argument("--antisense", action = "store_true", help = "compute antisense instead of sense values for each region")
    parser_query.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_query.add_argument('--res_to', type=int, default=None,
            help="Resolution to output at. Default same resolution as input. Note this \
            should be set no lower than the resolution of the input file. \
            Resolution changes are made by averaging over window before any calculations are done")
    parser_query.add_argument('--regions', type=str, help = "regions to grab data from")
    parser_query.add_argument('--coords', type = str, help = "When running in 'identity' mode do you want 'relative_start', 'relative_center', 'relative_end' or 'absolute' coords? Default = 'absolute'")
    parser_query.add_argument('--upstream', type=int, help = "bp upstream to add", default = 0)
    parser_query.add_argument('--downstream', type = int, help = "bp downstream to add", default = 0)
    parser_query.add_argument('--samp_names', type = str, nargs = "+", help = "sample names for each file")
    parser_query.add_argument('--summarize', type = str, help = "How to summarize data. 'identity' reports \
            each data point. 'single' gives a single number summary. Default = 'identity'", default = 'identity')
    parser_query.add_argument('--frac_na', type = float, help = "When reporting summaries for regions, how much of the region\
            can be NA before reporting the value as NA? default = 0.25", default = 0.25)
    parser_query.add_argument('--summary_func', type = str, help = "What function to use to summarize data when not using \
            'identity' summary. mean, median, max, min supported. Additionally, traveling ratio ('TR') and relative polymerase progression ('RPP'), \
            traveling ratio A window ('TR_A'), traveling ratio B window ('TR_B'),  Gini coefficient ('Gini') and 'summit_loc' (local peak identification) are supported. Default = 'mean'", default = "mean")
    parser_query.add_argument('--gzip', action = "store_true", help = "gzips the output if flag is included")
    parser_query.add_argument('--TR_A_center', type = int, help = "center of window A in fixed traveling ratio. In relative bp to region start")
    parser_query.add_argument('--wsize', type = int, default = 50, help = "Size of half window in bp for calcs that use windows. Default = 50 bp")
    parser_query.set_defaults(func=query_main)

    # compare verb
    parser_compare = subparsers.add_parser("compare", help = "Calculate operations between two BigWig files")
    parser_compare.add_argument("infile1", type=str, help = "input file 1")
    parser_compare.add_argument("infile2", type=str, help = "input file 2")
    parser_compare.add_argument("outfile", type=str, help = "file to output to")
    parser_compare.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input files")
    parser_compare.add_argument('--pseudocount', default = 0, type = int,
            help = "Add value to all unmasked regions before ratio calculations. \
                    Default = 0 ")
    parser_compare.add_argument('--operation', type=str, default="log2ratio",
            help="Default is log2ratio i.e. log2(infile1) - log2(infile2). Other \
                    options include: add, subtract, divide, recipratio (invert ratios less than 1. Center at 0)")
    parser_compare.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_compare.set_defaults(func=compare_main)

    # multicompare verb
    parser_multicompare = subparsers.add_parser("multicompare", help = "Calculate operations between two groups of BigWig files or a summary operation on one group")
    parser_multicompare.add_argument("outfile", type=str, help = "file to output to")
    parser_multicompare.add_argument("--groupA", type=str, nargs ="+", help = "input file 1")
    parser_multicompare.add_argument("--groupB", type=str, nargs = "+", help = "input file 2")
    parser_multicompare.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input files")
    parser_multicompare.add_argument('--within_operation', type = str, default = "mean",
            help = "Default is mean. Options are median, max, and min.")
    parser_multicompare.add_argument('--between_operation', type=str, default="log2ratio",
            help="Default is log2ratio i.e. log2(infile1) - log2(infile2). Other \
                    options include: add, subtract, divide, recipratio (invert ratios less than 1. Center at 0)")
    parser_multicompare.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_multicompare.set_defaults(func=multicompare_main)

    # normfactor verb
    parser_normfactor = subparsers.add_parser("normfactor", help = "Calculate normalization factors based on spike-ins or internal controls")
    parser_normfactor.add_argument("outfile", type=str, help = "file to output to")
    parser_normfactor.add_argument("fragCountTable", type = str, help = "table of fragment counts per contig for each sample")
    parser_normfactor.add_argument("metaDataTable", type = str, help = "table of sample relationships")
    parser_normfactor.add_argument("--ext_bws", type=str, nargs ="+", help = "extracted raw counts per region")
    parser_normfactor.add_argument("--ext_bws_minus", type=str, nargs ="+", help = "extracted raw counts per region for minus strand")
    parser_normfactor.add_argument("--inp_bws", type=str, nargs ="+", help = "input raw counts per region")
    parser_normfactor.add_argument("--inp_bws_minus", type=str, nargs ="+", help = "input raw counts per region for minus strand")
    parser_normfactor.add_argument("--samples", type=str, nargs ="+", help = "sample names to determine size factors for. Should match extracted bws")
    parser_normfactor.add_argument("--spikecontigs", type=str, nargs ="+", help = "sample names to determine size factors for. Should match extracted bws")
    parser_normfactor.add_argument('--pseudocount', default = 1, type = float,
            help = "Add value to all unmasked regions before ratio calculations. \
                    Default = 1 ")
    parser_normfactor.add_argument('--expected_regions', type = str, help = "bed file of regions where signal is expected")
    parser_normfactor.add_argument('--upstream', type=int, help = "bp upstream to add to expected regions", default = 0)
    parser_normfactor.add_argument('--downstream', type = int, help = "bp downstream to add to expected regions", default = 0)
    parser_normfactor.add_argument('--res', type=int, default=1,
            help="Resolution to compute statistics at. Default 1bp. Note this \
            should be set no lower than the resolution of the input file")
    parser_normfactor.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output bigwig")
    parser_normfactor.add_argument("--antisense", action = "store_true", help = "Use antisense instead of sense values for each expected region")
    parser_normfactor.set_defaults(func=normfactor_main)


    # convert verb
    parser_convert = subparsers.add_parser("convert", help = "Convert a bigwig file into a bedGraph file")
    parser_convert.add_argument('infile', type = str,
            help = "file to convert from")
    parser_convert.add_argument('outfile', type = str, help = "file to convert to")
    parser_convert.add_argument('--inres', type = int, default = 1, help = "input resolution. Default 1")
    parser_convert.add_argument('--outres', type = int, default = 1000, help = "output resolution. Default 1000")
    parser_convert.add_argument('--summary_func', type = str, default = "mean", help = "function to collapse values to lower resolutions, default = mean")
    parser_convert.add_argument('--header', action = "store_true", help = "include a header?")
    parser_convert.add_argument('--outfmt', type = str, default = ".bedgraph", help = "what output format? default = '.bedgraph'. Use '.bw' for outputting a bigwig")
    parser_convert.add_argument('--dropNaNsandInfs', action="store_true",
            help = "Drop NaNs and Infs from output file")
    parser_convert.add_argument('--chrms', type = str, nargs = "+", default = ".bedgraph", help = "what output format? default = '.bedgraph'. Use '.bw' for outputting a bigwig")
    parser_convert.set_defaults(func=convert_main)

    
    args = parser.parse_args()
    args.func(args) 
