import numpy as np
from scipy import stats
import logging

def call_bool_peaks(array, consolidate = 0):
    """ Take a logical array and find the start and stop of ones across
    the array. Consolidate any peaks that are within consolidate.

    TODO: Deal with peaks that go over the end of the array

    Args:
        array (np.array) - logical array to call peaks from.
        consolidate (int) - number of bins to consolidate peaks over
    Returns:
        peak_indices (list of lists) - a list of [start, stop] indices in the
                                        array
    """
    # first find all the places where there is a change from 1 to 0 or vice
    # versa, here we pad the incoming array with a zero on each end, then take
    # the difference along the array and finally take the absolute value to flip
    # the places where the difference was 0 - 1
    changes = np.abs(np.diff(np.concatenate(([0], array.view(np.int8), [0]))))
    # changes is now a logical array, I then find all the indices where changes
    # happen and reshape them into an ndarray of start, stop locations
    start_stops = np.where(changes)[0].reshape(-1, 2)
    if start_stops.size == 0:
        logging.warning("No bp was considered to be within a peak.")
        return []

    # now lets consolidate any peaks that are within consolidate
    consolidate_peaks = [[start_stops[0][0], start_stops[0][1]]]
    consolidated_peaks = 0
    removed_peaks = 0
    for start, stop in start_stops[1:]:
        if start - consolidate_peaks[-1][1] < consolidate:
            consolidate_peaks[-1][1] = stop
            consolidated_peaks += 1
        else:
            if stop-start > consolidate:
                consolidate_peaks.append([start, stop])
            else:
                removed_peaks += 1
    logging.info("Consolidated %s peaks within %s bps"%(consolidated_peaks, consolidate))
    logging.info("Removed %s peaks < %s bps"%(removed_peaks, consolidate))
    return consolidate_peaks

def threshold_peaks(array, threshold, consolidate=1, below = False):
    """ Function to call peaks by a threshold and consolidate peaks that are
    nearby

    Args:
        array - 1 dimensional numpy array
        threshold - float, a hard threshold for a peak 
        consolidate - int, number of adjacent locations to consolidate
        below - bool, consider peaks below the threshold (above if false)
    Returns:
        list - start,end tuples of peaks
    """
    if below:
        peaks = array < threshold
    else:
        peaks = array > threshold
    return call_bool_peaks(peaks, consolidate = consolidate)

def normalize_1D(array, center, scale):
    """
    Function to normalize a 1 dimensional array using the following formula:

    (array - center)/ scale

    Args:
        array - 1 dimensional numpy array
        center - a single value for the center (commonly mean(array))
        scale - a single value for the scale (commonly std(array))
    Returns:
        outarray - numpy array normalized
    """

    return (array - center)/scale

def robustz_1D(array, ignore_nan = True):
    """
    Function to Robust Z normalize an array. I.e.

    RZ = array - median(array)/(1.4826* MAD)

    Args:
        array - 1 dimensional numpy array
        ignore_nan - boolean, ignore individual nans and propogate those nans to final array
    Returns:
        outarray - numpy array robustZ normalized
    """
    if ignore_nan:
        # default scale is specified here
        MAD = stats.median_abs_deviation(array, nan_policy='omit', scale = 'normal')
        median = np.nanmedian(array)
    else:
        MAD = stats.median_absolute_deviation(array)
        median = np.median(array)
    return normalize_1D(array, median, MAD)

def smooth_1D(array, wsize, kernel_type = "flat", edge = "mirror", sigma = None):
    """
    Function to smooth a 1D signal using a convolution

    Args:
        array - 1 dimensional numpy array
        wsize - size of half the window
        kernel - one of flat, gaussian
        edge - one of mirror, wrap
        sigma - width of the gaussian in standard deviations. Default is (wsize*2)/6
    Returns:
        outarray - smoothed numpy array of same size
    """
    tot_size = wsize*2 + 1
    if tot_size >= len(array):
        old_wsize = wsize
        wsize = (len(array) - 1)//2
        tot_size = wsize*2 + 1
        logging.warning("Window is larger than array. Truncating window to size of array from: %s to: %s"%(old_wsize, wsize))

    if kernel_type == "flat":
        kernel = np.ones(tot_size)
    elif kernel_type == "gaussian":
        x = np.arange(-wsize, wsize + 1)
        # make the weights span the kernel 
        if sigma is None:
            sigma = (wsize*2)/6
        kernel = np.exp(- (x**2)/(2. *sigma**2))
    else:
        raise ValueError("kernel_type must be one of flat or gaussian. Not %s"%kernel_type)

    if edge == "mirror":
        # pad with a mirror reflection of the ends
        s = np.r_[array[wsize:0:-1], array, array[-2:-wsize-2:-1]]
    elif edge == "wrap":
        # pad with a wrap around the horn
        s = np.r_[array[(-wsize):], array, array[:wsize]]
    else:
        raise ValueError("edge must be one of mirror or wrap. Not %s"%edge)
    kernel = kernel/np.sum(kernel)
    out = np.convolve(s, kernel, mode = 'valid')
    return out

def savgol_1D(array, wsize, polyorder=5, deriv=0, delta = 1.0, edge = "mirror"):
    """
    Function to smooth a 1D signal using a Savitzky-Golay filter

    Args:
        array - 1 dimensional numpy array
        wsize - size of half the window
        edge - one of mirror, wrap, nearest, constant
    Returns:
        outarray - smoothed numpy array of same size

    See also:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
    """
    from scipy import signal
    tot_size = wsize*2 + 1

    if tot_size >= len(array):
        raise ValueError("Window must be smaller than array. Window %s, array %s"%(tot_size, len(array)))
    if polyorder >= tot_size:
        raise ValueError("polyorder shouldn't be larger than window. Window %s, polyorder %s"%(tot_size, polyorder))
    return signal.savgol_filter(array, tot_size, polyorder, deriv, delta, mode = edge)

def Bspline_1D(array, knot_locs = None):
    """
    Function to fit a cubic B spline with set knot locations

    Args:
        array - 1 dimensional numpy array
        knot_locs - sorted 1D knot locations. Internal nots only
    Returns:
        outarray - spline smoothed output

    """
    from scipy import interpolate
    x = np.arange(0, len(array))
    this_filter = np.isfinite(array)
    t,c,k  = interpolate.splrep(x[this_filter], array[this_filter], t=knot_locs, per =1)
    out_spline = interpolate.BSpline(t,c,k, extrapolate = False)
    return out_spline(x)


def GLMGam_1D(array, starting_knots = 10):
    """
    Function to fit GAM to smooth the data. Uses Bsplines as smoother
    in the model

    Args:
        array - 1 dimensional numpy array
        starting_knots - number of knots to start. 10 is mgcv::gam default in R
    Returns:
        outarray - spline smoothed output

    """
    import statsmodels.api as sm
    import pandas as pd
    from statsmodels.gam.api import GLMGam, CyclicCubicSplines
    # x is coordinates
    x = np.arange(0, len(array))
    # filter out nans
    this_filter = np.isfinite(array)
    # make dataframe
    d_1d = pd.DataFrame({"x":x[this_filter], "y" : array[this_filter]})
    # make spline basis. k=10 is default for mgcv::gam in R
    bs = CyclicCubicSplines(d_1d[['x']], df = [starting_knots])
    # set up formula for fitting; do poisson to avoid values < 0 since should be
    # acting on raw count data
    gam_bs = GLMGam.from_formula('y ~ 1', data = d_1d, smoother = bs, family = sm.families.Poisson())
    # fit the data
    res_bs = gam_bs.fit()
    # now predict on the full coordinate space
    d_out = pd.DataFrame({"x":x})
    # return the predicted smooth curve
    out_spline = np.asarray(res_bs.predict(d_out, exog_smooth = d_out))
    return out_spline

def fractional_1D(array):
    """
    Function to scale a 1D signal based on the total signal

    Args:
        array - 1 dimensional numpy array
    Returns:
        outarray - (array/sum(array))*len(array)
    """
    return (array/np.nansum(array))*len(array)

def min_max_1D(array):
    """
    Function to scale a 1D signal to a 0,1 scale

    Args:
        array - 1 dimensional numpy array
    Returns:
        outarray - min maxed array of the same size

    """
    return (array - np.nanmin(array))/(np.nanmax(array) - np.nanmin(array))

def weighted_center(array, only_finite = True, normalize = False):
    """
    Find the weighted center of a 1D array

    Args:
        array - 1 dimensional numpy array
        only_finite - subset the array to only include finite data points?
        normalize - divide by 1 to give a relative center?
    """
    if only_finite:
        finite = np.isfinite(array)
    else:
        finite = np.ones(array.shape(), dtype = "bool")
    locs = np.arange(len(array))
    weighted_mean = (locs[finite] * array[finite]).sum() / array[finite].sum()

    if normalize:
        return weighted_mean / len(array)
    else:
        return weighted_mean

def relative_summit_loc(array, wsize = 50):
    smoothed = smooth_1D(array, wsize, kernel_type = "gaussian", edge = "mirror")
    peak = np.nanargmax(smoothed)
    return peak

def traveling_ratio(array, wsize = 50, wA = None, wB = None, out = "ratio"):
    # if peak isn't specified then dynamically find it
    if wA is None:
        wA = relative_summit_loc(array, wsize)
    # peak should at the very least be in the first half of the region
    if wA >= len(array) / 2:
        return np.nan
    wA_avg = np.nanmean(array[max(wA - wsize, 0):min(wA + wsize, len(array))])

    if wB is None: 
        # if not defined make wB halfway between the end of the array and the
        # center of window A
        wB = int((len(array) + wA)/2)

    # wB should be far enough away that windows don't overlap; But also the
    # window should not go off the end of the array
    if wB - wsize < wA + wsize or wB + wsize > len(array):
        return np.nan

    wB_avg = np.nanmean(array[(wB - wsize):(wB + wsize)])
    if out == "ratio":
        out_val = wB_avg/wA_avg
    elif out == "A":
        out_val = wA_avg
    elif out == "B":
        out_val = wB_avg
    else:
        raise ValueError("out Must be ratio, A, or B")

    return out_val

def to_lower_resolution(array, res_scale_factor = 5, summary_func = np.nanmean): 
    size = len(array)
    if res_scale_factor >= size:
        raise ValueError("Attempt to go to lower resolution than entire array array size %s res factor %s"%(size, res_scale_factor))
    outarray = []

    for index in np.arange(0, size, res_scale_factor):
        end = index + res_scale_factor
        if end >= size:
            end = size
        outarray.append(summary_func(array[index:end]))
    return np.array(outarray)
