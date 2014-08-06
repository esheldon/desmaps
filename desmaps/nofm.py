from __future__ import print_function
import numpy

def make_nofm(mags, nbin=200, min_mag=14., max_mag=30.):
    """
    make N(mag) from the input file.

    parameters
    ----------
    mags: array
        array of magnitudes
    nbin: int
        number of bins
    min_mag: float
        min magnitude in histogram
    max_mag: float
        max mag in histogram

    returns
    -------
    the output is a numpy array with fields.
    
    It has the follow fields, 
        'center' center of the mag bin
        'low'    lower part of the mag bin
        'high'   upper part of the mag bin
        'num'    number of objects that fell in the bin
    """
    import esutil as eu

    bs=eu.stat.Binner(mags)
    bs.dohist(nbin=nbin, min=min_mag, max=max_mag)
    
    bs.calc_stats()

    dtype=[('center','f8'),
           ('low','f8'),
           ('high','f8'),
           ('num','i8')]
    output=numpy.zeros(bs['center'].size, dtype=dtype)

    output['center'] = bs['center']
    output['low'] = bs['low']
    output['high'] = bs['high']
    output['num'] = bs['hist']

    return output

def make_nofm_file(input_file,
                   output_file,
                   mag_name='mag_auto',
                   band=2,
                   nbin=100,
                   min_mag=14.,
                   max_mag=25.):
    """
    generate n(m) from the input file and write to the output file

    see make_nofm for format of the output data

    parameters
    ----------
    input_file: string
        file holding mags
    output_file: string
        file to hold the histogram        
    mag_name: string
        name of mag column
    band: int
        band number
    nbin: int
        number of bins
    min_mag: float
        min magnitude in histogram
    max_mag: float
        max mag in histogram
    """
    import fitsio

    print("reading:",input_file)
    data=fitsio.read(input_file, lower=True)

    output = make_nofm(data[mag_name][:,band])

    print("writing output file:",output_file)
    fitsio.write(output_file, output, clobber=True)


