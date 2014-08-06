from __future__ import print_function
import numpy

_DEFAULT_NINTERP=200
_DEFAULT_NINT=100

def get_neff(nofm, compobj, min_mag, max_mag,
             ninterp=_DEFAULT_NINTERP, nint=_DEFAULT_NINT):
    """
    compute the effective number of galaxies given the input true N(m) as
    computed by desmap.make_nofm() and a desmap.Completeness() object

    parameters
    ----------
    nofm: array with fields
        An array with fields as computed using desmap.make_nofm
    compobj: Completeness
        A Completeness object for a given n-sigma magnitude limit.
    min_mag: float
        minimum mag over which to integrate
    max_mag: float
        maximum mag over which to integrate
    ninterp: int
        number of points to interpolate the nofm and to evaluate the
        completeness
    nint: int
        number of points to use in the Gauss-Legendre integration

    examples
    --------
    import desmap
    nofm=fitsio.read('mag-auto-i-hist.fits')
    maglim=23.7 # e.g. a 10-sigma limit
    compobj=desmap.Completeness(maglim)

    min_mag=16.0
    max_mag=maglim # integrate to the same 10-sigma limit
    neff = desmap.get_neff(nofm, compobj, min_mag, max_mag)
    """
    from esutil.integrate import QGauss
    #import biggles

    interp_mag = numpy.linspace(min_mag, max_mag, ninterp)

    nofm_interp = numpy.interp(interp_mag, nofm['center'], nofm['num'])

    #biggles.plot(interp_mag, nofm_interp, xtitle='m interp values', ytitle='n(m) interp')

    comp = compobj(interp_mag)
    
    qg=QGauss(nint)

    ntotal = qg.integrate(interp_mag, nofm_interp)
    neff = qg.integrate(interp_mag, nofm_interp*comp)

    #print("n_total:",ntotal)
    #print("n_eff:",neff)

    return neff

def tabulate_neff(nofm,
                  maglims,
                  min_int_mag=None,
                  max_int_mag=None,
                  ninterp=_DEFAULT_NINTERP,
                  nint=_DEFAULT_NINT):
    """
    tabulate neff for the set of n-sigma maglim, integrating
    n(m)*completeness(m)

    parameters
    ----------
    nofm: array with fields
        array produced by desmaps.make_nofm
    maglims: array
        set of n-sigma magnitude limits
    min_int_mag: array or None, optional
        minimum value over which to integrate, default 16
        for all.
    max_int_mag: array or None, optional
        max value over which to integrate mags.  If not input,
        maglim is used
    ninterp: int
        number of points to interpolate the nofm and to evaluate the
        completeness
    nint: int
        number of points to use in the Gauss-Legendre integration
    """
    from .completeness import Completeness 
    maglims=numpy.array(maglims, ndmin=1, dtype='f8')

    if min_int_mag is None:
        min_int_mag = numpy.array([16.]*maglims.size, dtype='f8')
    if max_int_mag is None:
        max_int_mag = maglims.copy()

    neff_vals=numpy.zeros(maglims.size)

    for i in xrange(maglims.size):
        compobj = Completeness(maglims[i])
        neff_vals[i] = get_neff(nofm,
                                compobj,
                                min_int_mag[i],
                                max_int_mag[i],
                                ninterp=ninterp,
                                nint=nint)

    dtype=[('maglim','f8'),
           ('neff','f8')]
    output=numpy.zeros(maglims.size, dtype=dtype)
    output['maglim'] = maglims
    output['neff'] = neff_vals
    return output
 

