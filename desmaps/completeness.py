class Completeness(object):
    """
    class to compute the completeness as a function of magnitude for
    the given n-sigma completeness limit

    completeness = (eff/2)*(1-erf((mag - m50)/sqrt(2*w)))

    parameters
    ----------
    maglim: float
        n-sigma completeness magnitude limit
    w: float
        width of the error function
    eff: float
        completeness at the bright end.  Doesn't matter for relative
        density calculations

    examples
    --------
    import desmaps
    maglim=23.7
    cobj=desmaps.Completeness(maglim)
    comp = cobj(mags)
    """
    def __init__(self, maglim, w=0.193, eff=0.97):
        self.maglim=maglim
        self.w=w
        self.eff=eff
        self.m50 = maglim2mag50(maglim)

    def __call__(self, mag):
        from numpy import sqrt
        from scipy.special import erf
        return (self.eff/2)*(1-erf((mag - self.m50)/sqrt(2*self.w)))


def maglim2mag50(maglim):
    """
    convert the magnitude n-sigma limit to a m50 completeness limite

    this currently works for i-band

    parameters
    ----------
    maglim: float
        n-sigma completeness magnitude limit

    returns
    -------
    m50: float
        magnitude where the sample is 50 percent complete
    """
    m50 = 23.81 + 0.78*(maglim-23)

    return m50

