from __future__ import print_function
import desmaps
import fitsio
import numpy

def example1():
    pars={'nofm_nbin':200,
          'nofm_min_mag': 14.0,
          'nofm_max_mag': 30.0,

          'maglim_min':19.0,
          'maglim_max':24.5,
          'nmaglim':100}

    catname='/astro/u/esheldon/masks/des/mask-cat/combined_masked_all_cat.fit'
    print("reading:",catname)
    catalog=fitsio.read(catname,lower=True)

    # nofm for i-band

    print("making nofm")
    nofm=desmaps.make_nofm(catalog['mag_auto'][:,2],
                           nbin=pars['nofm_nbin'],
                           min_mag=pars['nofm_min_mag'],
                           max_mag=pars['nofm_max_mag'])

    # 10-sigma magnitude limits
    maglims=numpy.linspace(pars['maglim_min'], pars['maglim_max'], pars['nmaglim'])

    # get effective number including completeness for the 10-sigma limit
    # will integrate [16,maglim] unless specify min_int_mag and max_int_mag

    print("tabulating neff vs maglim")
    neffdata=desmaps.tabulate_neff(nofm, maglims)

    # now we can interpolate a maglimit healpix map to get neff
    hmap='/astro/u/esheldon/masks/des/sva1-gold/sva1_gold_1.0_nside4096-64_nest_i_weights.fits'
    print("reading maglim map:",hmap)
    maglim_map,hdr = fitsio.read(hmap, columns='I', header=True)
    maglim_map=maglim_map.ravel()

    print("interpolating to get neff for maglim map")
    w,=numpy.where(maglim_map > 0)
    neff_interpolated = numpy.interp(maglim_map[w], maglims, neffdata['neff'])


    # convert to strange healpix standard file format
    print("creating output")
    nrows=maglim_map.size/1024
    neff_map = numpy.zeros(nrows, dtype=[('i','f4',1024)])
    neff_flat = neff_map['i'].ravel()
    neff_flat[w] = neff_interpolated.astype('f4')

    pars['ordering'] = hdr['ordering'].strip()
    pars['nside']    = hdr['nside']

    neff_file="sva1_gold_1.0_nside4096-64_nest_i_neff.fits"
    print("writing:",neff_file)
    fitsio.write(neff_file, neff_map, header=pars, clobber=True)

    # write some data for posterity
    nofm_file='mag-auto-i-hist.fits'
    print("nofm_file:",nofm_file)
    fitsio.write(nofm_file, nofm, header=pars)
    neff_vs_maglim_file='neff-vs-maglim.fits'
    print("writing neff vs maglim:",neff_vs_maglim_file)
    fitsio.write(neff_vs_maglim_file, neffdata, header=pars)

if __name__=="__main__":
    example1()

