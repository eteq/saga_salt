from __future__ import print_function

import os, sys, string, shutil

import numpy as np
from matplotlib import pyplot as plt

from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy.constants import c

sys.path.insert(1, os.path.join(os.path.dirname(__file__), 'zSALT'))
import zsalt


def load_salt_ascii_spec(fn):
    import specutils
    tab = Table.read(fn, format='ascii', names=['wl', 'spec', 'unc'])
    return specutils.Spectrum1D(spectral_axis=tab['wl']*u.angstrom,
                                flux=tab['spec']*u.count,
                                uncertainty=StdDevUncertainty(tab['unc']**0.5))


def plot_spectrum(spec1d):
    plt.plot(spec1d.spectral_axis, spec1d.flux)
    plt.plot(spec1d.spectral_axis, spec1d.uncertainty.array)
    plt.show()


def extract_target(objfile, objrowtuple, clean=True,  cal_file='', out_fn_path='.'):
    """
    `objfile` is expected to have everything up to sky-subtraction (i.e., rectified)
    """
    from zsalt.salt_extract import extract_spectra

    lower = min(objrowtuple)
    upper = max(objrowtuple)
    half_width = (upper - lower)/2
    center_row = lower + half_width

    # salt_extract(objfile, center_row, half_width, specformat='ascii',
    #              ext=ext, calfile=cal_file, convert=True,
    #              normalize=False, clean_spectra=clean)
    #hdu = pyfits.open(objfile)
    hdu = fits.open(objfile)
    out_fn = os.path.join(out_fn_path, hdu[0].header['OBJECT'] + '_spec.txt')
    try:
        extract_spectra(hdu, center_row, half_width, out_fn, smooth=False,
                        grow=10, clobber=True, specformat='ascii', convert=True,
                        calfile=cal_file, cleanspectra=clean)
    finally:
        hdu.close()

    return out_fn

def template_redshift(infile, z1=-0.000005, z2=.02, zstep=0.000001, plot=None,
                      templs=list(range(23,29)), verbose=True):
    from zsalt.redshift import xcor_redshift, loadtext, loadiraf, loadsdss

    if infile.count('fits'):
        hdu = fits.open(infile)
        spec = loadiraf(hdu)
    else:
        spec = loadtext(infile)

    best_cc = 0
    results = {'templates': templs, 'template_names':[], 'template_paths':[],
               'zs':[], 'vs':[], 'best_ccs':[], 'cc_arrs':[]}
    for i, t in enumerate(results['templates']):
        if verbose:
            print('Redshifting template', i, 'of', len(results['templates']))
        template_name = 'spDR2-0{}'.format(string.zfill(t, 2))
        results['template_names'].append(template_name)
        template_path = os.path.join(os.path.dirname(zsalt.__file__),
                                     'template', template_name + '.fit')
        results['template_paths'].append(template_path)

        template = loadsdss(fits.open(template_path))

        z_arr, cc_arr = xcor_redshift(spec, template, z1=z1, z2=z2, zstep=zstep)
        z = z_arr[cc_arr.argmax()]

        results['zs'].append(z)
        results['vs'].append(z*c.to(u.km/u.s))
        results['cc_arrs'].append(cc_arr)
        results['best_ccs'].append(cc_arr.max())
        if best_cc < results['best_ccs'][-1]:
            best_cc = results['best_ccs'][-1]
            best_z = results['zs'][-1]
            best_idx = i
            best_cc_arr = cc_arr

    results['vs'] = u.Quantity(results['vs'])

    if plot:
        inter = plt.isinteractive()
        plt.ioff()

        def plot_templ(ax1, ax2, spec, templidx):
            template = loadsdss(fits.open(results['template_paths'][templidx]))

            ax1.plot(z_arr, results['cc_arrs'][templidx])
            ax1.axvline(results['zs'][templidx], c='k', ls=':', alpha=.8)
            cflux = np.convolve(spec.flux, np.ones(10), mode='same')
            ax2.plot(spec.wavelength, cflux)
            nflux = np.interp(spec.wavelength, (1+results['zs'][templidx])*template.wavelength, template.flux)
            ax2.plot(spec.wavelength, nflux*cflux.mean()/nflux.mean(), alpha=.7)
            ax1.set_ylabel(results['template_names'][templidx])
            ax2.set_title('z={}, v={}'.format(results['zs'][templidx],
                                              results['vs'][templidx]))

        if plot == 'best':
            fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 5))
            plot_templ(ax1, ax2, spec, best_idx)
        elif plot == 'all':
            nplots = len(results['template_paths'])
            fig, axs = plt.subplots(nplots, 2, figsize=(10, 5*nplots))
            for templ_idx, (ax1, ax2) in enumerate(axs):
                plot_templ(ax1, ax2, spec, templ_idx)
                if templ_idx == best_idx:
                    ax1.set_ylabel(ax1.get_ylabel() + '(BEST)')
        else:
            raise NotImplementedError(plot)

        fig.tight_layout()

        plotfn = infile.replace('.txt', '') + '.png'
        fig.savefig(plotfn)

        if inter:
            plt.ion()

    return best_idx, results, (best_z, best_z*c.to(u.km/u.s), results['template_names'][best_idx])


PYSALT_TEMPL="""
pysalt
saltred
saltspec

saltcrclean('{scifn}',outpref='c', crtype='edge', thresh='10.0', mbox='25',bthresh='5.0', flux_ratio='0.2', bbox='25', gain='1.0', rdnoise='5.0', fthresh='5.0', bfactor='2',gbox='4', maxiter='5', multithread='no', clobber='yes', logfile='salt.log', verbose='yes')

specidentify(images='{arcfn}', linelist='{arclinesfn}', outfile='db.sol', guesstype='rss', guessfile="",automethod='Matchlines', function='legendre', order='3', rstep='200', rstart='middlerow', mdiff='20', thresh='3', niter='5', startext='0', inter='yes', clobber='no', logfile='salt.log', verbose='yes')

specrectify(images='c{scifn}',  outimages="", outpref='x', solfile='db.sol', caltype='line',  function='legendre', order='3', inttype='interp', w1='None', w2='None', dw='None', nw='None', blank='0.0',  clobber='yes', logfile='salt.log', verbose='yes')

# result: xc{scifn}
"""[1:]
POSTPYSALT_TEMPL="""
from __future__ import print_function

import sys
sys.path.insert(1, '{path_to_sagasalt}')

from saga_salt import extract_target, template_redshift

extracted_fn = extract_target('xc{scifn}', (LOWER_OBJ, UPPER_OBJ))
print('Extracted', extracted_fn)
best_idx, results, summary = template_redshift(extracted_fn, plot='all')
print('Redshift summary:, ', summary)
with open('redshift_results', 'w') as f:
    f.write(str(best_idx))
    f.write('\n\n')
    f.write(str(summary))
    f.write('\n\n')
    f.write(str(results))
"""[1:]
def prepare_reduction(sciarcs, pathtoarcs, basepath='.', clobber=False):
    """
    ``sciarcs`` a list of 2-tuples (scifn, arcfn).  Each one gets a reduction
    directory, with symlinks, and a
    text file with the relevant pysalt commands, and corresponding post-pysalt
    python script.
    """
    targetnms = []
    for scifn, arcfn in sciarcs:
        scih = fits.getheader(scifn)
        arch = fits.getheader(arcfn)
        if scih['OBJECT'] == 'ARC' or scih['CCDTYPE'] != 'OBJECT':
            raise ValueError('File ' + str(scifn) + ' is not a sci!')
        if arch['CCDTYPE'] != 'ARC':
            raise ValueError('File ' + str(arcfn) + ' is not an arc!')
        if arch['LAMPID'] == 'NONE':
            raise ValueError('Arc ' + str(arcfn) + ' has no lamps on!')
        targetnms.append(scih['OBJECT'])

        arcpath = os.path.join(pathtoarcs, arch['LAMPID'])
        if not os.path.exists(arcpath):
            print('Warning: ', arcpath, 'does not seem to exist', arcfn, 'probably will not work')

    for tnm, (scifn, arcfn) in zip(targetnms, sciarcs):
        target_dir = os.path.join(basepath, tnm)
        if os.path.exists(target_dir):
            if clobber:
                print(target_dir, 'already exists, removing...')
                shutil.rmtree(target_dir)
            else:
                print(target_dir, 'already exists, skpping...')
                continue

        scitarget = os.path.join(target_dir, os.path.basename(scifn))
        arctarget = os.path.join(target_dir, os.path.basename(arcfn))

        print('Creating', target_dir)
        os.mkdir(target_dir)
        print('Linking', scifn, 'to', scitarget)
        os.symlink(os.path.relpath(scifn, target_dir), scitarget)
        print('Linking', arcfn, 'to', arctarget)
        os.symlink(os.path.relpath(arcfn, target_dir), arctarget)

        pysalt_fn = os.path.join(target_dir, tnm+'_pysalt.txt')
        pyscript_fn = os.path.join(target_dir, tnm+'_post_pysalt.py')
        print('Generating', pysalt_fn, 'and',  pyscript_fn)

        arcpath = os.path.join(pathtoarcs, fits.getheader(arcfn)['LAMPID'] + '.txt')
        templ_substitutions = {'scifn': os.path.basename(scifn),
                               'arcfn': os.path.basename(arcfn),
                               'arclinesfn': os.path.abspath(arcpath),
                               'path_to_sagasalt': os.path.relpath(os.path.dirname(__file__), target_dir)}
        with open(pysalt_fn, 'w') as f:
            f.write(PYSALT_TEMPL.format(**templ_substitutions))
        with open(pyscript_fn, 'w') as f:
            f.write(POSTPYSALT_TEMPL.format(**templ_substitutions))
