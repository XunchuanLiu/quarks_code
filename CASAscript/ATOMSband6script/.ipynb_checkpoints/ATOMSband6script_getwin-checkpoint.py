import numpy as np
import astropy.io.fits as fits
import scipy.interpolate as it
import astropy.table as t
import os
import regions
import spectral_cube as sc
import matplotlib.pyplot as plt
plt.ion()


def setwin(field, sourceconfig,usestronglinemask=True,linesnr=5,showplot=False,expandlength=7):
    vlsr = sourceconfig['vlsr']
    if not os.path.exists(os.path.join('materials','ds9region_meanspec')):
        os.system('mkdir %s' %(os.path.join('materials','ds9region_meanspec')) )
    T=t.Table.read('ATMOSband6_sourceinfo.csv',comment=r'^#.*')
    index = np.where(T['field']==field)[0][0]
    scienceid = T['scienceid'][index]
    groupid = T['groupid'][index]

    ds9file = os.path.join('materials/ds9region_for_spwwin', '%s.reg' %field)
    for outid in sourceconfig['wideband_outids']:
        member = T['TM1'][index].split('/')[0]
        inid = sourceconfig['widebands'][outid]['TM1']
        filepath = os.path.join('../2021.1.00095.S',
                                    scienceid,
                                    groupid,
                                    member,
                                    'product',
                                    '%s.%s_sci.spw%s.cube.I.pbcor.fits' %(member, field, inid)
                                   )
        #os.system("ls %s" %filepath)
        outfitsfile = os.path.join( 'materials','ds9region_meanspec','%s_spw%s.fits' %(field,outid) )
        if not os.path.exists(outfitsfile):
                cube=sc.SpectralCube.read(filepath)
                subcube=cube.subcube_from_regions(regions.Regions.read(ds9file))
                subcube.allow_huge_operations=True
                meanspec=subcube.mean(axis=(1,2))
                meanspec = meanspec*(1*meanspec.unit/('mJy/beam')).cgs.value
                meanspec._set_unit ('mJy/beam')
                meanspec.write(outfitsfile)
                del meanspec
                del subcube
                del cube
        outwinfile = os.path.join( 'materials','ds9region_meanspec','%s_spw%s_win.fits' %(field,outid) )
        if not os.path.exists(outwinfile):
            hdu = fits.open(outfitsfile)[0]
            header = hdu.header
            freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
            freq = freq/1E6
            if usestronglinemask:
                stronglinewin= getwin(freq,
                                  fitsfile='ATOMSband6script_stronglines_win.fits',
                                  velo=vlsr)[0]
            else:
                stronglinewin = np.zeros_like(freq,dtype='bool')
           
            rms = np.nanstd(hdu.data[~stronglinewin])
            linewin = (np.abs(hdu.data)>linesnr*rms) 
            linewin = np.convolve(linewin*1., [1]*expandlength,mode='same')>0
            linewin = linewin | stronglinewin
            rms = np.nanstd(hdu.data[~linewin])
            linewin = (np.abs(hdu.data)>linesnr*rms) 
            linewin = np.convolve(linewin*1., [1]*expandlength,mode='same')>0
            linewin = linewin | stronglinewin
            outwinhdu = fits.PrimaryHDU( (linewin*1.).astype(hdu.data.dtype),
                             header=header
                            )
            outwinhdu.writeto(outwinfile)
            del outwinhdu
    if showplot:
        fig=plt.figure('plotwin')
        plt.clf()
        N = len(sourceconfig['wideband_outids'])
        i = 0
        for outid in sourceconfig['wideband_outids']:
            i = i+1
            ax=plt.subplot( (N-1)//2+1 ,2, i  )
            fitsfile = os.path.join( 'materials','ds9region_meanspec','%s_spw%s.fits' %(field,outid) )
            winfile = os.path.join( 'materials','ds9region_meanspec','%s_spw%s_win.fits' %(field,outid) )
            header = fits.getheader(fitsfile)
            freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
            freq = freq/1E6
            plt.step(freq, fits.getdata(fitsfile),color='C0')
            ax1 = ax.twinx()
            plt.step(freq, fits.getdata(winfile),color='C1')       
            plt.title(field)     
    del T

def parsewin(win,mergelength=11):
    wintemp = np.hstack([[False],win,[False]])
    nxor = wintemp[1:-1] & np.logical_xor(wintemp[:-2], wintemp[2:]   )
    assert np.sum(nxor)%2==0
    nedge = np.arange(len(nxor))[nxor].reshape((-1,2))
    nedge1 = [nedge[0]]
    for i in nedge[1:]:
        if i[0] > nedge1[-1][1]+mergelength:
            nedge1.append(i)
        else:
            nedge1[-1][1] = i[1]
    win2 = np.zeros_like(win)
    for i in nedge:
        win2[i[0]:i[1]+1] = True    
    return win2, nedge

def getwin(fin, #in MHz
          fitsfile='ATOMSband6script_win.fits',
          velo = 71, #vlsr
          threshold = 5,
          convlength = 31,
          mergelength = 11,
          inverse = False,
          ):
    hdum = fits.open(fitsfile)[0]
    win = hdum.data>threshold
    win=np.convolve(win,[1]*convlength,mode='same')>0
    headerm = hdum.header
    fwin = (np.arange(headerm['NAXIS1'])+1-headerm['CRPIX1'])*headerm['CDELT1']+headerm['CRVAL1']
    fwin = fwin/1E6
    fm = fwin*(1-velo/3E5)
    iters=it.interp1d(fm,win,bounds_error=False, fill_value=0)
    win1 = iters(fin)>0
    if inverse:
        win1 = ~win1
    return parsewin(win1,mergelength=mergelength)

def getwin_version2(fin,field,spwid,sourceconfig=None,inverse=False,mergelength=1):
    winfile = os.path.join( 'materials','ds9region_meanspec','%s_spw%s_win.fits' %(field,spwid) )
    if not os.path.exists(winfile):
        assert sourceconfig is not None
        setwin(field,sourceconfig)
    hdu = fits.open(winfile)[0]
    header = hdu.header
    freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
    freq = freq/1E6
    iters=it.interp1d(freq,fits.getdata(winfile),bounds_error=False, fill_value=1) #be careful, ACA has wider spw
    win = iters(fin)>1E-8
    if inverse:
        win = ~win
    return parsewin(win,mergelength=mergelength)
