import numpy as np
import astropy.io.fits as fits
import scipy.interpolate as it
import astropy.table as t
import os
import regions
import spectral_cube as sc
import matplotlib.pyplot as plt
plt.ion()

def getmeanspec(filepath, ds9file, outfitsfile=None, save=True,overwrite=True,
                convertunit=True,newunit='mJy/beam', meanaxis=(1,2)
               ):
    """
    get and save the meanspec of a cube
    """
    cube=sc.SpectralCube.read(filepath)
    subcube=cube.subcube_from_regions(regions.Regions.read(ds9file))
    subcube.allow_huge_operations=True
    meanspec=subcube.mean(axis=meanaxis)
    if convertunit:
        meanspec = meanspec*(1*meanspec.unit/(newunit)).cgs.value
        meanspec._set_unit (newunit)
    if save and (outfitsfile is not None):
        if  (not os.path.exists(outfitsfile)):
            meanspec.write(outfitsfile)
        elif overwrite:
            os.remove(outfitsfile)
            meanspec.write(outfitsfile)
        else:
            pass
    del subcube
    del cube  
    return meanspec  

def getwinfile(hdu,linesnr,outwinfile,stronglinewin=None,Niter=3):
    header = hdu.header
    freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
    freq = freq/1E6
    if stronglinewin is None:
        stronglinewin = np.zeros_like(freq,dtype='bool')

    thedata = hdu.data.copy()
    linewin = np.zeros_like(theta,dtype='bool')
    linewin = linewin | stronglinewin
    thex = np.arange(len(thedata))
    for i in np.arange(Niter):        
        rms = np.nanstd(thedata[ (~linewin) ])
        linewin = (np.abs(thedata)>linesnr*rms) 
        linewin = np.convolve(linewin*1., [1]*expandlength,mode='same')>0
        linewin = linewin | stronglinewin
        thedata = thedata-np.polyval(np.polyfit(thex[~linewin],thedata[~linewin],1),thex)
    outwinhdu = fits.PrimaryHDU( (linewin*1.).astype(hdu.data.dtype),
                             header=header
                            )
    outwinhdu.writeto(outwinfile,overwrite=True)
    del outwinhdu

def setwin(field, 
           sourceconfig,
           usestronglinemask=True,
           linesnr=5,
           Niter = 3,
           showplot=False,
           expandlength=7,
           rawdatadir='../2021.1.00095.S',
           regdir = 'materials/ds9region_for_spwwin',
           spwwindir='materials/ds9region_meanspec',
           usecleanedcubes = False,
           datadir = 'products/version1/',
           overwritewin=False,
           ):
    vlsr = sourceconfig['vlsr']
    if not os.path.exists(spwwindir):
        os.system('mkdir %s' %spwwindir )
    #T=t.Table.read('ATMOSband6_sourceinfo.csv',comment=r'^#.*')
    #index = np.where(T['field']==field)[0][0]
    #scienceid = T['scienceid'][index]
    #groupid = T['groupid'][index]
    scienceid = sourceconfig['scienceid']
    groupid = sourceconfig['groupid']

    ds9file = os.path.join(regdir, '%s.reg' %field)
    for outid in sourceconfig['wideband_outids']:
        member = None
        for seg in sourceconfig['inmsfilespath']['TM1'].split('/'):
            if 'member.uid' in seg:
                member = seg
        assert member is not None
        inid = sourceconfig['widebands'][outid]['TM1'][0]
        if not usecleanedcubes:
            filepath = os.path.join(rawdatadir,
                                    scienceid,
                                    groupid,
                                    member,
                                    'product',
                                    '%s.%s_sci.spw%s.cube.I.pbcor.fits' %(member, field, inid)
                                   )
        else:
            filepath = os.path.join(datadir,field,'concat_contsub_spw%s_cube.image.pbcor.fits' %outid)
        outfitsfile = os.path.join( spwwindir,'%s_spw%s.fits' %(field,outid) )
        if not os.path.exists(outfitsfile):
            getmeanspec(filepath, ds9file, outfitsfile)
        outwinfile = os.path.join( spwwindir,'%s_spw%s_win.fits' %(field,outid) )
        if (not os.path.exists(outwinfile)) or overwritewin:
            hdu = fits.open(outfitsfile)[0]
            if usestronglinemask:
                header = hdu.header
                freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
                freq = freq/1E6
                stronglinewin= getwin(freq,
                                  fitsfile='ATOMSband6script_stronglines_win.fits',
                                  velo=vlsr)[0]
            else:
                stronglinewin = None
            getwinfile(hdu=hdu, linesnr=linesnr, outwinfile=outwinfile, stronglinewin=stronglinewin, Niter=Niter)
    if showplot:
        fitsfiles = []
        winfiles = []
        for outid in sourceconfig['wideband_outids']:
            fitsfiles.append(os.path.join( spwwindir,'%s_spw%s.fits' %(field,outid) ))
            winfiles.append(os.path.join( spwwindir,'%s_spw%s_win.fits' %(field,outid) ))
        _showplot(fitsfiles,winfiles,title=field)   

def parsewin(win,mergelength=11):
    wintemp = np.hstack([[False],win,[False]])
    nxor = wintemp[1:-1] & np.logical_xor(wintemp[:-2], wintemp[2:]   )
    assert np.sum(nxor)%2==0
    nedge = np.arange(len(nxor))[nxor].reshape((-1,2))
    nedge1 = [nedge[0]]
    #print(win[:100],nedge)#debug#
    for i in nedge[1:]:
        if i[0] > nedge1[-1][1]+mergelength:
            nedge1.append(i)
        else:
            nedge1[-1][1] = i[1]
    win2 = np.zeros_like(win)
    for i in nedge1:
        win2[i[0]:i[1]+1] = True    
    return win2, nedge1

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
    iters=it.interp1d(fm,win,bounds_error=False, fill_value=1)
    win1 = iters(fin)>0
    if inverse:
        win1 = ~win1
    return parsewin(win1,mergelength=mergelength)

def getwin_version2(fin,field,spwid,sourceconfig=None,inverse=False,mergelength=1,spwwindir='materials/ds9region_meanspec'):
    winfile = os.path.join( spwwindir,'%s_spw%s_win.fits' %(field,spwid) )
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

def _showplot(fitsfiles,winfiles,title):
        fig=plt.figure('plotwin')
        plt.clf()
        N = len(fitsfiles)
        for i in np.arange(N):
            ax=plt.subplot( (N-1)//2+1 ,min(2,N), i+1  )
            fitsfile = fitsfiles[i]
            winfile = winfiles[i]
            header = fits.getheader(fitsfile)
            freq = (np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CDELT1']+header['CRVAL1']
            freq = freq/1E6
            plt.step(freq, fits.getdata(fitsfile),color='C0')
            ax1 = ax.twinx()
            plt.step(freq, fits.getdata(winfile),color='C1')       
            plt.title(title)  
