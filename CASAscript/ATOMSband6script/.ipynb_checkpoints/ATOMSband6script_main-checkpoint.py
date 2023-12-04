#########################################
#by Xunchuan Liu 
#initialized at 2023.03.25
#to clean the ATOMS band 6 data
#########################################

import sys
if '.' not in sys.path:
    sys.path.append('.')
import numpy as np
import ATOMSband6script_getwin as getwin
import os
import glob
import shutil
import astropy.table as t
import ATOMSband6script_initsourceconfig

def run(field,
        doconcatfile=False,
        docontclean = False,
        docubeclean = False,
        sourceconfig = None,
        infocsv='ATMOSband6_sourceinfo.csv'
        ncore = 50,
        solint = 'inf'，
        DELETE = False，
        datadir = './data/',
        mpiscriptdir = './mpicasarun/',
       )

    """
    Parameters:
        DELETE: default=FALSE, be careful, a True value will lead to deleting previous results.
                if you are not so sure, set it as False.
    """

    #initialize sourceconfig
    if sourceconfig is not None:
        sourceconfig = ATOMSband6script_initsourceconfig.get_sourceconfig(field, infocsv)
    vlsr = sourceconfig['vlsr']
    wideband_outids = sourceconfig['wideband_outids']
    inmsfiles = sourceconfig['inmsfiles']
    phasethreshold = sourceconfig['phasethreshold']
    ampthreshold = sourceconfig['ampthreshold']
    contthreshold = sourceconfig['contthreshold']
    linethreshold = sourceconfig['linethreshold'] 
    doselfcal = sourceconfig['doselfcal']
    linewinmethod = sourceconfig['linewinmethod'] 
    #solint = "40s"
    done = sourceconfig['done']

    if not os.path.exists('data'):
        os.system('mkdir data')
    if os.path.exists('data/%s' %field):
        if DELETE:  #be careful
            if done:
                print("This source has been cleaned!")
                print("If you indeed want to reclean it, set 'done' as 0 for %s in the .csv configure file" %field)
                return 
            else:
                os.system('rm -r data/%s' %field)
                os.system('mkdir data/%s' %field)
    else:
        os.system('mkdir data/%s' %field)
    if not os.path.exists(mpiscriptdir):
        os.system('mkdir %s' %mpiscriptdir)

    #split out the uv data from the original ms files
    if True and doconcatfile:
        for inmsfile in inmsfiles:
            selectspwids = sourceconfig['inmsfilespwids'][inmsfile]
            inmsfilepath = sourceconfig['inmsfilespath'][inmsfile]
            outputvis = os.path.join(datadir,field,inmsfile+'.ms')
            #print(inmsfile,selectspwids)
            split( vis=inmsfilepath, 
                   spw=selectspwids,      
                   field=field, 
                   outputvis=outputvis,
                   datacolumn='corrected')

    #save before flag
    if True and doconcatfile:
        for inmsfile in inmsfiles:
            inmsfilepath =  os.path.join(datadir,field,inmsfile+'.ms')
            flagmanager(vis=inmsfilepath,
                        mode='save',
                        versionname='before_cont_flags')   

    #concat file, get the spwids in the concat file for ith widebands
    if True and doconcatfile:
        inmsfilepaths =  [os.path.join(datadir,field,inmsfile+'.ms') for inmsfile in inmsfiles]
        concatunflagfilepath = os.path.join(datadir,field,'concat_unflag.ms')
        concat(vis=inmsfilepaths,
           concatvis=concatunflagfilepath)  

    concatunflagfilepath = os.path.join(datadir,field,'concat_unflag.ms')
    if True and os.path.exists(concatunflagfilepath):
        ms.open(concatunflagfilepath)
        spwinfo=ms.getspectralwindowinfo()
        allspwids = [int(i) for i in spwinfo.keys()]
        allspwids.sort()
        allspwids = ['%i'%i for i in allspwids]
        for i in wideband_outids:
            sourceconfig['widebandsinconcatfile'][i] = []
            outfreq = sourceconfig['widebandsfreq'][i]
            for j in allspwids:
                infreq =  spwinfo[j]['Chan1Freq'] + spwinfo[j]['ChanWidth']*spwinfo[j][ 'NumChan' ]/2.
                infreq = infreq/1E6
                if np.abs(infreq-outfreq)<1000:
                    sourceconfig['widebandsinconcatfile'][i].append('%i' % (spwinfo[j]['SpectralWindowId'])  )

    #get flagged concat ms file. unflag the input ms files.
    if True and doconcatfile:
        for inmsfile in inmsfiles:
            inmsfilepath =  os.path.join(datadir,field,inmsfile+'.ms')
            assert ms.open(inmsfilepath)
            wins = []
            for i in wideband_outids:
                f=ms.cvelfreqs(spwids=int(i), 
                               mode='channel', 
                               width=0, 
                               outframe='LSRK')/1E6
                if linewinmethod==1:
                    _,win=getwin.getwin(f, inverse=False, velo=vlsr)
                elif linewinmethod==2:
                    _,win=getwin.getwin_version2(f, field, i,sourceconfig, inverse=False)
                else:
                    raise Excpetion("not matched linewinmethod, just 1 or 2 can be accepted")
                win = ';'.join(['%i~%i' %(i[0],i[1]) for i in win])
                wins.append('%s:%s' %(i,win))
            winstr = ','.join(wins)
            flagdata(vis=inmsfilepath,mode='manual',spw=winstr,flagbackup=False)
        inmsfilepaths =  [os.path.join(datadir,field,inmsfile+'.ms') for inmsfile in inmsfiles]
        concatflagfilepath = os.path.join(datadir,field,'concat_flag.ms')
        concat(vis=inmsfilepaths,
           concatvis=concatflagfilepath)    
        for inmsfile in inmsfiles:
            inmsfilepath =  os.path.join(datadir,field,inmsfile+'.ms')
            flagmanager(vis=inmsfilepath,
                        mode='restore',
                        versionname='before_cont_flags')  

    #round 1 continuum tclean.
    if True and docontclean:
        round1threshold = (phasethreshold if doselfcal else contthreshold)
        concatflagfilepath = os.path.join(datadir,field,'concat_flag.ms')
        outimagename = os.path.join(datadir,field,'concat_flag_continuum')
        field = sourceconfig['field']
        tclean(vis = concatflagfilepath,
          field=field,
          imagename = outimagename,
          restart=True,
          specmode = 'mfs',
          deconvolver = 'mtmfs',
          nterms = 2,
          scales = [0,5,15,50,150],
          imsize = [900,900], 
          cell= '0.05arcsec', 
          gridder = 'mosaic',
          mosweight = True,
          usepointing = False,
          weighting = 'briggs', 
          robust = 0.5,
          niter = 10000, 
          threshold = round1threshold, 
          interactive = False,
          pbcor = True,
          pblimit = 0.2,
          nmajor=5,
          datacolumn = 'data', #the first round always starts from the data column
          #calcres=True,
          savemodel = 'modelcolumn',
          ) 

    #do phase and phase/amp selfcal of continuum
    refants = []
    if ('TM1' in inmsfiles) or ('TM2' in inmsfiles):
        refants.append('DA47')
    if 'ACA' in inmsfiles:
        refants.append('CM01')
    refant = ','.join(refants)
    caltable_round1 =   os.path.join(datadir,field,'phase_1.cal')
    caltable_round2 =   os.path.join(datadir,field,'phase_2.cal')
    if doselfcal and docontclean:
        field = sourceconfig['field']
        concatflagfilepath = os.path.join(datadir,field,'concat_flag.ms')

        gaincal(vis=concatflagfilepath,
            caltable=caltable_round1,
            solint=solint,
            calmode="p",
            refant=refant,
            gaintype="G")

        applycal(vis=concatflagfilepath,
             gaintable=[caltable_round1],
             interp="linear",
             applymode='calonly',
             )

        outimagename = os.path.join(datadir,field,'concat_flag_continuum_round2')
        tclean(vis = concatflagfilepath,
          field=field,
          imagename = outimagename,
          restart=True,
          specmode = 'mfs',
          deconvolver = 'mtmfs',
          nterms = 2,
          scales = [0,5,15,50,150],
          imsize = [900,900], 
          cell= '0.05arcsec', 
          gridder = 'mosaic',
          mosweight = True,
          usepointing = False,
          weighting = 'briggs', 
          robust = 0.5,
          niter = 10000, 
          threshold = ampthreshold, 
          interactive = False,
          pbcor = True,
          pblimit = 0.2,
          nmajor=6,
          savemodel = "modelcolumn",
          )  

        gaincal(vis=concatflagfilepath,
            caltable=caltable_round2,
            gaintable=[caltable_round1],
            solint=solint,
            calmode="ap",
            refant= refant,
            gaintype="G")    


        applycal(vis=concatflagfilepath,
             gaintable=[caltable_round1,caltable_round2],
             interp="linear",
             applymode='calonly',
             )

        outimagename = os.path.join(datadir,field,'concat_flag_continuum_round3')
        tclean(vis = concatflagfilepath,
          field=field,
          imagename = outimagename,
          restart=True,
          specmode = 'mfs',
          deconvolver = 'mtmfs',
          nterms = 2,
          scales = [0,5,15,50,150],
          imsize = [900,900], 
          cell= '0.05arcsec', 
          gridder = 'mosaic',
          mosweight = True,
          usepointing = False,
          weighting = 'briggs', 
          robust = 0.5,
          niter = 10000, 
          threshold = contthreshold, 
          interactive = False,
          pbcor = True,
          pblimit = 0.2,
          nmajor=8,
          savemodel = "modelcolumn",
          )     

    #move the continuum tclean result to *_final*
    if True and docontclean:
        if doselfcal:
            prestr = os.path.join(datadir,field,'concat_flag_continuum_round3*')
        else:
            prestr = os.path.join(datadir,field,'concat_flag_continuum.*')
        for s in glob.glob(prestr):
            thesplit =  os.path.split(s)[-1].split('.')
            base = thesplit[0]
            ext = '.'.join(thesplit[1:])
            outpath = os.path.join(datadir ,field,'.'.join( ['concat_flag_continuum_final' ,ext]))
            shutil.move(s,outpath)
        exportfits(os.path.join(datadir,field,'concat_flag_continuum_final.image.tt0.pbcor'),
                   os.path.join(datadir,field,'concat_flag_continuum_final.image.tt0.pbcor.fits'),
                   overwrite=True,
                )

    #apply the continuum selcal gaintable to all channels
    if True and docubeclean:
        concatunflagfilepath = os.path.join(datadir,field,'concat_unflag.ms')
        if doselfcal:
                applycal(vis=concatunflagfilepath,
                     gaintable=[caltable_round1,caltable_round2],
                     interp="linear",
                     applymode='calonly',
                     )

    #subtract the continuum in uv domain
    if True and docubeclean:
        concatunflagfilepath = os.path.join(datadir,field,'concat_unflag.ms')
        concatcontsubfilepath = os.path.join(datadir,field,'concat_contsub.ms')
        ms.open(concatunflagfilepath)
        wins = []
        for i in wideband_outids:    
            for j in sourceconfig['widebandsinconcatfile'][i]:
                f=ms.cvelfreqs(spwids=int(j), 
                               mode='channel', 
                               width=0, 
                               outframe='LSRK')/1E6
                if linewinmethod==1:
                    _,win=getwin.getwin(f, inverse=True, velo=vlsr)
                elif linewinmethod==2:
                    _,win=getwin.getwin_version2(f, field, i,sourceconfig, inverse=True)
                else:
                    raise Excpetion("not matched linewinmethod, just 1 or 2 can be accepted")
                win = ';'.join(['%i~%i' %(k[0],k[1]) for k in win])
                wins.append('%s:%s' %(j,win))
        winstr = ','.join(wins)          
        uvcontsub(vis=concatunflagfilepath,
                  outputvis=concatcontsubfilepath,
                  fitspec=winstr,
                  fitorder=1,
                  datacolumn='corrected' if doselfcal else 'data',
                 )

    if False and docubeclean:
       concatunflagfilepath = os.path.join(datadir,field,'concat_unflag.ms')
       tclean(vis = concatunflagfilepath,
          field=field,
          imagename =  os.path.join(datadir ,field,'concat_flag_continuum_final'),
          restart=True,
          specmode = 'mfs',
          deconvolver = 'mtmfs',
          nterms = 2,
          scales = [0,5,15,50,150],
          imsize = [900,900], 
          cell= '0.05arcsec', 
          gridder = 'mosaic',
          mosweight = True,
          usepointing = False,
          weighting = 'briggs', 
          robust = 0.5,
          niter = 0, 
          threshold = '', 
          interactive = False,
          pbcor = False,
          pblimit = 0.2,
          nmajor=5,
          calcres=False,
          calcpsf=False,
          savemodel = "modelcolumn",
          )   

    #split out ms file for each wideband
    if True and docubeclean:
        concatcontsubfilepath = os.path.join(datadir,field,'concat_contsub.ms')
        for i in wideband_outids: 
            contsubfilepath = os.path.join(datadir,field,'concat_contsub_spw%s.ms' %i)
            thespws = ','.join(sourceconfig['widebandsinconcatfile'][i])
            split(vis=concatcontsubfilepath, 
                  spw=thespws,      
                  outputvis=contsubfilepath,
                  datacolumn = 'data',
                  )


    #create the script for cleaning the cube under mpicasa
    #the script can be run through:
    #    mpicasa -n 50 casa --nogui --nologger -c thescript.py
    #Here, 50 is the number of processes which can be changed according the cpu cores, 
    cleanscripttemp = """tclean(vis='{vis}',
        restart = True,
        calcres = {calcres},
        calcpsf = {calcpsf},
        field = '{field}',
        imagename='{imagename}',
        datacolumn='corrected',
        spw='{spw}',
        specmode='cube', 
        width='{width}',
        start='{start}',
        outframe='LSRK',
        nchan={nchan},
        threshold='{threshold}', 
        imsize=[900, 900], 
        cell=['0.05arcsec'], 
        niter={niter}, 
        deconvolver='multiscale',
        scales =[0,5,15,45], 
        gridder='mosaic', 
        weighting='briggs',
        robust=0.5,
        pbcor=False, 
        pblimit=0.2, 
        #restoringbeam='common',
        usemask = 'auto-multithresh',
        sidelobethreshold = 2.0,
        noisethreshold = 4.25,
        minbeamfrac = 0.3,
        lownoisethreshold = 1.5, 
        #chanchunks=-1, 
        perchanweightdensity = True,
        interactive=False,
        parallel=True,
        cycleniter={cycleniter},
        nmajor={nmajor},
        negativethreshold=6,
        )
    """
    if True and docubeclean:
        for i in wideband_outids:
            contsubfilepath = os.path.join(datadir,field,'concat_contsub_spw%s.ms' %i)
            theimagename = os.path.join(datadir,field,'concat_contsub_spw%s_cube' %i)
            fmin = None; fmax = None
            ms.open(contsubfilepath)
            spwinfo=ms.getspectralwindowinfo()
            df = None
            for j in spwinfo:
                thef = ms.cvelfreqs(spwids=int(j), mode='channel', 
                             width=0, outframe='LSRK')/1E6   
                if df is None:
                    df = np.abs(thef[1]-thef[0])
                else:
                    df = min(  np.abs(thef[1]-thef[0])  ,df )
                if fmin is None:
                    fmin = thef.min()
                else:
                    fmin = max(fmin, thef.min()) 
                if fmax is None:
                    fmax = thef.max()
                else:
                    fmax = min(fmax, thef.max()) 
            df = max(df, 0.488281)
            thewidth = '%fMHz' %(df)
            thenchan = (fmax-fmin)//df
            assert (thenchan>3000) and (thenchan<5000) 
            cleanscript = cleanscripttemp.format(
                vis=os.path.join('../',contsubfilepath),
                calcres = 'True',
                calcpsf = 'True',
                field = sourceconfig['field'],
                imagename = os.path.join('../',theimagename),
                spw="",
                width=thewidth,
                start='%fMHz' %fmin,
                nchan = '%i' %thenchan,
                threshold = '30mJy',
                niter='6000000',
                cycleniter='500',
                nmajor='6',
            ) 
            scriptfile = os.path.join(mpiscriptdir,'mpiscript_spw%s' %i)
            with open(scriptfile,'w') as thefile:
                thefile.write(cleanscript)

    if True and docubeclean:
        mpiscriptrun=os.path.join(mpiscriptdir,'mpiscriptrun.sh')
        with open(mpiscriptrun,'w') as thef:
            for i in wideband_outids:
                scriptfile = 'mpiscript_spw%s' %i
                thef.write("mpicasa -n {ncore} casa --nogui --nologger -c {script}\n".format(script=scriptfile, ncore=ncore))
        #open a terminal to execute the mpiscripts, e.g.
        #os.system('gnome-terminal --working-directory="%s" -- bash -c "cat mpiscriptrun.sh"' %(os.path.abspath(mpiscriptdir), ) )
        #However, be aware that gnome-terminal --disable-factory may do not support to open a blocked new terminal
        #so, I use xfce4-terminal instead. It assumes that xfce4-terminal has been installed
        #add --hold if you want to keep the terminal around after the executed command finishes executing
        os.system('xfce4-terminal --working-directory="%s" --command "bash mpiscriptrun.sh"' %(os.path.abspath(mpiscriptdir), ) )

    if True and docubeclean:
        for i in wideband_outids:
            theimagename = os.path.join(datadir,field,'concat_contsub_spw%s_cube.image' %i)
            thepb = os.path.join(datadir,field,'concat_contsub_spw%s_cube.pb' %i)
            thepbcor = os.path.join(datadir,field,'concat_contsub_spw%s_cube.image.pbcor' %i)
            thepbcorfits = os.path.join(datadir,field,'concat_contsub_spw%s_cube.image.pbcor.fits' %i)
            impbcor(imagename=theimagename,
                    pbimage=thepb,
                    outfile=thepbcor,
                    overwrite=True,
                   )
            exportfits(imagename=thepbcor,
                       fitsimage=thepbcorfits,
                       overwrite=True, 
                      )
