import sys
if '.' not in sys.path:
    sys.path.append('.')
import os
from casatasks import *
import casatools.ms
ms = casatools.ms()
tb = casatools.table()
import glob
import shutil
import numpy as np
import astropy.table as t

Ts=t.Table.read(os.path.join(os.path.dirname(__file__),'sourcetable.csv'))

defaultTM1dirs = [
"member.uid___A001_X1590_X341e/calibrated/uid___A002_Xfcf230_X2e1.ms",
"member.uid___A001_X1590_X3417/calibrated/uid___A002_Xfd02fc_X1154.ms",
"member.uid___A001_X1590_X33d1/calibrated/uid___A002_X10536e5_Xe009.ms.split.cal",
"member.uid___A001_X1590_X33fb/calibrated/merge.ms",
"member.uid___A001_X1590_X3402/calibrated/uid___A002_Xfce28e_Xe8c.ms",
"member.uid___A001_X1590_X3410/calibrated/uid___A002_Xfd02fc_Xdc5.ms",
"member.uid___A001_X1590_X3425/calibrated/uid___A002_Xfcf230_X658.ms",
"member.uid___A001_X1590_X342c/calibrated/uid___A002_Xfd02fc_X532e.ms",
"member.uid___A001_X1590_X3433/calibrated/uid___A002_Xfd02fc_X621.ms.split.cal",
]

def splitdata(rowdataroot='../../2021.1.00095.S/',
        line='SO2_234187',
        freqrange=[234156,234223], #MHz in restfreq at LSR (has not been corrected for VLSR)
        TM1dirs=[],
        fields = '',
        overwrite=False,
       ):
    if (TM1dirs == '') or ( (type(TM1dirs)==list) and (len(TM1dirs)==0)):
        TM1dirs = defaultTM1dirs
    elif type(TM1dirs) == str:
        TM1dirs = [TM1dirs]
    assert type(fields) == str
    if len(fields) == 0:
        fields = []
    else:
        fields = [i.strip() for i in fields.strip().split(',')]

    for TM1dir in TM1dirs:
        TM1fulldir =  glob.glob(os.path.join(rowdataroot,'*/*',TM1dir))[0]
        #print(TM1fulldir)
        assert ms.open(TM1fulldir)
        fieldsOfMS = [] 
        summary = ms.summary()
        for k in summary.keys():
            if ('scan_' in k) and (k[:5] == 'scan_'):
                field = summary[k]['0']['FieldName']
                if (field not in fieldsOfMS) and (field[0]=='I'):
                    fieldsOfMS.append(field)
        if len(fields)>0:
            exefieldsOfMS = [i for i in fieldsOfMS if i in fields]
        else:
            exefieldsOfMS = fieldsOfMS
        #print(TM1fulldir,exefieldsOfMS)  
        
        for field in exefieldsOfMS:
            vlsr=Ts['vlsrnew'][list(Ts['source']).index(field)]
            _splitdata(field=field,TM1fulldir=TM1fulldir,
                       line=line,freqrange=freqrange,overwrite=overwrite,
                       vlsr=vlsr)

        ms.close()
    pass

def _splitdata( field,
                TM1fulldir,
                vlsr=0, #km/s
                line='SO2_234187',
                freqrange=[], #MHz
                overwrite=False,
               ):
    if not os.path.exists(line):
        os.mkdir(line)
    outdir = os.path.join(line,field)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    elif overwrite:
        shutil.rmtree(outdir) #be careful!
        os.mkdir(outdir)

    #######split out the .ms file##########
    if freqrange[0]>freqrange[1]:
        freqrange = [freqrange[1],freqrange[0]]
    fcenter = (freqrange[0]+freqrange[1])/2.

    splitoutfile = os.path.join(outdir,line+'.ms')
    if not os.path.exists(splitoutfile):
        assert ms.open(TM1fulldir)
        spwinfo=ms.getspectralwindowinfo()
        allspwids = [int(i) for i in spwinfo.keys()]
        allspwids.sort()
        allspwids = ['%i'%i for i in allspwids]
        wins = []
        for j in allspwids:
            namedex = ms.metadata().fieldnames().index(field)  
            fieldid = ms.metadata().fieldsforname()[namedex]
            infreq =  ms.cvelfreqs(spwids=int(j), 
                               fieldids=[fieldid],
                               mode='channel', 
                               width=0, 
                               outframe='LSRK')/1E6
            if (spwinfo[j]['NumChan']>3000) and (np.abs(infreq.mean()-fcenter)<1000):
                spwid = spwinfo[j]['SpectralWindowId']
                index1 = np.argmin(np.abs(infreq-freqrange[0]*(1-vlsr/2.998E5) ))
                index2 = np.argmin(np.abs(infreq-freqrange[1]*(1-vlsr/2.998E5) ))
                print('infreq',index1,index2,infreq[index1], infreq[index2])
                if index1>index2:
                    index1, index2 = index2, index1
                wins.append('%i:%i~%i' %(spwid,index1,index2))
        winstr = ','.join(wins)
        ms.close()
        assert tb.open(TM1fulldir)
        splitrawdatacolumn = 'data'
        if 'CORRECTED_DATA' in tb.colnames():
            splitrawdatacolumn = 'corrected'
        tb.close()
        print("split",TM1fulldir,field,winstr,vlsr)
        split(vis=TM1fulldir,
              outputvis=splitoutfile,
              field=field,
              spw=winstr,
              datacolumn=splitrawdatacolumn)
    ######################################

def cleanline(field,
              line='SO2_234187',
              freqrange=[234157.3,234222.3],
              linefreeflag=[234157.3, 234172.3, 234209.3, 234222.3],
              threshold='20mJy',
              mpiscriptdir = './mpicasarun/',
              ncore=50,
              docubeclean=True):
    print(field,'start')
    vlsr=Ts['vlsrnew'][list(Ts['source']).index(field)]
    datadir = os.path.join(line,field)
    msfile = os.path.join(datadir,line+'.ms')
    contsubmsfile = os.path.join(datadir,line+'_contsub.ms')
    assert os.path.exists(msfile)
    
    if not os.path.exists(contsubmsfile):
        assert ms.open(msfile)
        spwinfo=ms.getspectralwindowinfo()
        allspwids = [int(i) for i in spwinfo.keys()]
        allspwids.sort()
        allspwids = ['%i'%i for i in allspwids]
        linefreewins = []  
        for j in allspwids:
            infreq =  ms.cvelfreqs(spwids=int(j), 
                               mode='channel', 
                               width=0, 
                               outframe='LSRK')/1E6
            spwid = spwinfo[j]['SpectralWindowId']
            indexes = np.array([np.argmin(np.abs(infreq-f*(1-vlsr/2.998E5) )) for f in linefreeflag])
            if indexes[1]<indexes[0]:
                indexes = indexes[::-1]
            linefreewins.append(('%i:' %spwid)+';'.join(['%i~%i'%(i[0],i[1]) for i in indexes.reshape((-1,2))]) )
        linefreewinstr = ','.join(linefreewins)
        uvcontsub(vis=msfile,
                  outputvis=contsubmsfile,
                  fitspec=linefreewinstr,
                  fitorder=1,
                  datacolumn='data',
                 )

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
        scales =[0,5,15,50,150], 
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
    if docubeclean and not os.path.exists(os.path.join(datadir,line+'_contsub_cube.image.pbcor.fits')):
            theimagename = os.path.join(datadir,line+'_contsub_cube')
            fmin = None; fmax = None
            assert ms.open(contsubmsfile)
            spwinfo=ms.getspectralwindowinfo()
            df = None
            for j in spwinfo:
                thef = ms.cvelfreqs(spwids=spwinfo[j]['SpectralWindowId'], mode='channel', 
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
            ms.close()
            df = max(df, 0.488290)
            thewidth = '%fMHz' %(df)
            thenchan = (fmax-fmin)//df
            cleanscript = cleanscripttemp.format(
                vis=os.path.join('../',contsubmsfile),
                calcres = 'True',
                calcpsf = 'True',
                field = field,
                imagename = os.path.join('../',theimagename),
                spw="",
                width=thewidth,
                start='%fMHz' %fmin,
                nchan = '%i' %thenchan,
                threshold = threshold,
                niter='50000000',
                cycleniter='2000',
                nmajor='5',
            ) 

            if not os.path.exists(mpiscriptdir):
                os.mkdir(mpiscriptdir)
            scriptfile = os.path.join(mpiscriptdir,'mpiscript.py')
            with open(scriptfile,'w') as thefile:
                thefile.write(cleanscript)


            mpiscriptrun=os.path.join(mpiscriptdir,'mpiscriptrun.sh')
            with open(mpiscriptrun,'w') as thef:
                scriptfile = 'mpiscript.py'
                thef.write("mpicasa -n {ncore} casa --nogui --nologger -c {script}\n".format(script=scriptfile, ncore=ncore))
            os.system('xfce4-terminal --working-directory="%s" --command "bash mpiscriptrun.sh"' %(os.path.abspath(mpiscriptdir), ) )

           
            impbcor(imagename=os.path.join(datadir,line+'_contsub_cube.image'),
                    pbimage=os.path.join(datadir,line+'_contsub_cube.pb'),
                    outfile=os.path.join(datadir,line+'_contsub_cube.image.pbcor'),
                    overwrite=True,
                   )
            exportfits(imagename=os.path.join(datadir,line+'_contsub_cube.image.pbcor'),
                       fitsimage=os.path.join(datadir,line+'_contsub_cube.image.pbcor.fits'),
                       overwrite=True, 
                      )

    print(field,'conpleted')

def run(fields=[],ncore=50,dofields='onlynopipeline'): #dofields: 'all' 'onlynopipeline' or any others
    if len(fields)==0:
        fields_pipeline = ['I08303-4303', 'I08448-4343_1', 'I08448-4343_2', 'I08470-4243', 'I09002-4732', 'I09018-4816', 'I16524-4300', 
                           'I16547-4247', 'I16562-3959', 'I16571-4029', 'I17006-4215', 'I17136-3617', 'I17143-3700', 'I17158-3901', 
                           'I17160-3707', 'I17439-2845_1', 'I17439-2845_2', 'I17441-2822', 'I17455-2800', 'I17545-2357', 'I17589-2312', 
                           'I17599-2148', 'I18032-2032', 'I18075-2040', 'I18079-1756', 'I18089-1732', 'I18110-1854', 'I18116-1646', 
                           'I18139-1842', 'I18290-0924_1', 'I18290-0924_2', 'I18308-0503_1', 'I18308-0503_2', 'I18311-0809', 'I18411-0338', 
                           'I18434-0242', 'I18440-0148', 'I18445-0222', 'I18461-0113', 'I19078+0901', 'I19095+0930', 'I19097+0847']
        allfields = ['I08303-4303', 'I08448-4343_1', 'I08448-4343_2', 'I08470-4243', 'I09002-4732', 'I09018-4816', 'I09094-4803', 
                     'I16524-4300', 'I16547-4247', 'I16562-3959', 'I16571-4029', 'I17006-4215', 'I17008-4040', 'I17016-4124', 'I17136-3617', 
                     'I17143-3700', 'I17158-3901', 'I17160-3707', 'I17175-3544', 'I17204-3636', 'I17220-3609', 'I17233-3606', 'I17244-3536', 
                     'I17258-3637', 'I17269-3312_1', 'I17269-3312_2', 'I17271-3439_1', 'I17271-3439_2', 'I17278-3541', 'I17439-2845_1', 
                     'I17439-2845_2', 'I17441-2822', 'I17455-2800', 'I17545-2357', 'I17589-2312', 'I17599-2148', 'I18032-2032', 'I18056-1952', 
                     'I18075-2040', 'I18079-1756', 'I18089-1732', 'I18110-1854', 'I18116-1646', 'I18117-1753', 'I18134-1942', 'I18139-1842', 
                     'I18159-1648', 'I18182-1433', 'I18223-1243', 'I18228-1312', 'I18236-1205', 'I18264-1152', 'I18290-0924_1', 'I18290-0924_2', 
                     'I18308-0503_1', 'I18308-0503_2', 'I18311-0809', 'I18314-0720_1', 'I18314-0720_2', 'I18316-0602', 'I18317-0513', 
                     'I18317-0757', 'I18341-0727', 'I18411-0338', 'I18434-0242', 'I18440-0148', 'I18445-0222', 'I18461-0113', 'I18469-0132', 
                     'I18479-0005', 'I18507+0110', 'I18507+0121', 'I18517+0437', 'I18530+0215', 'I19078+0901', 'I19095+0930', 'I19097+0847']
        if dofields == 'all':
            fields = allfields
        if dofields == 'onlynopipeline':
            fields = [i for i in allfields if i not in fields_pipeline]
          
    for field in fields:
        threshold='20mJy'
        if field == 'I17441-2822':
            threshold='30mJy'
        cleanline(field,threshold=threshold,ncore=ncore)
