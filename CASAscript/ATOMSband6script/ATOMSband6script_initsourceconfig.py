import astropy.table as t
import numpy as np
import os
from casatools import ms
ms = ms() 

#function to get the spwid from the original ms files 
def getspwids_for_inmsfile(inmsfile,sourceconfig):
    selectspwids = []
    for i in sourceconfig['wideband_outids']:
        selectspwids.extend(sourceconfig['widebands'][i][inmsfile])  
    return ','.join(selectspwids) 

def get_sourceconfig(field,infocsv='ATMOSband6_sourceinfo.csv'):
    sourceconfig = {
    'field':'',
    'inmsfiles': [], #['TM1','TM2','ACA'],
    'inmsfilespath':{
        #'TM1':'../2021.1.00095.S/science_goal.uid___A001_X1590_X33cf/group.uid___A001_X1590_X341d/member.uid___A001_X1590_X341e/calibrated/uid___A002_Xfcf230_X2e1.ms',
        #'TM2':'../2021.1.00095.S/science_goal.uid___A001_X1590_X33cf/group.uid___A001_X1590_X341d/member.uid___A001_X1590_X3420/calibrated/uid___A002_Xf7519f_X1a70.ms',
        #'ACA':'../2021.1.00095.S/science_goal.uid___A001_X1590_X33cf/group.uid___A001_X1590_X341d/member.uid___A001_X1590_X3422/calibrated/uid___A002_X10275c0_X12202.ms',
    },
    'wideband_outids':['0','1','2','3'],
    'widebandsfreq':{
       '0':217900,
       '1':220296,
       '2':231340,       
       '3':233480,
    },
    'widebands':{
       '0':{},
       '1':{},
       '2':{},
       '3':{},
    },
    'inmsfilespwids':{},
    'widebandsinconcatfile':{
       '0':[],
       '1':[],
       '2':[],
       '3':[],
    },
    'continuumthreshold':2, #mJy
    'selfcal': True,
    'selfcal_threshold':[15,8], #mJy, first is to cali phase, the second is to cali amp&phase.
    'vlsr':71,
    'doselfcal':True,
    'phasethreshold':'25mJy',
    'ampthreshold':'10mJy',
    'contthreshold':'2mJy',
    'linethreshold':'30mJy',
    }
    sourceconfig['field'] = field

    T=t.Table.read(infocsv,comment=r'^#.*')
    assert field in T['field'] 
    index = np.where(T['field']==field)[0][0]
    scienceid = T['scienceid'][index]
    groupid = T['groupid'][index]
    for inmsfile in ['TM1','TM2','ACA']:
        if T[inmsfile][index].upper() != 'NONE':
            sourceconfig['inmsfiles'].append(inmsfile)
            filepath = os.path.join('../2021.1.00095.S',
                                    scienceid,
                                    groupid,
                                    T[inmsfile][index]
                                   )
            sourceconfig['inmsfilespath'][inmsfile] = filepath
    sourceconfig['vlsr'] = T['vlsr'][index]
    sourceconfig['doselfcal'] = T['doselfcal'][index]>0
    sourceconfig['phasethreshold'] = '%fmJy' %( T['phasethreshold'][index] )
    sourceconfig['ampthreshold'] = '%fmJy' %( T['ampthreshold'][index])
    sourceconfig['contthreshold'] = '%fmJy' %( T['contthreshold'][index])
    sourceconfig['linethreshold'] = '%fmJy' %( T['linethreshold'][index])
    sourceconfig['linewinmethod'] = T['linewinmethod'][index]
    sourceconfig['done'] = T['done'][index]>0
    sourceconfig['scienceid'] = scienceid
    sourceconfig['groupid'] = groupid

    #get the spwid from the original ms files
    for inmsfile in sourceconfig['inmsfiles']:
        inmsfilepath = sourceconfig['inmsfilespath'][inmsfile]
        assert ms.open(inmsfilepath)
        spwinfo=ms.getspectralwindowinfo()
        allspwids = spwinfo.keys()
        for i in sourceconfig['wideband_outids']:
            outfreq = sourceconfig['widebandsfreq'][i]
            for j in allspwids:
                if spwinfo[j]['NumChan']>3000: 
                    infreq =  spwinfo[j]['Chan1Freq'] + spwinfo[j]['ChanWidth']*spwinfo[j][ 'NumChan' ]/2.
                    infreq = infreq/1E6
                    if np.abs(infreq-outfreq)<1000:
                        if inmsfile not in sourceconfig['widebands'][i]:
                            sourceconfig['widebands'][i][inmsfile] = []
                        sourceconfig['widebands'][i][inmsfile].append( '%i' % (spwinfo[j]['SpectralWindowId']) )

    for inmsfile in sourceconfig['inmsfiles']:
        inmsfilepath = sourceconfig['inmsfilespath'][inmsfile]
        sourceconfig['inmsfilespwids'][inmsfile] = getspwids_for_inmsfile(inmsfile,sourceconfig)

    return sourceconfig

def getwideband(f, sourceconfig):
    select = False
    mindf = 0
    chose = '0'
    for i in sourceconfig['wideband_outids']:
        f0 = sourceconfig['widebandsfreq'][i]
        df = np.abs(f-f0)
        if not select:
            mindf = df
            chose = i
            select = True
        else:
            if df<mindf:
                mindf = df
                chose = i
    return chose
