#!/usr/bin/env python

from desispec.log import get_logger
from desispec.io.brick import Brick

import desibest

from astropy.table import Table, join

import pylab
import argparse
import numpy as np
import sys
import math
from astropy.io import fits
import os.path

def main() :

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--b', type = str, default = None, required=False,
                        help = 'path of DESI brick in b')
    parser.add_argument('--r', type = str, default = None, required=False,
                        help = 'path of DESI brick in r')
    parser.add_argument('--z', type = str, default = None, required=False,
                        help = 'path of DESI brick in z')
    parser.add_argument('--outfile', type = str, default = None, required=False,
                        help = 'path of output file')
### Future support of spectra only in the range [first:nres] 
#    parser.add_argument('--nres', type = int, default = None, required=False,
#                        help = 'max number of results to analyse')
#    parser.add_argument('--first', type = int, default = 0, required=False,
#                        help = 'first result to analyse')
    parser.add_argument('--pathtruth', type = str, default = None, required=False,
                        help = 'path of truth table if does not exist in bricks')
    parser.add_argument('--zbest', type = str, default = None, required=False,
                        help = 'zbest file')
    parser.add_argument('--type', type = str, default = "ELG", required=False,
                        help = 'select a given type')


    args = parser.parse_args()
    log=get_logger()

    log.info("starting")

    if args.b is None:
        args.b = "%s/data/brick-b-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
    if args.r is None:
        args.r = "%s/data/brick-r-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
    if args.z is None:
        args.z = "%s/data/brick-z-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))

    b_brick=Brick(args.b)
    r_brick=Brick(args.r)
    z_brick=Brick(args.z)
    brickname = args.b
    
    
    if args.zbest is None :
        args.zbest = "%s/data/zbest-training-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        log.info("will use default zbest file %s"%args.zbest)
        log.info(" ")
        
    try :
        zb_hdulist=fits.open(args.zbest)
    except :
        log.error("error when parsing %s file:"%args.zbest)
        print sys.exc_info()
        sys.exit(12)
        
    c=3.e5 # light celerity in km/s

#- Checking zbest structure and format

    log.info("Checking zbest structure and format")
    log.info("-----------------------------------")
    if(zb_hdulist[0].size == 0 and zb_hdulist[1].name == 'ZBEST'):
        log.info("zbest file has the required structure")
        log.info(" ")
    else:
        log.error("zbest file does not have the required structure")
        log.error("Check ZBEST structure at http://desidatamodel.readthedocs.org/en/latest/DESI_SPECTRO_REDUX/PRODNAME/bricks/BRICKNAME/zbest-BRICKNAME.html")
        print sys.exc_info()
        sys.exit(12)
    keys=['BRICKNAME', 'TARGETID', 'Z', 'ZERR','ZWARN',  'TYPE', 'SUBTYPE']
    for k in keys:
        try :
            zb_hdulist['ZBEST'].columns.names.index(k)
        except :
            log.error("Missing column %s in %s"%(k,args.zbest))
            log.error("Check ZBEST format at http://desidatamodel.readthedocs.org/en/latest/DESI_SPECTRO_REDUX/PRODNAME/bricks/BRICKNAME/zbest-BRICKNAME.html")
            print sys.exc_info()
            sys.exit(12)
        else:
            log.info("ZBEST hdu has the required columns %s"%k)
    log.info(" ")

#- checking for truth table

    log.info("Checking for truth table")
    log.info("-----------------------------------")
    try:
        b_hdulist = fits.open(args.b)
        truth_table_hdu=b_hdulist['_TRUTH']
        log.info("brick has truth table")
        log.info(" ")
    except KeyError :
        truth_table_hdu=None
    
    if (truth_table_hdu is None):
        try:
            truthfile = fits.open(args.pathtruth)
            truth_table_hdu = truthfile['_TRUTH']
            log.info("Found truth table %s"%args.pathtruth)
            log.info(" ")
        except:
            log.error("A truth table should be provided")
            sys.exit(12)

#- Get results from zbest hdu and infos from truth table
    truth = truth_table_hdu.data
    zbres=zb_hdulist['ZBEST'].data

    zb = zbres['Z']
    zt = truth['TRUEZ']
    zw = zbres['ZWARN']
    
#    for i in range(len(zbres['Z'])):
#        print '%s %s %s'%(zb[i],zt[i],zw[i])
             

#- joining zbest and truth tables  

    zb_zt = join(zbres, truth, keys='TARGETID')

#    print zb_zt['Z','TRUEZ','ZWARN']
    truez=zb_zt['TRUEZ']
    bestz=zb_zt['Z']
    errz=zb_zt['ZERR']
    zwarn = zb_zt['ZWARN']
    n=bestz.size
    if (n == 0):
        log.error("target ids in zbest file are not a subset of target ids in truth table")                                                                                                            
        sys.exit(12) 

#- Select objtype

    log.info("Target type")
    log.info("-----------------------------------")

    obj=dict()
    objtypes=['ELG','LRG','QSO','QSO_BAD','STAR','SKY']
    totobj = len(zbres['TYPE'])

### Future support of spectra only in the range [first:nres] 
#    first = args.first

#    if args.nres is not None :
#        nres = min(args.nres,totobj)
#    else :
#        nres=totobj

#    if (args.first >= totobj -1):
#        log.error("First spectrum to consider can not have index greater than %s"%(totobj-1))
#        first=0

    for o in objtypes:
### Future support of spectra only in the range [first:nres] 
#        index=np.where(truth_table_hdu.data['OBJTYPE'][first:first+nres]== '%s'%o)[0]+first
        first=0
        nres=totobj
        index=np.where(zb_zt['OBJTYPE'][first:first+nres]== '%s'%o)[0]+first
        index=index[:nres]

        obj[o]=len(index)
        log.info("%i %s found"%(obj[o],o))
        if (obj[o] != 0): 
            log.info(" ")
            tz = np.zeros(len(index))
            bz = np.zeros(len(index))
            zw = np.zeros(len(index))
            dv = np.zeros(len(index))

            for i,j in zip(index,range(len(index))):
                tz[j]=truez[i]
                bz[j]=bestz[i]
                zw[j]=zwarn[i]
            dv = c*(bz-tz)/(1+tz)
            dz=dv/c

            true_pos = np.where((np.abs(bz-tz)<0.05) & (zw==0))[0]
            true_neg = np.where((np.abs(bz-tz)>0.05) & (zw!=0))[0]
            false_pos = np.where((np.abs(bz-tz)>0.05) & (zw==0))[0]
            false_neg = np.where((np.abs(bz-tz)<0.05) & (zw!=0))[0]

            #- total
            total = len(true_pos)+len(true_neg)+len(false_pos)+len(false_neg)

            #- computes sample efficiency
            efficiency = float(len(true_pos))/float(total)
            
            #- computes purity
            purity = float(len(true_pos))/float((len(true_pos)+len(false_pos)))

            #- catastrophic failures
            cata_fail = float(len(false_pos))/float(total)

            #- figure of merit
            fom = efficiency*purity

#            verr = np.std(dv[np.where(zw==0)[0]])

            # precision
            zerr = np.std(dz[np.where(zw==0)[0]])

            #accuracy
            acc = np.abs(np.mean(dz[np.where(zw ==0)[0]]))
            
            log.info("=====================================")
            log.info("%s: Precision and accuracy (zwarn=0)"%o)
            log.info("=====================================")
            log.info("zerr: %f, bias: %f"%(zerr,acc))

            if (o == 'ELG'):
                true_oIIflux=zb_zt["OIIFLUX"]
                true_pos_oII = np.where((np.abs(bz-tz)<0.05) & (zw==0) & (true_oIIflux>8e-17))[0]
                true_neg_oII = np.where((np.abs(bz-tz)>0.05) & (zw!=0) & (true_oIIflux>8e-17))[0]
                false_pos_oII = np.where((np.abs(bz-tz)>0.05) & (zw==0) & (true_oIIflux>8e-17))[0]
                false_neg_oII = np.where((np.abs(bz-tz)<0.05) & (zw!=0) & (true_oIIflux>8e-17))[0]

                #- total                                                                              
                total_oII = len(true_pos_oII)+len(true_neg_oII)+len(false_pos_oII)+len(false_neg_oII)
                
                #- computes sample efficiency                                                         
                efficiency_oII = float(len(true_pos_oII))/float(total_oII)
            
                #- computes purity                                                             
                purity_oII = float(len(true_pos_oII))/float((len(true_pos_oII)+len(false_pos)))

                #- catastrophic failures
                cata_fail_oII = float(len(false_pos_oII))/float(total_oII)

                #- figure of merit
                fom_oII = efficiency_oII*purity_oII
               
                log.info("=====================================")
                log.info("%s: For OII > 8e-17 erg/s/cm2"%o)
                log.info("=====================================")
                log.info('Efficiency_oII: %d/%d=%f'%(len(true_pos_oII),total_oII,efficiency_oII))
                log.info('Purity_oII: %d/%d=%f'%(len(true_pos_oII),(len(true_pos_oII)+len(false_pos_oII)),purity_oII))
                log.info('Catastrophic failures_oII: %d/%d=%f'%(len(false_pos_oII),total_oII,cata_fail_oII))
                log.info('FOM_oII: %f x %f=%f'%(efficiency_oII,purity_oII,fom_oII))

            log.info("=====================================")
            log.info("%s: Total sample"%o)
            log.info("=====================================")
            log.info('Efficiency: %d/%d=%f'%(len(true_pos),total,efficiency))
            log.info('Purity: %d/%d=%f'%(len(true_pos),(len(true_pos)+len(false_pos)),purity))
            log.info('Catastrophic failures: %d/%d=%f'%(len(false_pos),total,cata_fail))
            log.info('FOM: %f x %f=%f'%(efficiency,purity,fom))
            log.info("=====================================")

# Requirements for each target class should be added for comparison

            # computes spectrum S/N                                                                                                                                                                     
            mean_ston=np.zeros(len(index))
            mean_ston_oII=np.zeros(len(index))
            for spec in index:
                flux=[b_brick.hdu_list[0].data[spec],r_brick.hdu_list[0].data[spec],z_brick.hdu_list[0].data[spec]]
                ivar=[b_brick.hdu_list[1].data[spec],r_brick.hdu_list[1].data[spec],z_brick.hdu_list[1].data[spec]]
                wave=[b_brick.hdu_list[2].data,r_brick.hdu_list[2].data,z_brick.hdu_list[2].data]
                for i in range(3):
                    mean_ston[spec] += np.sum(np.abs(flux[i])*np.sqrt(ivar[i]))/len(wave[i])

                # computes mean S/N in OII lines for ELG
                if (o == 'ELG'):
                    for i in range(3):
                        ok = np.where((wave[i]>3722) & (wave[i]< 3734))[0]
                        if (len(ok) == 0):
                            break
                        else:
                            mean_ston_oII[spec] += np.sum(np.abs(flux[i][ok])*np.sqrt(ivar[i][ok]))/len(ok)
                
            #- plots
            
            #- histogram z
            
            pylab.figure()
            ok=np.where(zw==0)
            mu = np.mean(dz[ok])
            sigma = np.std(dz[ok])
            n, bins, patches = pylab.hist(dz[ok], 20, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
            gauss = pylab.normpdf(bins, mu, sigma)
            l = pylab.plot(bins, gauss, 'k--', linewidth=1.5)
            pylab.xlabel("Delta z")
            pylab.ylabel("Num. of targest per bin")
            
            pylab.figure()
            nx = 1
            ny = 2
            ai = 1

            ok = np.where(zw==0)
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.set_xlabel("<S/N>")
            a.set_ylabel("Zbest-Ztrue (ZWARN=0)")

            not_ok = np.where(zw !=0)
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.errorbar(mean_ston[not_ok],dz[not_ok],errz[not_ok],fmt="ro")
            a.set_xlabel("<S/N>")
            a.set_ylabel("Zbest-Ztrue (all ZWARN)")

            if (o == 'ELG'):
                pylab.figure()
                nx = 1
                ny = 2
                ai = 1

                ok = np.where(zw==0)
                not_ok = np.where(zw !=0)

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(zb_zt['OIIFLUX'][ok],dz[ok],errz[ok],fmt="bo")
                a.errorbar(zb_zt['OIIFLUX'][not_ok],dz[not_ok],errz[not_ok],fmt="ro")
                a.set_xlabel("True [OII] flux")
                a.set_ylabel("Zbest-Ztrue (all ZWARN)")

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(zb_zt['OIIFLUX'][ok],dz[ok],errz[ok],fmt="bo")
                a.set_xlabel("True [OII] flux")
                a.set_ylabel("Zbest-Ztrue (ZWARN=0)")

            pylab.figure()
            nx=2
            ny=2
            ai=1
    
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(tz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.set_xlabel("Ztrue (all ZWARN)")
            a.set_ylabel("Zbest - Ztrue (all ZWARN)")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(bz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.set_xlabel("Zbest (all ZWARN)")
            a.set_ylabel("Zbest - Ztrue")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.set_xlabel("Ztrue (ZWARN=0)")
            a.set_ylabel("Zbest - Ztrue")

            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.set_xlabel("Zbest (ZWARN=0)")
            a.set_ylabel("Zbest - Ztrue")



            pylab.show()                

if __name__ == '__main__':
    main()
