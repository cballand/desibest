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

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,prog='parsezbest.py', usage='%(prog)s [options]\n\n parsezbest.py without options runs a demo')

    parser.add_argument('--b', type = str, default = None, required=False,
                        help = 'path of DESI brick in b')
    parser.add_argument('--r', type = str, default = None, required=False,
                        help = 'path of DESI brick in r')
    parser.add_argument('--z', type = str, default = None, required=False,
                        help = 'path of DESI brick in z')
    parser.add_argument('--outfile', type = str, default = "parsezbest_results.dat", required=False,
                        help = 'path of output file')
    parser.add_argument('--pathtruth', type = str, default = None, required=False,
                        help = 'path of truth table if does not exist in bricks')
    parser.add_argument('--zbest', type = str, default = None, required=False,
                        help = 'zbest file')

    args = parser.parse_args()
    log=get_logger()

    file=open(args.outfile,"w")                                                                                                                                                                        

    log.info("starting")

    if ((args.zbest is None) and (args.b is None) and (args.r is None) and (args.z is None)):
        args.zbest = "%s/data/zbest-training-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.b = "%s/data/brick-b-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.r = "%s/data/brick-r-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.z = "%s/data/brick-z-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
    elif ((args.zbest is None) or (args.b is None) or (args.r is None) or (args.z is None)):
            log.error("Either all files (b, r and z bricks and zbest) or none should be provided")
            sys.exit(12)

    log.info("Using zbest file %s"%args.zbest)
    log.info("Using %s b brick"%args.b)
    log.info("Using %s r brick"%args.r)
    log.info("Using %s z brick"%args.z)
    log.info(" ")

    file.write("Using zbest file %s\n"%args.zbest)
    file.write("Using %s b brick\n"%args.b)
    file.write("Using %s r brick\n"%args.r)
    file.write("Using %s z brick\n"%args.z)
    file.write("\n")
        

    b_brick=Brick(args.b)
    r_brick=Brick(args.r)
    z_brick=Brick(args.z)
    brickname = args.b
    
    try :
        zb_hdulist=fits.open(args.zbest)
    except :
        log.error("error when parsing %s file:"%args.zbest)
        print sys.exc_info()
        file.close()
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
        file.close()
        print sys.exc_info()
        sys.exit(12)
    keys=['BRICKNAME', 'TARGETID', 'Z', 'ZERR','ZWARN',  'TYPE', 'SUBTYPE']
    for k in keys:
        try :
            zb_hdulist['ZBEST'].columns.names.index(k)
        except :
            log.error("Missing column %s in %s"%(k,args.zbest))
            log.error("Check ZBEST format at http://desidatamodel.readthedocs.org/en/latest/DESI_SPECTRO_REDUX/PRODNAME/bricks/BRICKNAME/zbest-BRICKNAME.html")
            file.close()
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
            file.close()
            sys.exit(12)

#- Get results from zbest hdu and infos from truth table
    truth = truth_table_hdu.data
    zbres=zb_hdulist['ZBEST'].data

    zb = zbres['Z']
    zt = truth['TRUEZ']
    zw = zbres['ZWARN']
    
#- joining zbest and truth tables  

    zb_zt = join(zbres, truth, keys='TARGETID')

    truez=zb_zt['TRUEZ']
    bestz=zb_zt['Z']
    errz=zb_zt['ZERR']
    zwarn = zb_zt['ZWARN']
    n=bestz.size
    if (n == 0):
        log.error("target ids in zbest file are not a subset of target ids in truth table")                                                                                                            
        log.error("did you provide bricks that correspond to zbest file ?")
        file.close()
        sys.exit(12) 

#- Select objtype

    log.info("Target type")
    log.info("-----------------------------------")

    obj=dict()
    objtypes=['ELG','LRG','QSO','QSO_BAD','STAR']
    totobj = len(zb_zt['Z'])

    for o in objtypes:
        index=np.where(zb_zt['OBJTYPE'] == '%s'%o)[0]
        obj[o]=len(index)
        log.info("%i %s found"%(obj[o],o))
        file.write("%i %s found\n"%(obj[o],o))
        if (obj[o] != 0): 
            log.info(" ")
            file.write("\n")
            tz = np.zeros(len(index))
            bz = np.zeros(len(index))
            zw = np.zeros(len(index))
            dv = np.zeros(len(index))
            if (o == 'ELG'): 
                trfloii = np.zeros(len(index))

            for i,j in zip(index,range(len(index))):
                tz[j]=truez[i]
                bz[j]=bestz[i]
                zw[j]=zwarn[i]
                if (o == 'ELG'):
                    trfloii[j] = zb_zt["OIIFLUX"][i]
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
            log.info(" ")

            file.write("=====================================\n")
            file.write("%s: Precision and accuracy (zwarn=0)\n"%o)
            file.write("=====================================\n")
            file.write("zerr: %f, bias: %f\n"%(zerr,acc))
            file.write("\n")

            if (o == 'ELG'):
                true_pos_oII = np.where((np.abs(bz-tz)<0.05) & (zw==0) & (trfloii>8e-17))[0]
                true_neg_oII = np.where((np.abs(bz-tz)>0.05) & (zw!=0) & (trfloii>8e-17))[0]
                false_pos_oII = np.where((np.abs(bz-tz)>0.05) & (zw==0) & (trfloii>8e-17))[0]
                false_neg_oII = np.where((np.abs(bz-tz)<0.05) & (zw!=0) & (trfloii>8e-17))[0]

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
                log.info(" ")

                file.write("=====================================\n")
                file.write("%s: For OII > 8e-17 erg/s/cm2\n"%o)
                file.write("=====================================\n")
                file.write('Efficiency_oII: %d/%d=%f\n'%(len(true_pos_oII),total_oII,efficiency_oII))
                file.write('Purity_oII: %d/%d=%f\n'%(len(true_pos_oII),(len(true_pos_oII)+len(false_pos_oII)),purity_oII))
                file.write('Catastrophic failures_oII: %d/%d=%f\n'%(len(false_pos_oII),total_oII,cata_fail_oII))
                file.write('FOM_oII: %f x %f=%f\n'%(efficiency_oII,purity_oII,fom_oII))
                file.write("\n")


            log.info("=====================================")
            log.info("%s: Total sample"%o)
            log.info("=====================================")
            log.info('Efficiency: %d/%d=%f'%(len(true_pos),total,efficiency))
            log.info('Purity: %d/%d=%f'%(len(true_pos),(len(true_pos)+len(false_pos)),purity))
            log.info('Catastrophic failures: %d/%d=%f'%(len(false_pos),total,cata_fail))
            log.info('FOM: %f x %f=%f'%(efficiency,purity,fom))
            log.info("=====================================")
            log.info(" ")

            file.write("=====================================\n")
            file.write("%s: Total sample\n"%o)
            file.write("=====================================\n")
            file.write('Efficiency: %d/%d=%f\n'%(len(true_pos),total,efficiency))
            file.write('Purity: %d/%d=%f\n'%(len(true_pos),(len(true_pos)+len(false_pos)),purity))
            file.write('Catastrophic failures: %d/%d=%f\n'%(len(false_pos),total,cata_fail))
            file.write('FOM: %f x %f=%f\n'%(efficiency,purity,fom))
            file.write("=====================================\n")
            file.write("\n")


# Requirements for each target class should be added for comparison

            # computes spectrum S/N                                                                                                                                                                     
            mean_ston=np.zeros(totobj)
            mean_ston_oII=np.zeros(totobj)
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
            
            #- histograms
            
            pylab.figure()
            ok=np.where(zw==0)
            mu = np.mean(dz[ok])
            sigma = np.std(dz[ok])
            n, bins, patches = pylab.hist(dz[ok], 20, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
            gauss = pylab.normpdf(bins, mu, sigma)
            l = pylab.plot(bins, gauss, 'k--', linewidth=1.5, label="%s"%o)
            pylab.xlabel("(zb-zt)/(1+zt)")
            pylab.ylabel("Num. of %s targest per bin"%o)
            

            pylab.figure()
            ok=np.where(zw==0)
            mu = np.mean(dv[ok])
            sigma = np.std(dv[ok])
            n, bins, patches = pylab.hist(dv[ok], 20, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
            gauss = pylab.normpdf(bins, mu, sigma)
            l = pylab.plot(bins, gauss, 'k--', linewidth=1.5, label="%s"%o)
            pylab.xlabel("Delta v = c(zb-zt)/(1+zt) [km/s]")
            pylab.ylabel("Num. of %s targest per bin"%o)


            pylab.figure()
            nx = 1
            ny = 2
            ai = 1

            ok = np.where(zw==0)
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.set_xlabel("%s <S/N>"%o)
            a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0)")

            not_ok = np.where(zw !=0)
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.errorbar(mean_ston[not_ok],dz[not_ok],errz[not_ok],fmt="ro")
            a.set_xlabel("%s <S/N> "%o)
            a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")

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
                a.set_xlabel("%s True [OII] flux"%o)
                a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(zb_zt['OIIFLUX'][ok],dz[ok],errz[ok],fmt="bo")
                a.set_xlabel("%s True [OII] flux"%o)
                a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0)")

            pylab.figure()
            nx=2
            ny=2
            ai=1
    
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(tz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.set_xlabel("%s zt (all ZWARN)"%o)
            a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(bz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.set_xlabel("%s zb (all ZWARN)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.set_xlabel("%s zt (ZWARN=0)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")

            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.set_xlabel("%s zb (ZWARN=0)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")



            pylab.show()                


if __name__ == '__main__':
    main()
