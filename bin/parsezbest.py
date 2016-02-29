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

def desi_requirement(obj):
    req=dict()
    objtype=['ELG','LRG','QSO']
    if obj not in objtype:
        return
    
    if (obj == 'ELG'):
        req['SIG_Z']=0.0005
        req['SIG_V']=150.
        req['BIAS_Z']=0.0002
        req['BIAS_V']=60.
        req['CATA_FAIL_THRES']=1000.
        req['CATA_FAIL_MAX']=0.05
        req['EFFICIENCY']=0.9

    if (obj == 'LRG'):
        req['SIG_Z']=0.0005
        req['SIG_V']=150.
        req['BIAS_Z']=0.0002
        req['BIAS_V']=60.
        req['CATA_FAIL_THRES']=1000.
        req['CATA_FAIL_MAX']=0.05
        req['EFFICIENCY']=0.95
    
    if (obj == 'QSO'):
        req['SIG_Z']=0.0025
        req['SIG_V']=750.
        req['BIAS_Z']=0.0004
        req['BIAS_V']=120.
        req['CATA_FAIL_THRES']=1000.
        req['CATA_FAIL_MAX']=0.05
        req['EFFICIENCY']=0.9

    return req

def main() :
    """
    parsezbest.py computes common metrics and makes plots for analyzing results of zdc1 redshift challenge (zbest file).

    Metrics:
    + dz = (zbest-ztrue)/(1+ztrue)
    + dv = c*dz
    + pull = (zbest-ztrue)/zerr
    + precision: sigma_z = std(dz), sigma_v = std(dv), nmad_z, nmad_v
    + accuracy (bias): mu_z = mean(dz), mu_v = mean(dv)
    + efficiency
    + purity
    + % of catastrophic failures
    + FOM = purity*efficiency

    Plots:
    + Histograms dz, dv, pull
    + dz as a function of zt and zb for:
        - zwarn=0
        - zwarn=0 without catastrophic failures
        - zwarn!=0
    + dz as a function of average S/N per wavelegnth bin for:
        - zwarn=0                                                                                                                                                                                     
        - zwarn=0 without catastrophic failures
        - zwarn!=0      
    
    Color code:
    + zwarn=0: blue filled circles
    + zwarn !=0: red filled circles
    + catastrophic failures: green filled circles

    """
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

    #- if no arguments is passed, parsezbest runs a demo
    if ((args.zbest is None) and (args.b is None) and (args.r is None) and (args.z is None)):
        args.zbest = "%s/data/zbest-training-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.b = "%s/data/brick-b-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.r = "%s/data/brick-r-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
        args.z = "%s/data/brick-z-elg-100-zztop.fits"%(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(desibest.__file__)))))
    elif ((args.zbest is None) or (args.b is None) or (args.r is None) or (args.z is None)):
            log.error("Either all files (b, r and z bricks and zbest) or none should be provided")
            sys.exit(12)

    try:
        b_brick=Brick(args.b)
    except:
        log.error("can not open brick %s"%args.b)
        print sys.exc_info()
        file.close()
        sys.exit(12)
    try:
        r_brick=Brick(args.r)
    except:
        log.error("can not open brick %s"%args.r)
        print sys.exc_info()
        file.close()
        sys.exit(12)
    try:
        z_brick=Brick(args.z)
    except:
        log.error("can not open brick %s"%args.z)
        print sys.exc_info()
        file.close()
        sys.exit(12)

    try :
        zb_hdulist=fits.open(args.zbest)
    except :
        log.error("can not open file %s:"%args.zbest)
        print sys.exc_info()
        file.close()
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

    brickname = args.b
    
    c=3.e5 # light celerity in vacuum [km/s]

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
        log.error("did you provide the bricks that correspond to zbest file ?")
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

            #load requirements for target
            req= desi_requirement(o)
            
            tz = np.zeros(len(index))
            bz = np.zeros(len(index))
            zw = np.zeros(len(index))
            dv = np.zeros(len(index))
            ez = np.zeros(len(index))
            if (o == 'ELG'): 
                trfloii = np.zeros(len(index))

            dz=0.
            dv=0.
            pull=0.

            for i,j in zip(index,range(len(index))):
                tz[j]=truez[i]
                bz[j]=bestz[i]
                zw[j]=zwarn[i]
                ez[j]=errz[i]
                if (o == 'ELG'):
                    trfloii[j] = zb_zt["OIIFLUX"][i]
            dv = c*(bz-tz)/(1+tz)
            dz=dv/c
            pull =(bz-tz)/ez

            true_pos = np.where((np.abs(dz)<0.0033) & (zw==0))[0]
            true_neg = np.where((np.abs(dz)>0.0033) & (zw!=0))[0]
            false_pos = np.where((np.abs(dz)>0.0033) & (zw==0))[0]
            false_neg = np.where((np.abs(dz)<0.0033) & (zw!=0))[0]

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

            # precision
            ok = np.where(zw==0)[0]
            zerr = np.std(dz[ok])
            verr = np.std(dv[ok])

            #accuracy
            zacc = np.abs(np.mean(dz[ok]))
            vacc = np.abs(np.mean(dv[ok]))

            ok_no_cata = np.where((zw==0) & (np.abs(dz)<0.0033))[0]
            if (len(ok) != len(ok_no_cata)):
                zerr_no_cata = np.std(dz[ok_no_cata])
                verr_no_cata = np.std(dv[ok_no_cata])
                zacc_no_cata = np.abs(np.mean(dz[ok_no_cata]))
                vacc_no_cata = np.abs(np.mean(dv[ok_no_cata]))

            # NMAD
            nmad_z = np.median(np.abs(dz[ok]-np.median(dz[ok])))
            nmad_z *= 1.4826
            nmad_v = np.median(np.abs(dv[ok]-np.median(dv[ok])))
            nmad_v *= 1.4826

            #pull
            mu_pull = np.mean(pull[ok_no_cata])
            sigma_pull = np.std(pull[ok_no_cata])

            # zwarn
            zw0 = len(np.where(zw == 0)[0])
            zw_non0 = len(np.where(zw != 0)[0])
#            assert (zw0 + zw_non0) != len(index), "zw=0 + zw!=0 not equal to number of objects"

            log.info("=====================================")
            log.info("%s: Precision and accuracy (zwarn=0)"%o)
            log.info("=====================================")
            log.info("sigma_z: %f, mu_z: %f"%(zerr,zacc))
            log.info("NMAD_z: %f"%nmad_z)
            log.info("sigma_v: %f, mu_v: %f"%(verr,vacc))
            log.info("NMAD_v: %f"%nmad_v)
            if req is not None:
                if (zerr>req['SIG_Z']):
                    log.info("sigma_z & sigma_v do not meet DESI requirements on precision for %s"%o)
                if (zacc>req['BIAS_Z']):
                    log.info("mu_z & mu_v do not meet DESI requirements on bias for %s"%o)
            log.info(" ")
            if (len(ok) != len(ok_no_cata)):
                log.info("=====================================")
                log.info("%s: Precision and accuracy "%o)
                log.info("zwarn=0 without catastrophic failures")
                log.info("=====================================")
                log.info("sigma_z: %f, mu_z: %f"%(zerr_no_cata,zacc_no_cata))
                log.info("sigma_v: %f, mu_v: %f"%(verr_no_cata,vacc_no_cata))
                if req is not None:
                    if (zerr_no_cata>req['SIG_Z']):
                        log.info("sigma_z & sigma_v do not meet DESI requirements on precision for %s"%o)
                    if (zacc_no_cata>req['BIAS_Z']):
                        log.info("mu_z & mu_v do not meet DESI requirements on bias for %s"%o)
                log.info(" ")

            file.write("=====================================\n")
            file.write("%s: Precision and accuracy (zwarn=0)\n"%o)
            file.write("=====================================\n")
            file.write("sigma_z: %f, mu_z: %f\n"%(zerr,zacc))
            file.write("NMAD_z: %f\n"%nmad_z)
            file.write("sigma_v: %f, mu_v: %f\n"%(verr,vacc))
            file.write("NMAD_v: %f\n"%nmad_v)
            if req is not None:
                if (zerr>req['SIG_Z']):
                    file.write("sigma_z & sigma_v do not meet DESI requirements on precision for %s\n"%o)
                if (zacc>req['BIAS_Z']):
                    file.write("mu_z & mu_v do not meet DESI requirements on bias for %s\n"%o)
            file.write("\n")
            if (len(ok) != len(ok_no_cata)):
                file.write("=====================================\n")
                file.write("%s: Precision and accuracy \n"%o)
                file.write("zwarn=0 without catastrophic failures\n")
                file.write("=====================================\n")
                file.write("sigma_z: %f, mu_z: %f\n"%(zerr_no_cata,zacc_no_cata))
                file.write("sigma_v: %f, mu_v: %f\n"%(verr_no_cata,vacc_no_cata))
                if req is not None:
                    if (zerr_no_cata>req['SIG_Z']):
                        file.write("sigma_z & sigma_v do not meet DESI requirements on precision for %s\n"%o)
                    if (zacc_no_cata>req['BIAS_Z']):
                        file.write("mu_z & mu_v do not meet DESI requirements on bias for %s\n"%o)
                file.write("\n")



            log.info("=====================================")
            log.info("%s: Pull (zwarn=0, no cata. fail.)"%o)
            log.info("=====================================")
            log.info("mu: %f, sigma: %f"%(mu_pull,sigma_pull))
            log.info(" ")

            file.write("=====================================\n")
            file.write("%s: Pull (zwarn=0, no cata. fail.)\n"%o)
            file.write("=====================================\n")
            file.write("mu: %f, sigma: %f\n"%(mu_pull,sigma_pull))
            file.write("\n")

            if (o == 'ELG'):
                true_pos_oII = np.where((np.abs(dz)<0.0033) & (zw==0) & (trfloii>8e-17))[0]
                true_neg_oII = np.where((np.abs(dz)>0.0033) & (zw!=0) & (trfloii>8e-17))[0]
                false_pos_oII = np.where((np.abs(dz)>0.0033) & (zw==0) & (trfloii>8e-17))[0]
                false_neg_oII = np.where((np.abs(dz)<0.0033) & (zw!=0) & (trfloii>8e-17))[0]

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
                if req is not None:
                    if (efficiency_oII < req['EFFICIENCY']):
                        log.info("Efficiency_oII does not meet DESI requirements for %s"%o)
                log.info('Purity_oII: %d/%d=%f'%(len(true_pos_oII),(len(true_pos_oII)+len(false_pos_oII)),purity_oII))
                log.info('Catastrophic failures_oII: %d/%d=%f'%(len(false_pos_oII),total_oII,cata_fail_oII))
                if req is not None:
                    if (cata_fail_oII>req['CATA_FAIL_MAX']):
                        log.info("Catastrophic failure rate does not meet DESI requirements for %s"%o)
                log.info('FOM_oII: %f x %f=%f'%(efficiency_oII,purity_oII,fom_oII))
                log.info(" ")

                file.write("=====================================\n")
                file.write("%s: For OII > 8e-17 erg/s/cm2\n"%o)
                file.write("=====================================\n")
                file.write('Efficiency_oII: %d/%d=%f\n'%(len(true_pos_oII),total_oII,efficiency_oII))
                if req is not None:
                    if (efficiency_oII < req['EFFICIENCY']):
                        file.write("Efficiency_oII does not meet DESI requirements for %s"%o)
                file.write('Purity_oII: %d/%d=%f\n'%(len(true_pos_oII),(len(true_pos_oII)+len(false_pos_oII)),purity_oII))
                file.write('Catastrophic failures_oII: %d/%d=%f\n'%(len(false_pos_oII),total_oII,cata_fail_oII))
                if req is not None:
                    if (cata_fail_oII>req['CATA_FAIL_MAX']):
                        file.write("Catastrophic failure rate does not meet DESI requirements for %s"%o)
                file.write('FOM_oII: %f x %f=%f\n'%(efficiency_oII,purity_oII,fom_oII))
                file.write("\n")


            log.info("=====================================")
            log.info("%s: Total sample"%o)
            log.info("=====================================")
            log.info("zwarn = 0: %d"%zw0)
            log.info("zwarn !=0: %d"%zw_non0)
            log.info('Efficiency: %d/%d=%f'%(len(true_pos),total,efficiency))
            if req is not None:
                if (efficiency < req['EFFICIENCY']):
                    log.info("Efficiency does not meet DESI requirements for %s"%o)
            log.info('Purity: %d/%d=%f'%(len(true_pos),(len(true_pos)+len(false_pos)),purity))
            log.info('Catastrophic failures: %d/%d=%f'%(len(false_pos),total,cata_fail))
            if req is not None:
                if (cata_fail>req['CATA_FAIL_MAX']):
                    log.info("Catastrophic failure rate does not meet DESI requirements for %s"%o)
            log.info('FOM: %f x %f=%f'%(efficiency,purity,fom))
            log.info("=====================================")
            log.info(" ")

            file.write("=====================================\n")
            file.write("%s: Total sample\n"%o)
            file.write("=====================================\n")
            file.write('Efficiency: %d/%d=%f\n'%(len(true_pos),total,efficiency))
            if req is not None:
                if (efficiency < req['EFFICIENCY']):
                    file.write("Efficiency does not meet DESI requirements for %s\n"%o)
            file.write('Purity: %d/%d=%f\n'%(len(true_pos),(len(true_pos)+len(false_pos)),purity))
            file.write('Catastrophic failures: %d/%d=%f\n'%(len(false_pos),total,cata_fail))
            if req is not None:
                if (cata_fail>req['CATA_FAIL_MAX']):
                    file.write("Catastrophic failure rate does not meet DESI requirements for %s\n"%o)
            file.write('FOM: %f x %f=%f\n'%(efficiency,purity,fom))
            file.write("=====================================\n")
            file.write("\n")


            # computes spectrum S/N                                                                                                                                                                     
            mean_ston=np.zeros(len(index))
            mean_ston_oII=np.zeros(len(index))
            for spec,sp in zip(index,range(len(index))):
                flux=[b_brick.hdu_list[0].data[spec],r_brick.hdu_list[0].data[spec],z_brick.hdu_list[0].data[spec]]
                ivar=[b_brick.hdu_list[1].data[spec],r_brick.hdu_list[1].data[spec],z_brick.hdu_list[1].data[spec]]
                wave=[b_brick.hdu_list[2].data,r_brick.hdu_list[2].data,z_brick.hdu_list[2].data]
                for i in range(3):
                    mean_ston[sp] += np.sum(np.abs(flux[i])*np.sqrt(ivar[i]))/len(wave[i])

                # computes mean S/N in OII lines for ELG
                if (o == 'ELG'):
                    for i in range(3):
                        ok = np.where((wave[i]>3722) & (wave[i]< 3734))[0]
                        if (len(ok) == 0):
                            break
                        else:
                            mean_ston_oII[sp] += np.sum(np.abs(flux[i][ok])*np.sqrt(ivar[i][ok]))/len(ok)
                
            #- plots

            ok=np.where(zw==0)[0]
            cata = np.where((zw == 0) & (np.abs(dz)>0.0033))[0]
            not_ok = np.where(zw !=0)[0]
#            ok_no_cata = np.where((zw == 0) & (np.abs(dz)<0.0033))[0]

            #- histograms
            
            pylab.figure()
            n, bins, patches = pylab.hist(dz[ok_no_cata], 30, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.75)
            if (o != 'QSO'):
                muz = np.mean(dz[ok_no_cata])
                sigmaz = np.std(dz[ok_no_cata])
                gauss = pylab.normpdf(bins, muz, sigmaz)
                l = pylab.plot(bins, gauss, 'k--', linewidth=1.5, label="mu=%2.0f *1e-6, sig=%2.0f *1e-6"%(muz/1e-6,sigmaz/1e-6))
                pylab.legend()
            pylab.xlabel("(zb-zt)/(1+zt) (ZWARN=0 without catastrophic failures)")
            pylab.ylabel("Num. of %s targets per bin"%o)
            

            pylab.figure()
            n, bins, patches = pylab.hist(dv[ok_no_cata], 30, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
            if (o != 'QSO'): 
                muv = np.mean(dv[ok_no_cata])
                sigmav = np.std(dv[ok_no_cata])
                gauss = pylab.normpdf(bins, muv, sigmav)
                l = pylab.plot(bins, gauss, 'k--', linewidth=1.5, label="mu=%2.0f, sig=%2.0f"%(muv,sigmav))
                pylab.legend()
            pylab.xlabel("Delta v = c(zb-zt)/(1+zt) [km/s] (ZWARN=0 without catastrophic failures)")
            pylab.ylabel("Num. of %s targets per bin"%o)

            #- pull distribution

            pylab.figure()
            n, bins, patches = pylab.hist(pull[ok_no_cata], 30, normed=1, histtype='stepfilled')
            pylab.setp(patches, 'facecolor', 'c', 'alpha', 0.75)
            mu_pull = np.mean(pull[ok_no_cata])
            sigma_pull = np.std(pull[ok_no_cata])
            gauss = pylab.normpdf(bins, mu_pull, sigma_pull)
            l = pylab.plot(bins, gauss, 'k--', linewidth=1.5, label="mu=%2.3f, sig=%2.3f"%(mu_pull,sigma_pull))
            pylab.legend()
            mu=0.
            sig=1.
            gauss1 = pylab.normpdf(bins, mu, sig)
            l1 = pylab.plot(bins, gauss1, 'k--', linewidth=1.5, color='r', label="mu=0., sig=1.")
            pylab.legend()
            pylab.xlabel("Pull = (zb-<zt>)/zerr (ZWARN=0 without catastrophic failures)")
            pylab.ylabel("Num. of %s targets per bin"%o)

            #- other plots

            pylab.figure()
            nx = 1
            if len(cata) !=0:
                ny = 3
            else:
                ny = 2
            ai = 1

            # catastrophic failures in green
            # zwarn != 0 in red
            # zw =0 no catastrophic in blue

            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.errorbar(mean_ston[cata],dz[cata],errz[cata],fmt="go")
            a.set_xlabel("%s <S/N>"%o)
            a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0)")

            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(mean_ston[ok],dz[ok],errz[ok],fmt="bo")
            a.errorbar(mean_ston[not_ok],dz[not_ok],errz[not_ok],fmt="ro")
            a.errorbar(mean_ston[cata],dz[cata],errz[cata],fmt="go")
            a.set_xlabel("%s <S/N> "%o)
            a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")

            if len(cata) !=0:
                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(mean_ston[ok_no_cata],dz[ok_no_cata],errz[ok_no_cata],fmt="bo")
                a.set_xlabel("%s <S/N> "%o)
                a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0, no cata. fail.)")


            if (o == 'ELG'):
                pylab.figure()
                nx=1
                if len(cata) !=0:
                    ny=3
                else:
                    ny=2
                ai=1

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(zb_zt['OIIFLUX'][ok],dz[ok],errz[ok],fmt="bo")
                a.errorbar(zb_zt['OIIFLUX'][cata],dz[cata],errz[cata],fmt="ro")
                a.set_xlabel("%s True [OII] flux"%o)
                a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(zb_zt['OIIFLUX'][ok],dz[ok],errz[ok],fmt="bo")
                a.errorbar(zb_zt['OIIFLUX'][not_ok],dz[not_ok],errz[not_ok],fmt="ro")
                a.errorbar(zb_zt['OIIFLUX'][cata],dz[cata],errz[cata],fmt="go")
                a.set_xlabel("%s True [OII] flux"%o)
                a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0)")

                if len(cata) != 0:
                    a=pylab.subplot(ny,nx,ai); ai +=1
                    a.errorbar(zb_zt['OIIFLUX'][ok_no_cata],dz[ok_no_cata],errz[ok_no_cata],fmt="bo")
                    a.set_xlabel("%s True [OII] flux"%o)
                    a.set_ylabel("(zb-zt)/(1+zt) (ZWARN=0, no cata. fail.)")

            pylab.figure()
            nx=2
            if len(cata) !=0:
                ny=3
            else:
                ny=2
            ai=1
    
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(tz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.errorbar(tz[cata],dz[cata],errz[cata],fmt="o",c="g")
            a.set_xlabel("%s zt (all ZWARN)"%o)
            a.set_ylabel("(zb-zt)/(1+zt) (all ZWARN)")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(bz[not_ok],dz[not_ok],errz[not_ok],fmt="o",c="r")
            a.errorbar(bz[cata],dz[cata],errz[cata],fmt="o",c="g")
            a.set_xlabel("%s zb (all ZWARN)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")
            
            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(tz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(tz[cata],dz[cata],errz[cata],fmt="o",c="g")
            a.set_xlabel("%s zt (ZWARN=0)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")

            a=pylab.subplot(ny,nx,ai); ai +=1
            a.errorbar(bz[ok],dz[ok],errz[ok],fmt="o",c="b")
            a.errorbar(bz[cata],dz[cata],errz[cata],fmt="o",c="g")
            a.set_xlabel("%s zb (ZWARN=0)"%o)
            a.set_ylabel("(zb-zt)/(1+zt)")

            if len(cata) !=0:
                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(tz[ok_no_cata],dz[ok_no_cata],errz[ok_no_cata],fmt="o",c="b")
                a.set_xlabel("%s zt (ZWARN=0, no cata. fail.)"%o)
                a.set_ylabel("(zb-zt)/(1+zt)")

                a=pylab.subplot(ny,nx,ai); ai +=1
                a.errorbar(bz[ok_no_cata],dz[ok_no_cata],errz[ok_no_cata],fmt="o",c="b")
                a.set_xlabel("%s zb (All ZWARN, no cata. fail.)"%o)
                a.set_ylabel("(zb-zt)/(1+zt)")


            pylab.show()                


if __name__ == '__main__':
    main()
