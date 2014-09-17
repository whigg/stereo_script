#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 21:52:23 2013

@author: ben

dem_cull.py
parameters:
    input_file, output_file, smooth_scale, smooth_tol, slope_tol, R_tol, dliate_by
"""

import argparse
from osgeo import gdal, gdalconst
from im_subset import im_subset
import numpy as np
import scipy.ndimage as snd
import matplotlib.pyplot as plt
import sys, os, time
blocksize=4096
DOPLOT=None


def smooth_corrected(z, mask, w_smooth):
     mask1=snd.gaussian_filter(np.float32(mask), w_smooth, mode="constant", cval=0)
     ztemp=np.nan_to_num(z)
     ztemp[mask==0]=0.0
     zs=snd.gaussian_filter(ztemp, w_smooth, mode="constant", cval=0)
     zs[mask1>0]=zs[mask1>0]/mask1[mask1>0]
     return zs, mask1

parser = argparse.ArgumentParser(description='cull out spurious values from a DEM')
parser.add_argument('input_file')
parser.add_argument('output_file')
parser.add_argument('--smooth_scale','-s', type=float, default=16.)
parser.add_argument('--smooth_tol','-t', type=float, default=2.)
parser.add_argument('--slope_tol', '-m', type=float, default=10*np.pi/180.)
parser.add_argument('--R_tol','-r', type=float, default=5.)
parser.add_argument('--erode_By','-b', type=float, default=1.)
parser.add_argument('--decimate_by','-d', type=float, default=1.)
parser.add_argument('--error_RMS_scale','-e', type=float, default=64)
parser.add_argument('--geolocation_error','-g', type=float, default=5)
parser.add_argument('--facet_tol', '-f', type=float, default=None)
args=parser.parse_args()

ds=gdal.Open(args.input_file);
driver = ds.GetDriver()
band=ds.GetRasterBand(1)
noData=band.GetNoDataValue()
if noData is None:
    noData=0.

dec=args.decimate_by

nX=band.XSize;
nY=band.YSize;
xform_in=np.array(ds.GetGeoTransform())
dx=xform_in[1]

xform_out=xform_in.copy()
xform_out[1]=xform_in[1]*dec
xform_out[5]=xform_in[5]*dec
if np.mod(dec,2)==0:   # shift output origin by 1/2 pixel if dec is even
    xform_out[0]=xform_in[0]+xform_in[1]/2
    xform_out[3]=xform_in[3]+xform_in[5]/2

nX_out=np.int(nX/dec)
nY_out=np.int(nY/dec)
if args.error_RMS_scale is not None:
    out_bands=[1,2]
    w_error=np.int(args.error_RMS_scale/dx/dec)
else:
    out_bands=[1]

if os.path.isfile(args.output_file):
    print "output_file %s exists, deleting" % args.output_file
    os.remove(args.output_file)

outDs = driver.Create(args.output_file, nX_out, nY_out, len(out_bands), gdalconst.GDT_Float32)

argDict=vars(args)
for key in argDict:
    if argDict[key] is not None:
        print "\t%s is %s" %(key, str(argDict[key]))
        outDs.SetMetadataItem("dem_filter_"+key, str(argDict[key]))

if args.smooth_scale is not None:
    w_smooth=args.smooth_scale/dx

if args.erode_By is not None:
    N_erode=np.ceil(args.erode_By/dx)
    xg,yg=np.meshgrid(np.arange(0, N_erode)-N_erode/2, np.arange(0, N_erode)-N_erode/2)
    k_erode=(xg**2 + yg**2) <= N_erode/2.

if args.facet_tol is not None:
    xxg, yyg=np.meshgrid(np.arange(-8., 9), np.arange(-8., 9.))
    opening_kernel=(xxg**2+yyg**2 <= 25)
    closing_kernel=(xxg**2+yyg**2 <= 64)

pad=np.max([1, np.int((w_smooth+N_erode)/dec)]);

stride=np.int(blocksize/dec)
in_sub=im_subset(0, 0, nX, nY, ds, pad_val=0, Bands=[1])

last_time=time.time()

for out_sub in im_subset(0, 0,  nX_out,  nY_out, outDs, pad_val=0, Bands=out_bands, stride=stride, pad=pad):
    in_sub.setBounds(out_sub.c0*dec, out_sub.r0*dec, out_sub.Nc*dec, out_sub.Nr*dec, update=True)
    z=in_sub.z[0,:,:]
    mask=np.ones_like(in_sub.z[0,:,:])
    mask[np.isnan(in_sub.z[0,:,:])]=0
    mask[in_sub.z[0,:,:]==noData]=0
    out_temp=np.zeros([len(out_bands), stride, stride])

    if np.all(mask.ravel()==0):
        out_temp=out_temp+np.NaN
        out_sub.z=out_temp
        out_sub.setBounds(out_sub.c0+pad, out_sub.r0+pad, out_sub.Nc-2*pad, out_sub.Nr-2*pad)
        out_sub.writeSubsetTo(out_bands, out_sub)
        continue

    if (args.R_tol is not None) | (args.facet_tol is not None):
        lap=np.abs(snd.laplace(in_sub.z[0,:,:], mode='constant', cval=0.0))

    if args.R_tol is not None:
        mask[lap>args.R_tol]=0

    if args.facet_tol is not None:
        mask1=mask.copy()
        mask1[lap < args.facet_tol]=0
        mask1=snd.binary_closing(snd.binary_opening(mask1, structure=opening_kernel), structure=closing_kernel)
        #mask1=snd.binary_erosion(mask1, structure=simplify_kernel);
        mask[mask1==0]=0

    if args.smooth_scale is not None:
        zs, mask2 = smooth_corrected(z, mask, w_smooth)
        mask[np.abs(in_sub.z[0,:,:]-zs)>args.smooth_tol]=0.

    if args.slope_tol is not None:
        gx, gy=np.gradient(zs, dx, dx)
        mask[gx**2+gy**2 > args.slope_tol**2]=0
        z[mask==0]=0

    if args.erode_By is not None:
        mask=snd.binary_erosion(mask, k_erode)
        z[mask==0]=0

    if args.decimate_by is not None:  # smooth again and decimate
        zs, mask1=smooth_corrected(z, mask, w_smooth)
        zs[mask1 < .25]=0
        mask[mask1 < .25]=0
        if args.error_RMS_scale is not None:
            r2=(z-zs)**2
            r2[mask==0]=0
            r2, dummy=smooth_corrected(r2, mask, w_smooth)
            r2[mask==0]=0
            r2=r2[np.int(dec/2+1/2)::dec, np.int(dec/2+1/2)::dec]
        zs=zs[np.int(dec/2.+0.5)::dec, np.int(dec/2.+0.5)::dec]  # this is now the same res as the output image, includes pad
        z=zs
        mask=mask[np.int(dec/2+0.5)::dec, np.int(dec/2+0.5)::dec]
        z[mask==0]=0

    if args.geolocation_error is not None:
        gx, gy=np.gradient(zs, dec*dx, dec*dx)
        mask1=snd.binary_erosion(mask, np.ones((3,3)))
        gxs, mask_x=smooth_corrected(gx, mask1, 4)
        gys, mask_y=smooth_corrected(gy, mask1, 4)
        edge_mask=(mask==1) & (mask1==0) & (mask_x >.25)
        gx[mask1==0]=0
        gy[mask1==0]=0
        gx[ edge_mask ] = gxs[edge_mask]
        gy[ edge_mask ] = gys[edge_mask ]
        e2_geo=(gx**2+gy**2)*args.geolocation_error**2
        e2_geo[mask==0]=0
    else:
        e2_geo=np.zeros_like(zs)


    if args.error_RMS_scale is not None:
        zss, mask2=smooth_corrected(z, mask, w_error)
        r2s, dummy=smooth_corrected(e2_geo+(zss-z)**2, mask, w_error)
        r2, dummy=smooth_corrected(r2, mask, w_error)
        error_est=np.sqrt(r2s+r2)
        error_est[mask==0]=np.NaN
        out_temp[1,:,:]=error_est[pad:-pad, pad:-pad]
    else:
        out_temp[1,:,:]=e2_geo[pad:-pad, pad:-pad]

    z[mask==0]=np.NaN
    out_temp[0,:,:]=z[pad:-pad, pad:-pad]

    if DOPLOT is not None:
        plt.figure
        plt.subplot(2,2,1)
        plt.imshow(in_sub.z[0,:,:]); plt.colorbar()
        plt.subplot(2,2,2)
        plt.imshow(mask)
        plt.subplot(2,2,3)
        plt.imshow(error_est); plt.colorbar()
        plt.show()

    out_sub.z=out_temp
    out_sub.setBounds(out_sub.c0+pad, out_sub.r0+pad, out_sub.Nc-2*pad, out_sub.Nr-2*pad)
    out_sub.writeSubsetTo(out_bands, out_sub)
    delta_time=time.time()-last_time
    sys.stdout.write("\r\b %d out of %d, last dt=%f" %(out_sub.count, out_sub.xy0.shape[0], delta_time))
    sys.stdout.flush()
    last_time=time.time()

outDs.SetGeoTransform(tuple(xform_out))
for b in out_bands:
    outDs.GetRasterBand(b).SetNoDataValue(np.NaN)
outDs.SetProjection(ds.GetProjection())
outDs=None

