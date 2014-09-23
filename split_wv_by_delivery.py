#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:04:35 2014

@author: ben
"""


import re, sys
from osgeo import ogr,  osr
import os, glob,  subprocess, shutil
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import  datetime
import xml.etree.ElementTree as ET
ll_refSys=osr.SpatialReference()
ll_refSys.ImportFromEPSG(4326)
#print "done importing"
DEBUG=False

def read_xml_data(filename):
    # use the Elementtree package to read data from the Worldview XMLs
    IMD=dict()
    tree=ET.parse(filename) 
    IMD_root=tree.getroot().find('IMD')
    ProdOrder=IMD_root.find('PRODUCTORDERID').text
    
    IMD['orderID']=ProdOrder.split('_')[0]
    BAND_P=IMD_root.find('BAND_P')
    IMD['lonlat']=np.zeros([4, 2])
    for count, corner in enumerate(('UL', 'LL', 'LR', 'UR')):
        IMD['lonlat'][count,1]= np.float(BAND_P.find(corner+'LAT').text)
        IMD['lonlat'][count,0]= np.float(BAND_P.find(corner+'LON').text)
    IMAGE=IMD_root.find('IMAGE')
    IMD['GSD']=np.float(IMAGE.find('MEANPRODUCTGSD').text)
    EPH=tree.find('EPH')
    IMD['STARTTIME']=EPH.find('STARTTIME').text
    return IMD
    
def xform_pts(xform, A,B=None,C=None):
    # utility function to transform an array of coordinates from one CS to another
    if B is None:
        if A.shape[1]==3:
            C=A[:,2]
        B=A[:,1]
        A=A[:,0]
    if C is None:
        C=np.zeros_like(A)
    qrs= [xform.TransformPoint(A[k], B[k], C[k]) for k in range(0, A.shape[0])]
    qrs=np.array(qrs)
    return qrs

def bVec_2d(lat, lon, xform):
    # utility function to calculate the 2-d basis vectors for a lat/lon pair in a specified CS
    x0=xform_pts(xform, np.array([[lon, lat, 0]]))
    x1=xform_pts(xform, np.array([[lon, lat-.001, 0]]))
    N_hat=(x1-x0).ravel();
    N_hat=N_hat[0:2]/np.sqrt(N_hat[1]**2+N_hat[0]**2)
    x1=xform_pts(xform, np.array([[lon-.001, lat, 0]]))
    E_hat=(x1-x0).ravel()
    E_hat=E_hat[0:2]/np.sqrt(E_hat[1]**2+E_hat[0]**2)
    return N_hat, E_hat


def calc_at_proj4_for_ll(lonlat):
    # utility that returns a proj4 string for a coordinate system whose y axis is
    # aligned with the long axis of a worldview pair, given longitude-latitude points
    # defining the data extent.  Assumes that the first two coordinates are on the
    # left-hand boundary of the pair
    lonlat0=lonlat[1,:]
    ps_standard_sys=osr.SpatialReference()
    if lonlat0[1] < 0:
        # Antarctica:
        # standard system and mask cs are AA polar stereo
        ps_standard_sys.ImportFromEPSG(3031)
    else:
        # Arctic:
        # standard PS system is GL polar stereo, mask is in lonlat
        ps_standard_sys.ImportFromEPSG(3413)
    ll_sys=osr.SpatialReference()
    ll_sys.ImportFromEPSG(4326)
    ll2ps_xform=osr.CoordinateTransformation(ll_sys, ps_standard_sys)    
    
    xy1=xform_pts(ll2ps_xform, lonlat)
    E_hat, N_hat=bVec_2d(lonlat0[1], lonlat0[0], ll2ps_xform)

    #xyz1=xform_pts(ll2ecf, lonlat1)
    dN=np.dot(N_hat, xy1[1,0:2]-xy1[0,0:2])
    dE=np.dot(E_hat, xy1[1,0:2]-xy1[0,0:2])
    norm=np.sqrt(dN**2.+dE**2.)
    dN=dN/norm
    dE=dE/norm
    if lonlat0[1] < 0 :
        if DEBUG:
            print "South!!!"
        lon0=lonlat0[0]-np.arctan2(dN, dE)*180./np.pi
    else:
        if DEBUG:
            print "North!!!"
        lon0=lonlat0[0]+np.arctan2(dN, dE)*180./np.pi
    # make sure the standard longitude isn't close to the center of the map.  Rotate by 180 degrees if it is    
    if np.abs(np.mod(lon0+180.-lonlat0[0], 360.)) < 10. or np.abs(np.mod(lon0+180.-lonlat0[0], 360.)) > 350.:
        lon0=np.mod(lon0+180., 360.)
    if lon0 > 180.:
        lon0=lon0-360.
    if lon0 < -180.:
        lon0=lon0+360.
    proj4_str='+proj=sterea +lat_0=%f +lon_0=%f +k_0=1.0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' % (lonlat0[1], lon0)
    return proj4_str, lon0, lonlat0[1], ps_standard_sys

def Intersect_Remainder(ID_data, Remainder):
    # Given:
    #  ID_data : a 2-element list of data structures defining the geometry of the data for each ID in each order
    #  Remainder : A geometry defining the remaining region to be mapped
    # 
    # find the order-pair intersection that covers the largest area of Remainder
    # Return:
    #  bestOrders: a 2-element list of the orders in ID1 and ID2 whose intersection covers the largest area
    #  bestA : the area of the intersection
    #  Intersection_poly : a geometry defining the intersection
    #  Remainder : a geometry defining the remaining region to be mapped
    #  Remainder_A: the area of Remainder
    
    keys1=ID_data[0]['Order_dict'].keys()
    keys2=ID_data[1]['Order_dict'].keys()
    A0=np.zeros([len(keys1), len(keys2)])
    for Ord1, key1 in enumerate(keys1):
        for Ord2, key2 in enumerate(keys2):
            P=ID_data[0]['Order_dict'][key1]['poly']
            P=P.Intersection(ID_data[1]['Order_dict'][key2]['poly'])
            if P is not None:
                P=P.Intersection(Remainder)
                # make sure the intersection has a valid area
                if (np.uint16(P.GetGeometryType())==ogr.wkbPolygon or np.uint16(P.GetGeometryType())==ogr.wkbMultiPolygon or np.uint16(P.GetGeometryType())==ogr.wkbGeometryCollection) :        
                    A0[Ord1, Ord2]=P.Area()
                else:
                    A0[Ord1, Ord2]=0.0
            else:
                P=None
                A0[Ord1, Ord2]=0.
    if np.any(A0.ravel()>0):
        bestOrderSub=np.unravel_index(np.argmax(A0), A0.shape)
        best_A=A0[bestOrderSub[0], bestOrderSub[1]]
        bestOrders=[keys1[bestOrderSub[0]], keys2[bestOrderSub[1]]]
        Intersection_poly=ID_data[0]['Order_dict'][bestOrders[0]]['poly'].Intersection(ID_data[1]['Order_dict'][bestOrders[1]]['poly']).Intersection(Remainder)
        Remainder=Remainder.Difference(Intersection_poly)
        if Remainder is not None:
            Remainder_A=Remainder.Area()
        else:
            Remainder_A=0.
    else:
        bestOrders=None
        best_A=0
        Intersection_poly=None
        Remainder=None
        Remainder_A=0
    return bestOrders, best_A, Intersection_poly, Remainder, Remainder_A

def getGeom_xy(geom):
    # Utility that returns a list lists of coordinates for the polygons in geom
    x=list()
    y=list()
    for i in range(geom.GetGeometryCount()):
        r=geom.GetGeometryRef(i)
        if r.GetPointCount()==0:
            # this is a multigeometry, recurse into it
            xi, yi=getGeom_xy(r)
            for count,xii in enumerate(xi):
                x.append(xii)
                y.append(yi[count])
            continue
        x.append([r.GetX(j) for j in range(r.GetPointCount())])
        y.append([r.GetY(j) for j in range(r.GetPointCount())])
    return x, y


def make_geom_shapefile(shpName, geom, field_dict, out_ref_sys=None, in_ref_sys=None):
    # output a shapefile for geom with fields matching field_dict
    if out_ref_sys is None:
        # if nothing else is specified, use latlon
        out_ref_sys = osr.SpatialReference()
        out_ref_sys.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

    if in_ref_sys is not None:
        xform=osr.CoordinateTransformation(in_ref_sys, out_ref_sys)
    else:
        xform=None
    # Define shape and layer
    driver = ogr.GetDriverByName('ESRI Shapefile')
 
    # delete the shape file, etc if it exists
    old_files=glob.glob(shpName+'.*')
    if len(old_files) > 0:
        for old_file in old_files:
            os.remove(old_file)

    shapeData=driver.CreateDataSource(shpName+'.shp')
    layer = shapeData.CreateLayer(field_dict['name'], out_ref_sys, ogr.wkbMultiPolygon)
    # copy the field_dict into 
    for key in field_dict.keys():
        this_val=field_dict[key]
        if isinstance(this_val,float):
            field_defn= ogr.FieldDefn(key , ogr.OFTReal )
        if isinstance(this_val,int):
            field_defn= ogr.FieldDefn(key,  ogr.OFTInteger)
        if isinstance(this_val,basestring):
            field_defn = ogr.FieldDefn( key, ogr.OFTString )
            field_defn.SetWidth( 60 )
        layer.CreateField ( field_defn ) 
             
    xx, yy=getGeom_xy(geom)
    myPoly = ogr.Geometry(type=ogr.wkbPolygon)            
    for xi, yi in zip(xx, yy):
        myRing = ogr.Geometry(type=ogr.wkbLinearRing)
        for xp, yp in zip(xi, yi):
            if xform is not None:
                xy=xform.TransformPoint(xp, yp)
                xp=xy[0]
                yp=xy[1]
            myRing.AddPoint(xp, yp)
        myRing.CloseRings()
        myPoly.AddGeometry(myRing)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(myPoly)        
    for key in field_dict.keys():
        feature.SetField(key, field_dict[key])
    layer.CreateFeature(feature)
    feature.Destroy()
        
    shapeData.Destroy()
    return
    
def main():    
    if len(sys.argv)==1:
        top_dir=os.getcwd()
    else:
        top_dir=sys.argv[1]

    # extract the IDs making up the pair from the name of top_dir 
    dir_base=os.path.basename(top_dir)
    IDs=re.compile('WV.*_.*_(.*)_(.*)').search(dir_base).groups()

    # find the xmls corresponding to each ID
    xmls_1=glob.glob(top_dir+'/*'+IDs[0]+'*.xml')
    xmls_2=glob.glob(top_dir+'/*'+IDs[1]+'*.xml')
    if os.path.isdir(top_dir+'/bad_fn'):
        xmls_1.extend(glob.glob(top_dir+'/bad_fn/*'+IDs[0]+'*.xml'))
        xmls_2.extend(glob.glob(top_dir+'/bad_fn/*'+IDs[1]+'*.xml'))
    print "# IDs: %s %s" % IDs

    # Extract the outline for the first xml of the first ID
    IMD0=read_xml_data(xmls_1[0])

    # define an along-track CS based on the first xml's outline
    proj4_str, lon_proj_ctr, lat_proj_ctr, ps_standard_sys = calc_at_proj4_for_ll(IMD0['lonlat'])
    AT_sys=osr.SpatialReference()
    AT_sys.ImportFromProj4(proj4_str)
    ll2AT_xform=osr.CoordinateTransformation(ll_refSys, AT_sys)    

    # read the polygons and image information from each xml file, store in dict IMD
    # Build a polygon for each ID and each order from the union of all subscenes in the ID and order
    # Build a polygon for each ID from the union of all subscenes in the ID
    ID_data=list()
    for ID_count, xmls in enumerate([xmls_1, xmls_2]):
        IMD=list()
        I_poly=ogr.Geometry(ogr.wkbPolygon)
        Order_dict=dict()
        for xml_count, xml_file in enumerate(xmls):
            # read this xml file
            IMD.append(read_xml_data(xml_file))
            IMD[xml_count]['xy']=xform_pts(ll2AT_xform, IMD[xml_count]['lonlat'])
            IMD[xml_count]['name']=xml_file
            ring=ogr.Geometry(ogr.wkbLinearRing)    
            for point in IMD[xml_count]['xy']:
                ring.AddPoint(point[0], point[1])
            ring.CloseRings()
            IMD[xml_count]['poly']=ogr.Geometry(ogr.wkbPolygon)  
            IMD[xml_count]['poly'].AddGeometry(ring)

            # add the current polygon to the geometry for this ID
            I_poly=I_poly.Union(IMD[xml_count]['poly'])
            if not Order_dict.has_key(IMD[xml_count]['orderID']):
                print "#found ID: %s" % IMD[xml_count]['orderID']
                # this is the first xml from this order
                Order_dict[IMD[xml_count]['orderID']]=dict()
                Order_dict[IMD[xml_count]['orderID']]['xml_list']=list()
                Order_dict[IMD[xml_count]['orderID']]['poly_list']=list()
                Order_dict[IMD[xml_count]['orderID']]['poly']=IMD[xml_count]['poly']
            else:
                # add the polygon for this file to the geometry for the order
                Order_dict[IMD[xml_count]['orderID']]['poly']=Order_dict[IMD[xml_count]['orderID']]['poly'].Union(IMD[xml_count]['poly'])
            Order_dict[IMD[xml_count]['orderID']]['xml_list'].append(xml_file) 
            Order_dict[IMD[xml_count]['orderID']]['poly_list'].append(IMD[xml_count]['poly'])
        ID_data.append({'ID': IDs[ID_count], 'poly': I_poly,'Order_dict': Order_dict,'IMD': IMD })

    # find the overlap between the geometries for the two IDs
    CompletePoly=ID_data[0]['poly']
    CompletePoly=CompletePoly.Intersection(ID_data[1]['poly'])
    # the complete overlap is the area we want to map.  It becomes 'remainder' to initialize the first iteration
    Remainder=CompletePoly
    Remainder_A=CompletePoly.Area()
    print "# total area = %6.2f km^2" % (Remainder_A/1.e6)

    # bulid a list of intersections between orders
    BestOrderList=list()
    IntPolyList=list()
    if DEBUG:
        xx0,yy0=getGeom_xy(Remainder)
        for i in range(len(xx0)):
            plt.plot(xx0[i], yy0[i],'k')
    #print '# Order1\tOrder2\tArea(km^2)'
    xml_str=''
    set_count=0
    # iterate until the remaining area is less than 1 km^2  (the 'buffer' statement removes small islands and non-polygon geometries)
    while Remainder.Buffer(-500).Area() > 1.e6:
        # find the best pair of orders intersecting the remaining area
        bestOrders, best_A, Intersection_poly, Remainder, Remainder_A=Intersect_Remainder(ID_data, Remainder)
        BestOrderList.append(bestOrders)
        IntPolyList.append(Intersection_poly)
        # make a list of the ID_data for each pair's data
        oDL=[ID_data[0]['Order_dict'][bestOrders[0]], ID_data[1]['Order_dict'][bestOrders[1]]]
        # erode the intersection poly by 500 m to avoid intersections with small scraps (lines, small islands)
        Intersection_poly=Intersection_poly.Buffer(-500)
        order_xml_list=list()
        order_poly_list=list()
        # loop over the IDs
        for oD in oDL:
            this_xml_list=list()
            this_poly_list=list()
            # loop over the subscenes in each ID
            for this_poly, this_xml in zip(oD['poly_list'], oD['xml_list']):
                # build the intersection of the subscenes that contribute to the current overlap area
                temp_poly=Intersection_poly.Intersection(this_poly)
                if temp_poly.Area() > 0:
                    this_xml_list.append(os.path.splitext(os.path.basename(this_xml))[0])
                    this_poly_list.append(this_poly)
            order_xml_list.append(this_xml_list)
            order_poly_list.append(this_poly_list)
        # report the order names and the inersection area (as commented text)
        print "# Order1=%s\tOrder2=%s\tArea= %6.1f km^2" % (bestOrders[0], bestOrders[1], best_A/1.e6)

        # write out a set of shapefiles for this intersection (can be used to clip the data later)
        shp_geom_list=[Intersection_poly.Buffer(1000.), oDL[0]['poly'], oDL[1]['poly']]
        shp_geom_fields=[{'name':'Intersection', 'OrderID': bestOrders[0]+'-'+bestOrders[1]}, {'name':IDs[0], 'OrderID': bestOrders[0]}, {'name':IDs[1], 'OrderID': bestOrders[1]}]
        for shp_geom, shp_fields in zip(shp_geom_list, shp_geom_fields):            
            make_geom_shapefile(bestOrders[0]+'-'+bestOrders[1]+'-'+shp_fields['name'], shp_geom, shp_fields, out_ref_sys=ll_refSys, in_ref_sys=AT_sys)
            
        # build a string listing the xml files in each order for each intersection
        set_count=set_count+1
        xml_str=xml_str+"OrderSet_"+str(set_count)+" "+bestOrders[0]+" "+bestOrders[1]
        for ox in order_xml_list:
            for ofile in ox:
                xml_str=xml_str+" "+ofile
        xml_str=xml_str+"\n"
        
        if DEBUG:
            plt.figure()
            for i in range(len(xx0)):
                plt.plot(xx0[i], yy0[i],'kx-')
            xx,yy=getGeom_xy(ID_data[0]['Order_dict'][bestOrders[0]]['poly'])
            for i in range(len(xx)):
                plt.plot(xx[i], yy[i],'r')
            xx,yy=getGeom_xy(ID_data[1]['Order_dict'][bestOrders[1]]['poly'])
            for i in range(len(xx)):
                plt.plot(xx[i], yy[i],'b')    
            xx,yy=getGeom_xy(Intersection_poly)
            for i in range(len(xx)):
                plt.plot(xx[i], yy[i],'m*-')  
            
            xx,yy=getGeom_xy(Remainder)
            for i in range(len(xx)):
                plt.plot(xx[i], yy[i],'g*-')    
            plt.axis('equal')

    # write out a table of all the xml files in each order
    print "#Order1 Order2 files"
    print xml_str
    if DEBUG:
        plt.show()          
    
            
    return
    
if __name__ == '__main__':
    main()

