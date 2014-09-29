#! /usr/bin/env python

# ben - Sep 26 2014
# script to parse the output of split_wv_delivery_by_order and generate a set of directories
# and links that will allow straightforward processing of stereo data.
# given:
# data_dir : the directory containing ntf and xml files for the pair
# pair : the worldview pair 
# proc_dir : a directory to hold processing files and WV output
# (optional) split_file : the output of split_wv_delivery_by_order
# make (if DNE):
# proc_dir 
# data and mosaic directories for each ID and order:
# proc_dir/(ID1)_(order1)
# proc_dir/(ID1)/(order1)/data
#      links to ntf and xml files
# proc_dir/(ID1)/(order1)/mosaic
#      wv_correct output and mosaics 
# proc_dir/(ID1)_(order2)  ...etc
#
# processing directories for each pair of orders:
# proc_dir/(pair)_(order1)_(order2)
#       links to mosaics 


import os, argparse, re, sys

parser=argparse.ArgumentParser(description='generate links for a stereo pair')
parser.add_argument('data_dir', type=string, help='data directory root')
parser.add_argument('pair', type=string, help='pair to parse')
parser.add_argument('proc_dir', type=string, help='directory root for the processing.  The pair directory will be created within this directory.', default='.')
parser.add_argument('split_file', type=string, help='optional path to the output of split_wv_xml', default=None)

args=parser.parse_args()

if args.pair is None:
    args.pair=os.path.basename(os.path.getcwd())

if not os.path.isdir(pair):
    os.makedirs(pair)

data_source_dir=args.data_dir+"/"+pair

if not os.path.islink(pair+"/base"):
    os.symlink(data_source_dir, pair+"/base")

if args.split_file is None:
    split_file=pair+"/base/split_info.txt"
if not os.path.isfile(split_file):
    os.system("split_wv_by_delivery " +data_source_dir+" > "+split_file)

order_re=re.compile('# OrderSet_(\d+): Order1=(\S+).*Order2=(\S+)');
orderpair_list=list()
ID_order_list=list()
for line in open( data_source_dir+"/split_info.txt" ):
    fields=line.split(' ')
    if line.startswidth("#"):
        temp=order_re.match(line)
        if temp is not None:
            orderpair=temp.group(1)
            order1=temp.group(2)
            order2=temp.group(3)
            orderpair_list.append([orderpair, order1, order2])
    else:
        ID=fileds[0]
        order=fields[1]
        ID_order_list.append([ID, order, fields[2:]]]


# for each ID in the list, and each order in the ID, make a directory: 
# data_dir/ID1/order1
# within that directory, make a /data directory, populate it with symlinks to the xmls and ntfs for that order. 
for temp in ID_order_list:
    this_ID_dir=args.proc_dir+"/"+pair+"/", temp[0]
    if not os.isdir(this_ID_dir):
        os.makedirs(this_ID_dir)
    this_order_dir=this_ID_dir+"/"+temp[1]
    if not os.isdir(this_order_dir):
        os.makedirs(this_order_dir)
    this_data_dir=this_order_dir+"/data"
    if not os.isdir(this_data_dir):
        os.makedirs(this_data_dir)
    for this_fname in temp[1]:
        for extension in ('.xml','.ntf'):
            this_link_to=this_data_dir+'/'+this_fname+ext
            this_link_from=data_source_dir+'/'+this_fname+ext
            if not os.path.islink(this_link_to):
                os.symlink(this_link_from, this_link_to)
# Will need to run a script to call wv_correct and dg_mosaic to make the mosaics


# make one processing directory for each pair of orders,
# link to the mosaicked tif and xml files (even though they don't yet exist)
# 
for temp in orderpair_list:
    this_pair_dir=args.proc_dir+pair+"/"+pair"_"+temp[1]+"_"+temp[2]
    if not os.path.isdir(this_pair_dir):
        os.makedirs(this_pair_dir)
    for ind in [0, 1]:
        for ext in [".r100.tif",".r100.xml"]:
            this_link_to=this_pair_dir+"/"+ID_order_list[ind][0]+ext
            this_link_from=data_souce_dir+"/"+ID_order_list[ind+1]+"/"+ID_order_list[ind][0]+ext
            if not os.path.islink(this_link_to):
                os.symlink(this_link_from, this_link_to)

