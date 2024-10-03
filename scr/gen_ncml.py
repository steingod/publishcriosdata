#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate NCML files for aggregation of CF-NetCDF files containing data from Frost.
"""
__author__ = "Øystein Godøy"
#__copyright__ = "Copyright Info"
__credits__ = ["Øystein Godøy",]
__license__ = """
    This file is part of extractfromfrost.

    This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This file is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
    """
__version__ = "0.0.0"
__maintainer__ = "Øystein Godøy"
__email__ = "steingod@met.no"

import os
import sys
import argparse
import logging
import logging.handlers
import lxml.etree as ET
from netCDF4 import Dataset
import numpy
import pytz
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--dest",dest="destination",
            help="Destination which to add NCML to", required=True)
    parser.add_argument("-l","--log",dest="logdir",
            help="Destination where to put logfiles", required=True)
    parser.add_argument("-o","--overwrite",action='store_true',
            help="Overwrite if NCML is existing")
    args = parser.parse_args()

    """
    if args.cfgfile is False:
        parser.print_help()
        parser.exit()
    """

    return args

def initialise_logger(outputfile = './log'):
    # Check that logfile exists
    logdir = os.path.dirname(outputfile)
    if not os.path.exists(logdir):
        try:
            os.makedirs(logdir)
        except:
            raise IOError
    # Set up logging
    mylog = logging.getLogger()
    mylog.setLevel(logging.INFO)
    #logging.basicConfig(level=logging.INFO, 
    #        format='%(asctime)s - %(levelname)s - %(message)s')
    myformat = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(myformat)
    mylog.addHandler(console_handler)
    file_handler = logging.handlers.TimedRotatingFileHandler(
            outputfile,
            when='w0',
            interval=1,
            backupCount=7)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(myformat)
    mylog.addHandler(file_handler)

    return(mylog)

def traverse_structure(myfolder):

    for item in os.listdir(myfolder):
        mydir = '/'.join([myfolder,item])
        if not os.path.isdir(mydir):
            continue
        if not item.startswith('SN'):
            mylog.warn('Apparently this is not a station folder.')
            continue
        mylog.info('Processing folder: %s', mydir)
        myncmlfile = '/'.join([mydir,item])+'-aggregated.ncml'
        try:
            create_ncml(myncmlfile, mydir)
        except Exception as e:
            mylog.error('Something failed.')
            raise

def create_ncml(myncmlfile, aggdir):
    # Check if NCML already exist
    if os.path.isfile(myncmlfile) and not args.overwrite:
        mylog.warning("%s already exists, won't do anything", myncmlfile)
        return

    # Create NCML
    ns_map = {None:'http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2'}
    root = ET.Element(ET.QName('netcdf'), nsmap = ns_map)
    aggel = ET.SubElement(root, ET.QName('aggregation'))
    aggel.set('dimName','time')
    aggel.set('type','joinExisting')
    aggel.set('recheckEvery','1 day')
    """
    Removed in favor of listing fo files
    scanel = ET.SubElement(aggel, ET.QName('scan'))
    scanel.set('location', aggdir)
    scanel.set('suffix','.nc')
    """
    # Set up specific files to include, assuming data stored in years
    # Set a time stamp into the future...
    mystarttime = 4000000000 
    for item in sorted(os.listdir(aggdir)):
        curdir = '/'.join([aggdir,item])
        # Check content of yearly folder
        if os.path.isdir(curdir):
            # Process files
            for item2 in sorted(os.listdir(curdir)):
                myfile = '/'.join([curdir,item2])
                # Open NetCDF and check content
                if myfile.endswith('.nc'):
                    myncds = Dataset(myfile)
                    tmp = myncds.variables['time'][:]
                    if min(tmp) < mystarttime:
                        mystarttime = min(tmp)
                    tmpstring = ' '.join(str(num) for num in tmp)
                    netcdf = ET.SubElement(aggel, ET.QName('netcdf'))
                    netcdf.set('location',myfile)
                    netcdf.set('coordValue', tmpstring)
                    myncds.close()
    mymeta = ET.SubElement(root, ET.QName('attribute'))
    mymeta.set('name', 'time_coverage_start')
    mymeta.set('value', datetime.fromtimestamp(mystarttime).astimezone(pytz.timezone('UTC')).strftime("%Y-%m-%dT%H:%M:%S%z"))

    # Dump NCML file
    et = ET.ElementTree(root)
    et.write(myncmlfile, xml_declaration=True, encoding='UTF-8', pretty_print=True)

if __name__ == '__main__':
    
    # Parse command line arguments
    try:
        args = parse_arguments()
    except:
        raise SystemExit('Command line arguments didn\'t parse correctly.')

    # Parse configuration file
    #cfgstr = parse_cfg(args.cfgfile)

    # Initialise logging
    #output_dir = cfgstr['output']
    mylog = initialise_logger(args.logdir)
    mylog.info('Configuration of logging is finished.')

    try:
        traverse_structure(args.destination)
    except Exception as e:
        mylog.error('Something failed %s', e)
