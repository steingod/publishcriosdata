#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import requests
import json
from json.decoder import JSONDecodeError
import numpy as np
import xarray as xr
import datetime
import dateutil
import argparse
import yaml
import logging
import logging.handlers
from calendar import monthrange
import uuid

def parse_arguments():
    """
    Parsing of command line arguments. Using a mix of positional and optional arguments.
    """
    parser = argparse.ArgumentParser(
            description='Download and save SPICE data. If station is specified, output has to be specified as well. If start_ts is not specified, the software enters update mode and will check for updates since the last file was created (not in single station mode).', 
            epilog=' The connection to the remote server is specified in the configuration file. Positional arguments will override potential station lists provided in the configuration file. If no end time is provided, the current time is used. If station is specified in the command line, the destination path and filename is also required.')
    parser.add_argument("-c","--cfg",dest="cfgfile",
            help="Configuration file", required=True)
    parser.add_argument("-s", "--station", dest='station',
            help='Station code in the API. Only used for testing retrieval of one station, else configuration file is used.', required=False)
    parser.add_argument('start_ts', nargs="?", default=argparse.SUPPRESS,
            help='Start time in UTC (YYYY-MM-DDTHH:MM:SSZ).')
    parser.add_argument('end_ts', nargs="?", default=argparse.SUPPRESS,
            help='End time in UTC (YYYY-MM-DDTHH:MM:SSZ). Can be excluded and will then be set to now.')
    parser.add_argument('-o', '--output', dest='output', 
            help='Filename (including path) to where the CF-NetCDF file will be saved.', required=False)

    args = parser.parse_args()

    if args.cfgfile is False:
        parser.print_help()
        parser.exit()

    if args.station and not args.output:
        parser.print_help()
        parser.exit()

    return args

def parse_cfg(cfgfile):
    """
    Read the configuration file. This is based on YAML and have a structure like:
    criosapicfg:
        sttype: spice|aws
        endpoint: <URL>
        stations: <list of stations>
            stid: <station id>
                attributes following ACDD
    output:
        logfile: <location>
        destdir: <destination folder>
        default_metadata: <are overriden by information under each station>
            attributes following ACDD
    """
    with open(cfgfile, 'r') as ymlfile:
        try:
            cfgstr = yaml.full_load(ymlfile)
        except Exception as e:
            print('Error parsing YAML ', e)
            raise Exception("Couldn't parse the YAML file")

    return cfgstr

def initialise_logger(outputfile = './log'):
    """
    Set up logging to a rotated logfile.
    """
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

def split_period(mylog, start_ts, end_ts):
    """
    Split the requested time period into monthly sections for better granularity in the file system.
    """

    mylog.info("Splitting requested time period into monthly segments")
    
    print(start_ts + ' - ' + end_ts)

    # Convert from string to something useful
    try:
        mystart = datetime.datetime.strptime(start_ts, "%Y-%m-%dT%H:%M:%SZ")
    except Exception as e:
        mylog.warning("Conversion of start time string to datetime failed: ", e)
        raise Exception("Error in time specification")
    try:
        myend = datetime.datetime.strptime(end_ts, "%Y-%m-%dT%H:%M:%SZ")
    except Exception as e:
        mylog.warning("Conversion of end time string to datetime failed: ", e)
        raise Exception("Error in time specification")

    # Get month span
    mydelta = dateutil.relativedelta.relativedelta(myend, mystart)
    numberofperiods = mydelta.years*12+mydelta.months+1

    # Loop over time periods to create
    if numberofperiods > 1:
        endtime = datetime.datetime(mystart.year,mystart.month,monthrange(mystart.year,mystart.month)[1],23,59,59,0)
        periods = [{'start': start_ts, 'end': endtime.strftime("%Y-%m-%dT%H:%M:%SZ")}]
        for mymonth in range(2,numberofperiods+1):
            starttime = endtime+datetime.timedelta(0,1)
            endtime += datetime.timedelta(monthrange(starttime.year, starttime.month)[1],0)
            myperiod = {'start': starttime.strftime("%Y-%m-%dT%H:%M:%SZ"), 'end': endtime.strftime("%Y-%m-%dT%H:%M:%SZ")}
            periods.append(myperiod)
    else:
        periods = [{'start': start_ts, 'end': end_ts}]
                
    return periods

def retrievedata(mylog, station_type, endpoint, station, start_ts, end_ts):
    """
    Connect to the remote API, download data and return these as JSON
    """

    mylog.info("Now in retrievespicedata")
    mylog.info("Collecting data for %s - %s", start_ts, end_ts)

    # Set up the request depending on station type
    if station_type == 'aws':
        return
        myrequest = endpoint + '?id='+station+'&startTs='+start_ts+'&endTs='+end_ts+'&meta=yes&awsId=yes&time=yes&setf=no'
    elif station_type == 'spice':
        myrequest = endpoint + '?spice='+station+'&startTs='+start_ts+'&endTs='+end_ts+'&meta=yes&spiceId=yes&time=yes'
    else:
        mylog.error('The station_type is not supported %s', station_type)
        raise Exception('Unsupported station type')

    #print(myrequest)
    try:
        response = requests.get(myrequest)
    except Exception as e:
        raise Exception("Couldn't connect to the end point.")

    # Check the returned document
    if response.text == '0 results':
        mylog.warning('This station does not seem to have any data for the time period specified.')
        return None

    # Load data as JSON for further processing
    try:
        json_data = json.loads(response.text)
    except JSONDecodeError:
        mylog.warning('Could not load data as JSON')
        raise Exception('Could not load data as JSON')

    return json_data

def transformdata(mylog, station_type, myjson, md):
    """
    Select transformation approach for the returned JSON payload.
    """
    mylog.info('Now transforming data for station type %s', station_type)
    if station_type == 'spice':
        try:
            mydata = spice2ds(mylog, myjson, md)
        except Exception as e:
            mylog.warning('Something failed when converting SPICE JSON to XARRAY Dataset: %s', e)
            raise Exception("Could not tranform the SPICE data retrieved.")
    elif station_type == 'aws':
        try:
            mydata = aws2ds(mylog, myjson, md)
        except Exception as e:
            mylog.warning('Something failed when converting AWS JSON to XARRAY Dataset: %s', e)
            raise Exception("Could not tranform the AWS data retrieved.")
    else:
        raise Exception('Station type is not recognised')

    return mydata

def aws2ds(mylog, json_data, md):
    """
    Convert AWS data to XARRAY dataset
    """

    return

def spice2ds(mylog, json_data, md):
    """
    Convert SPICE data to XARRAY dataset
    """

    # Variable specifications and conversions
    variables = (
            'surface_snow_thickness', 
            'air_temperature', 
            'relative_humidity', 
            'air_pressure'
            )
    long_names = {
            'surface_snow_thickness': 'snow thickness measured at station',
            'air_temperature': 'air temperature measured at station',
            'relative_humidity': 'relative humidity measured at station',
            'air_pressure': 'air pressure measured at station'
            }
    data_arrays = {}

    for variable in variables:
        # Must read time in second precision to get correct values, and then convert to nanosecond precision to
        # silence xarray warning about non-nanosecond precision.
        time = np.array([item['t'] for item in json_data[variable]], dtype='datetime64[s]').astype('datetime64[ns]')
        data = np.array([item['v'] for item in json_data[variable]], dtype=float)

        data_arrays[variable] = xr.DataArray(data=data,
                dims=['time'],
                coords=dict(time=time),
                attrs=dict(long_name=long_names[variable],
                    standard_name=variable,
                    units=json_data['metadata'][variable]['v']['units'],
                    coverage_content_type='physicalMeasurement'))

        try:
            _FillValue = json_data['metadata'][variable]['v']['_FillValue']
        except KeyError:
            pass
        else:
            # Replace any null values with _FillValue if it exists.
            data_arrays[variable] = data_arrays[variable].where(data_arrays[variable].notnull(), _FillValue)
            data_arrays[variable].attrs['_FillValue'] = _FillValue

    ds = xr.Dataset(data_vars=data_arrays)
    ds['time'].attrs = {'standard_name': 'time', 'long_name': 'time of observation'}

    # Positions of the stations. SPICE37 is currently not set up, so it has -999 as a placeholder for lat, lon. 
    lat = md['md']['geospatial_lat_min']
    lon = md['md']['geospatial_lon_min']

    ds = ds.assign(lat=lat)
    ds['lat'].attrs = {'units': 'degree_north', 'standard_name': 'latitude'}

    ds = ds.assign(lon=lon)
    ds['lon'].attrs = {'units': 'degree_east', 'standard_name': 'longitude'}

    # Prepare attributes
    creationtime = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    # TODO: check whether attributes are OK.
    ds.attrs = {}
    # id is set when writing to file
    if 'title' in md['md'].keys():
        ds.attrs['title'] = md['md']['title']
    else:
        ds.attrs['title'] = f'Measurement data from {md["stid"]} station'
    if 'summary' in md['md'].keys():
        ds.attrs['summary'] = md['md']['summary']
    else:
        ds.attrs['summary'] = ('Data from SPICE stations includes snow depth measurements, surface temperature, relative humidity, and atmospheric pressure. These stations and their data are a part of the CRIOS project.')
    if 'keywords' in md['md'].keys():
        ds.attrs['keywords'] = md['md']['keywords']
    else:
        ds.attrs['keywords'] = ('GCMDSK:EARTH SCIENCE > CRYOSPHERE > SNOW/ICE > SNOW DEPTH,'  +
        'GCMDSK:EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC TEMPERATURE > SURFACE TEMPERATURE > AIR TEMPERATURE,' +
        'GCMDSK:EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC WATER VAPOR > WATER VAPOR INDICATORS > HUMIDITY > RELATIVE HUMIDITY,' +
        'GCMDSK:EARTH SCIENCE > ATMOSPHERE > ATMOSPHERIC PRESSURE > ATMOSPHERIC PRESSURE MEASUREMENTS')
    if 'keywords_vocabulary' in md['md'].keys():
        ds.attrs['keywords_vocabulary'] = md['md']['keywords_vocabulary']
    else:
        ds.attrs['keywords_vocabulary'] = 'GCMDSK: GCMD Science Keywords'
    if 'geospatial_lat_min' in md['md'].keys():
        ds.attrs['geospatial_lat_min'] = md['md']['geospatial_lat_min']
    else:
        ds.attrs['geospatial_lat_min'] = str(lat)
    if 'geospatial_lat_max' in md['md'].keys():
        ds.attrs['geospatial_lat_max'] = md['md']['geospatial_lat_max']
    else:
        ds.attrs['geospatial_lat_max'] = str(lat)
    if 'geospatial_lon_min' in md['md'].keys():
        ds.attrs['geospatial_lon_min'] = md['md']['geospatial_lon_min']
    else:
        ds.attrs['geospatial_lon_min'] = str(lon)
    if 'geospatial_lon_max' in md['md'].keys():
        ds.attrs['geospatial_lon_max'] = md['md']['geospatial_lon_max']
    else:
        ds.attrs['geospatial_lon_max'] = str(lon)
    ds.attrs['Conventions'] = 'CF-1.11, ACDD-1.3'
    if 'creator_type' in md['md'].keys():
        ds.attrs['creator_type'] = md['md']['creator_type']
    else:
        ds.attrs['creator_type'] = 'person'
    if 'creator_institution' in md['md'].keys():
        ds.attrs['creator_institution'] = md['md']['creator_institution']
    else:
        ds.attrs['creator_institution'] = 'University of Silesia'
    if 'creator_name' in md['md'].keys():
        ds.attrs['creator_name'] = np.array(md['md']['creator_name'].encode("utf-8"))
    else:
        ds.attrs['creator_name'] = np.array('Łukasz Małarzewski'.encode("utf-8"))
    if 'creator_email' in md['md'].keys():
        ds.attrs['creator_email'] = md['md']['creator_email']
    else:
        ds.attrs['creator_email'] = 'lukasz.malarzewski@us.edu.pl'
    if 'creator_url' in md['md'].keys():
        ds.attrs['creator_url'] = md['md']['creator_url']
    else:
        ds.attrs['creator_url'] = 'https://us.edu.pl/instytut/inoz/en/osoby/malarzewski-lukasz/'
    if 'project' in md['md'].keys():
        ds.attrs['project'] = md['md']['project']
    else:
        ds.attrs['project'] = 'CRIOS'
    if 'license' in md['md'].keys():
        ds.attrs['license'] = md['md']['license']
    else:
        ds.attrs['license'] = 'https://spdx.org/licenses/CC-BY-4.0 (CC-BY-4.0)'
    ds.attrs['iso_topic_category'] = 'climatologyMeteorologyAtmosphere'
    ds.attrs['activity_type'] = 'In Situ Land-based station'
    if 'operational_status' in md['md'].keys():
        ds.attrs['operational_status'] = md['md']['operational_status']
    else:
        ds.attrs['operational_status'] = 'Not available'
    ds.attrs['time_coverage_start'] = np.datetime_as_string(ds.time[0], timezone='UTC', unit='s')
    ds.attrs['time_coverage_end'] = np.datetime_as_string(ds.time[-1], timezone='UTC', unit='s')
    ds.attrs['history'] = creationtime + ': created file'
    ds.attrs['date_created'] = creationtime
    ds.attrs['featureType'] = 'timeSeries'

    #print(json.dumps(ds.attrs, indent=2))
    float_type = {'dtype': 'float32'}
    integer_type = {'dtype': 'int32'}

    encoding = {'time': integer_type,
                'surface_snow_thickness': float_type,
                'air_temperature': float_type,
                'relative_humidity': float_type,
                'air_pressure': float_type,
                'lat': float_type,
                'lon': float_type}

    return {'ds': ds, 'enc': encoding}

def createMETuuid(infile):
    """
    Create a unique identifier for the dataset. Here we use the approach from mdharvest with data centre MET/ADC since data are published there
    """
    # Prepare creation of UUID
    filename = "https://arcticdata.met.no/ds/"+os.path.basename(infile)+"-"
    filename += datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    # Create UUID
    myuuid = uuid.uuid5(uuid.NAMESPACE_URL,filename)
    # Add namespace
    myidentifier = 'no.met.adc:'+str(myuuid)

    return myidentifier

def ds2netcdf(mylog, ds, filename):
    mylog.info('Now dumping data to file')
    # First check if file already exists
    oldhist = None
    oldid = None
    if os.path.isfile(filename):
        mylog.info("%s exists, will retrieve existing identifier", filename)
        # The open the file and retrieve existing identifier
        try:
            oldds = xr.open_dataset(filename)
        except Exception as e:
            mylog.warning("Could not open old dataset: %s", e)
            mylog.warning("Overwriting...")
        if 'id' in oldds.attrs.keys():
            oldid = oldds.attrs['id']
        else:
            mylog.info('No id found, adding one')
        if 'history' in oldds.attrs.keys():
            oldhist = oldds.attrs['history']
        oldds.close()
    # Add the identifier
    if oldid:
        ds['ds'].attrs['id'] = oldid
    else:
        ds['ds'].attrs['id'] = createMETuuid(filename)
    # Fix the history
    if oldhist:
        ds['ds'].attrs['history'] = oldhist+'\n'+datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')+': Updated'
    # Dump the data
    try:
        ds['ds'].to_netcdf(filename, encoding=ds['enc'], unlimited_dims='time', mode="w")
    except Exception as e:
        mylog.warning('Something failed when creating NetCDF file: %s', e)
        raise Exception('Something failed when creating NetCDF file')
    mylog.info("%s is created", filename)

    return

if __name__ == '__main__':

    # parse command line arguments
    try:
        myargs = parse_arguments()
    except:
        raise SystemExit('Command line arguments didn\'t parse correctly.')

    # Parse configuration file
    if myargs and 'cfgfile' in vars(myargs):
        try:
            cfgstr = parse_cfg(myargs.cfgfile)
        except Exception as e:
            raise Exception('Could not parse the configuration file.')

    # Initialise logging
    if 'output' in cfgstr.keys() and 'logfile' in cfgstr['output'].keys():
        mylog = initialise_logger(cfgstr['output']['logfile'])
        mylog.info('Configuration of logging is finished.')
    else:
        raise Exception('No output section or logfile defined in configuration file.')

    # If end time is not provided, set it now
    if 'end_ts' not in vars(myargs):
        end_ts = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
    else:
        end_ts = myargs.end_ts

    # Loop over stations in the configuration file or do the single stations requested in command line
    sttype = cfgstr['criosapicfg']['sttype']
    endpoint =cfgstr['criosapicfg']['endpoint'] 
    if myargs.station:
        """
        Handling single station requests in one go
        """
        # Retrieve data
        try:
            mydata = retrievedata(mylog, sttype, endpoint, myargs.station, myargs.start_ts, end_ts)
        except Exception as e:
            raise SystemExit("Retrival and conversion didn't work", e)
        # Transform to Dataset
        try:
            myds = transformdata(mylog, sttype, mydata)
        except Exception as e:
            mylog.warning("Something failed in transformdata")
        # Write to file
        try:
            ds2netcdf(mylog, myds, args.output)
        except Exception as e:
            mylog.warning('Something failed when dumping to NetCDF file: %s', e)
            raise Exception('Something failed when dumping to NetCDF')
    else:
        if 'stations' in cfgstr['criosapicfg']:
            # If start time is not provided, get it from last generated files
            # FIXME
            # Split requested time period into monthly periods
            try:
                myperiods = split_period(mylog, myargs.start_ts, end_ts)
            except Exception as e:
            mylog.error("Inconcistency in time specification: %s", e)
            raise SystemExit("Inconcistency in time specification")

            for mystation in cfgstr['criosapicfg']['stations'].keys():
                mymd = {
                        'stid': mystation,
                        'md': cfgstr['criosapicfg']['stations'][mystation]
                        }
                mylog.info("")
                mylog.info("Processing station: %s", mystation)
                # Loop time periods
                for p in myperiods:
                    mylog.info("Processing period: %s - %s", p['start'], p['end'])
                    # Retrieve data
                    try:
                        mydata = retrievedata(mylog, sttype, endpoint, mystation, p['start'], p['end'])
                    except Exception as e:
                        mylog.warning("Something failed in retrievedata")
                    if mydata == None:
                        continue
                    # Transform to Dataset
                    try:
                        myds = transformdata(mylog, sttype, mydata, mymd)
                    except Exception as e:
                        mylog.warning("Something failed in transformdata")
                    # Write to file
                    outfile = cfgstr['output']['destdir']+'/'+mystation+'-'+datetime.datetime.strptime(p['start'],'%Y-%m-%dT%H:%M:%SZ').strftime('%Y%m')+'.nc'
                    try:
                        ds2netcdf(mylog, myds, outfile)
                    except Exception as e:
                        mylog.warning('Something failed when dumping to NetCDF file: %s', e)
                        raise Exception('Something failed when dumping to NetCDF')

