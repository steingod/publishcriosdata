# SPICE

The fetch_spice.py script allows you to download SPICE data and save it as a netCDF file. In order to run the script you must have installed the necessary dependencies. A requirements.txt file is provided to make this easy. The following is an example of how the script can be run from the command line:

`python fetch_spice.py SPICE38 2023-09-01T00:00:00Z 2023-09-30T23:59:59Z test.nc`

Running the script with the -h flag shows what all of the different parameters are:
```
python fetch_spice.py -h
usage: fetch_spice.py [-h] {SPICE34,SPICE35,SPICE36,SPICE37,SPICE38} start_ts end_ts path

Download and save SPICE data.

positional arguments:
  {SPICE34,SPICE35,SPICE36,SPICE37,SPICE38}
                        SPICE station code. Can be SPICE34, SPICE35, SPICE36, SPICE37 or SPICE38.
  start_ts              Start time in UTC (YYYY-MM-DDTHH:MM:SSZ).
  end_ts                End time in UTC (YYYY-MM-DDTHH:MM:SSZ).
  path                  Path (including filename) to where netCDF file will be saved.

options:
  -h, --help            show this help message and exit
```
