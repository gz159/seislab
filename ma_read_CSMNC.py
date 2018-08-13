# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 14:29:21 2018
@author: maqiang,maqiang@iem.ac.cn
Read strong motion record for
China Strong Motion Networks Center(CSMNC),
Institute Engeering Mechanics,
China Earthquake Administrtion
中国地震局工程力学研究所
CSMNC，"new" format files have recorder time and convert to UTC.
"""
# stdlib imports
from datetime import datetime, timedelta
import re
#import os.path

# third party
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.core.trace import Stats
import numpy as np

TEXT_HDR_ROWS = 16
TIMEFMT = '%Y-%m-%d %H:%M:%S.%f'
COLS_PER_LINE = 8
HDR1 = 'STATION:'  #line6
HDR2 = 'COMP. '    #line10



SRC = ('China Strong Motion Networks Center(CSMNC), '
       'Institute Engeering Mechanics, China Earthquake Administrtion')

def get_channel_name(sample_rate, is_acceleration=True,
                     is_vertical=False, is_north=True):
    """Create a SEED compliant channel name.

    SEED spec: http://www.fdsn.org/seed_manual/SEEDManual_V2.4_Appendix-A.pdf

    Args:
        sample_rate (int): Sample rate of sensor in Hz.
        is_acceleration (bool): Is this channel from an accelerometer.
        is_vertical (bool): Is this a vertical channel?
        is_north (bool): Is this channel vaguely pointing north or the channel
                         you want to be #1?
    Returns:
        str: Three character channel name according to SEED spec.

    """
    band = 'H'  # High Broad Band
    if sample_rate < 80 and sample_rate >= 10:
        band = 'B'

    code = 'N'
    if not is_acceleration:
        code = 'H'  # low-gain velocity sensors are very rare

    if is_vertical:
        number = 'Z'
    else:
        number = '2'
        if is_north:
            number = '1'

    channel = band+code+number
    return channel


def is_CSMNC(filename):
    """Check to see if file is a Chinese CSMNC file.

    Args:
        filename (str): Path to possible GNS V1 data file.
    Returns:
        bool: True if CSMNC, False otherwise.
    """
    try:
        with open(filename, 'rt') as f:  
            lines = [next(f) for x in range(TEXT_HDR_ROWS)]
#            print(lines[5],lines[9])
            if lines[5].startswith(HDR1) and lines[9].startswith(HDR2):
                return True
    except Exception:
        return False
    return False

def read_CSMNC(filename):
    """Read  CSMNC strong motion file.
    Args:
        filename (str): Path to data file.
        kwargs (ref): Other arguments will be ignored.
    Returns:
        Stream: Obspy Stream containing three channels of acceleration data
            (cm/s**2).
    """
    if not is_CSMNC(filename):
        raise Exception('%s is not a valid CSMNC file' % filename)

    # Parse the header portion of the file
    with open(filename, 'rt') as f:
        lines = [next(f) for x in range(TEXT_HDR_ROWS)]
        
    CShdr = {}
    hdr = {}
    coordinates = {}
    standard = {}
    standard['units'] = 'acc'
    
    CShdr['record_ro'] = lines[0].split()[0]
    CShdr['origin_time'] = lines[1].split()[1:3]
    CShdr['event_location'] = lines[2].split()[0]+" "\
                              +lines[2].split()[1]+" "\
                              +lines[2].split()[2]
    hdr['location'] = CShdr['event_location']
    CShdr['e_lat'] = float(re.match(r'\d+\.?\d*', lines[3].split()[1]).group())
    CShdr['e_long'] = float(re.search(r'\d+\.?\d*', lines[3].split()[2]).group())
    CShdr['e_depth'] = lines[3].split()[4]
    CShdr['e_mag'] = float(re.search(r'\d+\.?\d*', lines[4].split()[1]).group())
    #hdr['E_Mag_Type'] = re.match(r'Ms', lines[4].split()[1]).group()
    CShdr['station_name'] = lines[5].split()[1]
    hdr['station'] = CShdr['station_name']
    standard['station_name'] = ''
    CShdr['s_lat'] = float(re.match(r'\d+\.?\d*', lines[5].split()[2]).group())
    coordinates['latitude'] = CShdr['s_lat']
    CShdr['s_long'] = float(re.match(r'\d+\.?\d*', lines[5].split()[3]).group())
    coordinates['longitude'] = CShdr['s_long']
    coordinates['elevation'] = np.nan
    CShdr['site_condition']= lines[6].split()[2]  #SOIL ROCK 
    CShdr['instrument_type'] = lines[7].split()[2]
    CShdr['observating_point'] = lines[8].split()[2] #GROUND HOLE STRUCTURE
    CShdr['comp'] = lines[9].split()[1]
    CShdr['record_type'] = lines[10].split()[0]+" "+lines[10].split()[1]
    CShdr['unit'] = lines[10].split()[3]
    CShdr['no_points'] = int(lines[11].split()[3])
    hdr['npts'] = CShdr['no_points']
    CShdr['time_delta'] = float(lines[11].split()[8])
    hdr['delta'] = CShdr['time_delta']
    hdr['sampling_rate'] = float(int(1/hdr['delta']))
    CShdr['pga']= float(lines[12].split()[2])
    CShdr['pga_time']= float(lines[12].split()[4])
    CShdr['duration']= float(lines[12].split()[7])
    CShdr['pre_time']= float(lines[13].split()[2])
    if len(lines[13].split()) == 9:
        hdr['record_time'] = lines[13].split()[6]+" "+lines[13].split()[7]
        rt = 'NEW'
    else:
        hdr['record_time']=np.nan
        rt = 'OLD'
    CShdr['network'] = lines[14].split()[0]
    hdr['network'] = CShdr['network']
 #  ----------------------------------------------------------------------   
    # CSMNC files have directions listed as NS, EW, or UD,
    dir_string = CShdr['comp']
    if dir_string in ['NS']:
        hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                          is_acceleration=True,
                                          is_vertical=False,
                                          is_north=True)
    elif dir_string in ['EW']:
        hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                          is_acceleration=True,
                                          is_vertical=False,
                                          is_north=False)
    elif dir_string in ['UD']:
        hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                          is_acceleration=True,
                                          is_vertical=True,
                                          is_north=False)
    else:
        raise Exception('Could not parse direction %s' )
    #  ------------------------------------------------------------ 
    # The CSMNC data logger adds a 20s time delay CShdr['Pre_Time']
    # "new" format have recorder time
    if rt == 'NEW':
        pret = CShdr['pre_time']
        timestr = hdr['record_time']
        sttime = datetime.strptime(timestr, TIMEFMT) #- timedelta(seconds=pret)
        # Shift the time to UTC (Chinese time is 8 hours ahead)
        sttime = sttime - timedelta(seconds=8 * 3600.)
        print(sttime)
        hdr['starttime'] = sttime
    else:
        hdr['starttime'] =  datetime.strptime('9999-01-01 00:00:00.00', TIMEFMT)
    #  -----------------------------------------------------------------  
        
    # read in the data - there is a max of 8 columns per line
    # the code below handles the case when last line has
    # less than 8 columns
    if hdr['npts'] % COLS_PER_LINE != 0:
        nrows = int(np.floor(hdr['npts'] / COLS_PER_LINE))
        nrows2 = 1
    else:
        nrows = int(np.ceil(hdr['npts'] / COLS_PER_LINE))
        nrows2 = 0
    data = np.genfromtxt(filename, skip_header=TEXT_HDR_ROWS,
                         max_rows=nrows, filling_values=np.nan)
    data = data.flatten()
    if nrows2:
        skip_header = TEXT_HDR_ROWS + nrows
        data2 = np.genfromtxt(filename, skip_header=skip_header,
                              max_rows=nrows2, filling_values=np.nan)
        data = np.hstack((data, data2))
        nrows += nrows2
    #  -----------------------------------------------------------------    
    # fill out the rest of the standard dictionary
    standard['horizontal_orientation'] = np.nan
    standard['instrument_period'] = np.nan
    standard['instrument_damping'] = np.nan
    standard['process_time'] = ''
    standard['process_level'] = 'V1'
    standard['sensor_serial_number'] = ''
    standard['instrument'] = ''
    standard['comments'] = ''
    standard['structure_type'] = ''
    standard['corner_frequency'] = np.nan
    standard['units'] = 'acc'
    standard['source'] = SRC
    standard['source_format'] = 'CSMNC'
    hdr['coordinates'] = coordinates
    hdr['standard'] = standard
    #hdr['CSMNC'] = CShdr
    # create a Trace from the data and metadata
    trace = Trace(data.copy(), Stats(hdr.copy()))

    # to match the max values in the headers,
    # we need to detrend/demean the data (??)
#    trace.detrend('linear')
#    trace.detrend('demean')

    stream = Stream(trace)
    return stream
test=read_CSMNC("testcsmncdata.dat")
print(test)
test.plot()
