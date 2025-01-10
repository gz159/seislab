import os
from glob import glob
import re
from subprocess import call
from sys import stdout
from os.path  import exists, isdir
from shutil   import copy, rmtree, move

def RM(name):
    if exists(name):
        if isdir(name): rmtree(name)
        else:           os.unlink(name)
    return

def list_specific_files(directory, filename):
    file_pattern = os.path.join(directory,filename)
    files = glob(file_pattern)
    return files

# Using PZs file to create a temporary meta file for mseed2sac
def create_temp_mseed_metafile(pzs_file):

    metafile = open('tmp_config.meta','w')
    headline = "#Network|Station|Location|Channel|"
    headline += "Latitude|Longitude|Elevation|Depth|"
    headline += "Azimuth|Dip|SensorDescription|Scale|"
    headline += "ScaleFreq|ScaleUnits|SampleRate|StartTime|EndTime\n"
    metafile.write(headline)

    with open(pzs_file, 'r') as file:
        for line in file:
            if "NETWORK" in line:
                 network = line.strip().split()[3]
            if "STATION" in line:
                 station = line.strip().split()[3]
            if "LOCATION" in line:
                 location = line.strip().split()[3]
            if "CHANNEL" in line:
                 channel = line.strip().split()[3]
            if "START" in line:
                 start = line.strip().split()[3]
            if "END" in line:
                 end = line.strip().split()[3]
            if "LATITUDE" in line:
                 latitude = float(line.strip().split()[3])
            if "LONGITUDE" in line:
                longitude = float(line.strip().split()[3])
            if "ELEVATION" in line:
                 elevation = float(line.strip().split()[3])
            if "DEPTH" in line:
                 depth = float(line.strip().split()[3])
            if "DIP" in line:
                 dip = float(line.strip().split()[3])
            if "AZIMUTH" in line:
                 azimuth = float(line.strip().split()[3])
            if "SAMPLE RATE" in line:
                 sample_rate = float(line.strip().split()[4])
            if "INSTTYPE" in line:
                 insttype= line.strip().split()[3]
            if "SENSITIVITY" in line:
                 sensitivity = float(line.strip().split()[3])
                 scaleunit=re.sub(r'^\(|\)$', '', line.strip().split()[4].strip())

        formatted_line = '%s|%s|%s|%s|'
        formatted_line += '%.2f|%.2f|%.1f|%.1f|'
        formatted_line += '%.1f|%.1f|%s|%.1f|'
        formatted_line += '%.1f|%s|%.1f|%s|%s\n'
        formatted_line  = formatted_line%(network,station,location,channel,
                                        latitude,longitude,elevation,depth,
                                        azimuth,dip,insttype,sensitivity,
                                        0.1,scaleunit,sample_rate,start,end)
        metafile.write(formatted_line)

    metafile.close()
    return metafile.name

#-------------------------Main Program----------
from datetime import datetime

# earthquake information
eymd = '2025-01-07'
jday = datetime(int(eymd.strip().split('-')[0]),int(eymd.strip().split('-')[1]),
                int(eymd.strip().split('-')[2])).timetuple().tm_yday
etime = '01:05:16'
elat = '28.50'
elon = '87.45'
edepth = '10'
ename = 'xizang'

# mseed file directory
directory = 'mseed'
filename = '*.mseed'
mseed_flist = list_specific_files(directory,filename)
mseed_flist.sort()

# instrument response file directory
directory = 'sacpz'
filename = 'SAC_PZs*'
pzs_flist = list_specific_files(directory,filename)
pzs_flist.sort()

# directory to keep observed SAC and instrument response file for wphase inverion
cwd =  os.getcwd()
destdir = 'DATA_org'
if exists(destdir): RM(destdir)
os.makedirs(destdir)


# for each mseed file, find the corresponding PZs file and create a temporary meta file for mseed2sac
# then run mseed2sac to convert mseed to SAC format,
# and move the SAC and corresponding PZs file to the destination directory
for mseed_file in mseed_flist:
    for pzs_file in pzs_flist:
        if mseed_file.strip().split('/')[1].split('.')[0] in pzs_file.strip().split('/')[1]:
            tmp_meta_file = create_temp_mseed_metafile(pzs_file.strip())
            cmd = 'mseed2sac -m %s -E %s,%s,%s/%s/%s/%s/%s %s'
            cmd = cmd%(tmp_meta_file,eymd.strip().split('-')[0],jday,etime,elat,elon,edepth,ename,mseed_file.strip())
            call(cmd, shell=True, stdout=stdout)
            os.unlink(tmp_meta_file)
            for p in glob('*.SAC'):
                if len(p) < 1:
                    continue
                else:
                    try:
                        move(p,destdir)
                        copy(pzs_file.strip(),destdir)
                    except Exception as e:
                        print(f"Error:{e}")
