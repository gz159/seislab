from obspy import read

msfn = read('data.mseed')
length = len(msfn)
for i in range(length):
    wmsfn = msfn[i].stats.network + '_' + msfn[i].stats.station + '_' + msfn[i].stats.channel + '.mseed'
    msfn[i].write(wmsfn, format='MSEED')