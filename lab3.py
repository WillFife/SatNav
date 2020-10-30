"""
Main lab3 script. This is the code to produce the results wanted
in the "Analysis" portion of the lab procedures.
"""

from satloc import *
import matplotlib.pyplot as plt


# load the days worth of data
heads = ['RRT_week','RRT_sec','ORT_week',
         'ORT_sec','ORT_fracsec','dopp',
         'carrier_phase','pseudorange','C/N0',
         'validflag','error','status','type','svid'] 
channel = pd.read_csv(DIR + 'channel.csv', names=heads)

# load txau receiver position
txau = np.array([-743773.730, -5460645.972, 3200347.862]).reshape((3,1)) * u.meter

# slice to get every 30 seconds
data_30s = channel.loc[0::30, :]

# grab the GPS L1 CA
gps_L1_CA = data_30s[data_30s['type']==0]

# grab the gpsWeek
gpsWeek = int(gps_L1_CA.loc[0, 'ORT_week'])

data_2hrs = data_30s.loc[0:5200, :]

"""
# loop through each gps second to get the elevation
N_rows = len(gps_L1_CA)
t      = gps_L1_CA.loc[:, 'ORT_sec'].to_numpy(dtype=np.float64)
svids  = gps_L1_CA.loc[:, 'svid'].to_numpy()
c_n0   = gps_L1_CA['C/N0'].to_numpy()
elevs  = np.zeros(N_rows) * u.degree
sec_n  = t[0]
for ii in range(N_rows):
    gpsSec = t[ii]
    svid   = svids[ii]
    if t[ii] % 7200 == 0:
        sec_n = t[ii]
        r, el = channel2navsol(gpsWeek, int(gpsSec), svid, rx_ecef=txau, createf=True)
        elevs[ii] = el.to(u.degree)
    else:
        r, el = channel2navsol(gpsWeek, int(gpsSec), svid, sec_n=int(sec_n), rx_ecef=txau)
        elevs[ii] = el.to(u.degree)

plt.figure()
plt.title('C/N0 vs elevation')
plt.xlabel('elevation {deg}')
plt.ylabel('C/N0 {dB-Hz}')
plt.grid()
plt.plot(elevs, c_n0, color='black', marker='o', linestyle='')
plt.show()
"""