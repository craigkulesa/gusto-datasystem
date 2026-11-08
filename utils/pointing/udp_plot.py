#%matplotlib inline
from GUSTOTelemetryData import GUSTOTelemetryData
import matplotlib.pyplot as plt
import glob

plt.rcParams['axes.formatter.useoffset'] = False

upath = '/Users/volkertolls/Projects/GUSTO/Data/Iridium/udps/'
upath = '/data/scratch/GUSTO/gusto-datasystem/utils/pointing/UDP/'

filter = 'udp_0[2-9]*.dat'
fss = sorted(glob.glob(upath+filter))


fig, ax1 = plt.subplots(1, 1, sharex=True, figsize=(12,8))


for fs in fss:
    gtd = GUSTOTelemetryData(fs)
    l, b, utime = gtd.getlbPosition()
    ax1.plot(l, b, '*', linewidth=1.5)

ax1.set_ylabel('Galactic Latitude (deg)')
ax1.set_xlabel('Galactic Longitude (deg)')

fig.tight_layout(pad=0.1)
#plt.savefig("plot.png", dpi=130.0, pad_inches=0.1, bbox_inches = 'tight')

plt.xlim([330,359])
plt.ylim([-3.2,2.2])
plt.gca().invert_xaxis()
plt.savefig('udp_test.png')
