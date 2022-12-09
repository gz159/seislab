def cross_spec_bf(ccfd, pos, f, fs, s1=1e-4, s2=5e-4, b1=0, b2=360):
    w = 2 * np.pi * f; tmp = -1j * w
    nf = len(ccfd[0]); fn = int(f*nf/fs)
    nb = 121; ns = 51
    baz = np.linspace(b1, b2, nb) * np.pi / 180
    bazz = baz - np.pi
    slow = np.linspace(s1, s2, ns)
    p = np.zeros((ns, nb), dtype=complex)
    for si, ss in enumerate(slow):
        for bj, bb in enumerate(bazz):
            shift = ss * (np.sin(bb)*pos[:, 0] + np.cos(bb)*pos[:, 1]) 
            p[si, bj] = sum(ccfd[:, fn]*np.exp(shift*tmp))
    return baz, slow, p


import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from obspy import read

def main():
    # Seismic data can be downloaded from "https://github.com/geophydog/Seismic_Data_Examples/Syn_Single_Seismic_Event"
    st = read('DATA/*.SAC')
    fd = []; coord = []
    for i, tr in enumerate(st):
        tr.data -= np.mean(tr.data)
        tr.detrend(); tr.data /= abs(tr.data).max()
        k = tr.stats.station
        x = tr.stats.sac.stlo
        y = tr.stats.sac.stla
        coord.append([x, y])
        tmp = np.fft.fft(tr.data)
        fd.append(list(tmp))
    fd = np.array(fd); coord = np.array(coord)
    n = len(fd[:, 0]); pos = []; ccfd = []
    for i in range(n-1):
        x1 = coord[i, 0]; y1 = coord[i, 1]
        for j in range(i+1, n):
            x2 = coord[j, 0]; y2 = coord[j, 1]
            pos.append([x2-x1, y2-y1])
            ccfd.append(fd[i]*np.conjugate(fd[j]))
    pos = np.array(pos) * 1e3
    ccfd = np.array(ccfd)
    s1 = 1e-4; s2 = 5e-4
    b1 = 0; b2 = 360
    T = np.linspace(5, 50, 10); fs = 1/ tr.stats.delta
   
    plt.figure(figsize=(18, 9.5))
    for ti, TT in enumerate(T):
        print('Scaning period %.1f s ...'%TT)
        f = 1 / TT
        baz, slow, p = cross_spec_bf(ccfd, pos, f, fs,
                                     s1=s1, s2=s2, b1=b1, b2=b2)
        PP = abs(p)
        bazz  = np.linspace(baz[0],  baz[-1], 361)
        sloww = np.linspace(slow[0], slow[-1], 201)
        inp = interpolate.interp2d(baz, slow, PP, kind='cubic')
        PP = inp(bazz, sloww)
        PP = PP / PP.max()
        ax1 = plt.subplot(2, 5, ti+1, projection='polar')
        plt.pcolormesh(bazz, sloww*1e3, PP, cmap='CMRmap')
        cbar = plt.colorbar(shrink=0.4, pad=0.2)
        cbar.set_label(r'Normalized $\theta-S$ Spectra', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        ax1.grid(color='gray', ls='--', lw=2)
        ax1.tick_params(axis='y', colors='gray', labelsize=12)
        ax1.set_theta_zero_location('N')
        ax1.set_theta_direction(-1)
        ax1.plot(baz, np.ones(len(baz))*slow.min()*1e3, c='k', lw=2)
        ax1.plot(baz, np.ones(len(baz))*slow.max()*1e3, c='k', lw=2)
        ax1.set_xlabel('Slowness (s/km)', fontsize=13)
        plt.xticks(fontsize=13); plt.yticks(fontsize=13)
        plt.title('%.1f s'%TT, fontsize=15, pad=15)
    plt.tight_layout() 
    #plt.savefig('Cross_Spectral_Beamforming.png')
    plt.show()

if __name__ == '__main__':
    main()