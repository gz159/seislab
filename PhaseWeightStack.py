import numpy as np
import matplotlib.pyplot as plt
from obspy import read
from scipy.signal import hilbert

#You need to install the Python library stockwell first
from stockwell import st

def tf_PWS(stream, v, f1, f2):
    dt = stream[0].stats.delta
    b = stream[0].stats.sac.b
    t = np.arange(stream[0].stats.npts) * dt + b
    df = 1 / (stream[0].stats.sac.e-b)
    fn1 = int(f1/df); fn2 = int(f2/df)
    for i, tr in enumerate(stream):
        d = tr.data
        s = st.st(d, fn1, fn2)
        if i < 1:
            f = np.linspace(f1, f2, len(s[:, 0]))
            T, F = np.meshgrid(t, f)
            c = np.zeros_like(s)
        ph = s / abs(s) * np.exp(2j*np.pi*F*T)
        c += ph
    c /= len(stream)
    stc = stream.copy()
    stc.stack()
    ds = st.st(stc[0].data, fn1, fn2)
    pws = st.ist(ds*abs(c)**v, fn1, fn2)
    tr.data = pws
    return tr

def PWS(st, v, sm=False, sl=15):
    m = len(st)
    n = st[0].stats.npts
    dt = st[0].stats.delta
    t = np.arange(n) * dt
    c = np.zeros(n, dtype=complex)
    for i, tr in enumerate(st):
        h = hilbert(tr.data)
        c += h/abs(h)
    c = abs(c/m)
    if sm:
        operator = np.ones(sl) / sl
        c = np.convolve(c, operator, 'same')
    stc = st.copy()
    stc.stack()
    tr = stc[0]
    tr.data = tr.data*c**v
    return tr

# Gte the coherence of input stream.
def get_coh(st, v, sm=False, sl=20):
    m = len(st)
    n = st[0].stats.npts
    dt = st[0].stats.delta
    t = np.arange(n) * dt
    ht = np.zeros((m, n), dtype=complex)
    c = np.zeros(n)
    for i, tr in enumerate(st):
        ht[i] = hilbert(tr.data)
    pha = ht / abs(ht)
    for i in range(n):
        c[i] = abs( sum(pha[:, i]) )
    # Smooth the coherence if necessary.
    if sm:
        c = np.convolve(c/m, np.ones(sl)/sl, 'same') ** v
    else:
        c = ( c/m ) ** v
    return t, c

v = 2
st = read()
st.filter('bandpass', freqmin=2.5, freqmax=5, corners=4, zerophase=True)
d1 = st[0].data.copy()
st[0].data = d1 + np.random.randn(len(d1))*d1.max()*0.05
st[1].data = d1 + np.random.randn(len(d1))*d1.max()*0.10
st[2].data = d1 + np.random.randn(len(d1))*d1.max()*0.20

t, c = get_coh(st, v)

plt.figure(figsize=(8, 6))
plt.subplot(411)
plt.plot(t, st[0].data, lw=1, color='r', label='+%5 noise')
plt.legend(loc='upper right')
plt.subplot(412)
plt.plot(t, st[1].data, lw=1, color='g', label='+%10 noise')
plt.legend(loc='upper right')
plt.subplot(413)
plt.plot(t, st[2].data, lw=1, color='b', label='+%20 noise')
plt.legend(loc='upper right')
plt.subplot(414)
plt.plot(t, c, lw=1, color='k')
plt.tight_layout()
plt.show()