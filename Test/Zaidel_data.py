__author__ = 'JosephHughes'

import numpy as np

ncol = 200
dl = 500
dz = 5.
z0 = 20


ib = np.ones(ncol, np.int)
ib[0] = -1
ib[-1] = -1

z = z0
botm = np.empty(ncol, np.float)
icnt = 0
for idx in xrange(ncol):
    icnt += 1
    if icnt > 40:
        z -= dz
        icnt = 0
    botm[idx] = z

hsp1 = np.empty(ncol, np.float)
hsp1.fill(10.)
hsp1[0] = 23.

hsp2 = np.empty(ncol, np.float)
hsp2.fill(1.)
hsp2[0] = 23.

f = open('Zaidel_data.out', 'w')

for i in ib:
    f.write('{0:3d}'.format(i))
f.write('\n\n')

for v in botm:
    f.write('{0:5.0f}'.format(v))
f.write('\n\n')

for v in hsp1:
    f.write('{0:5.0f}'.format(v))
f.write('\n\n')

for v in hsp2:
    f.write('{0:5.0f}'.format(v))
f.write('\n\n')

f.close()