from evtk.hl import pointsToVTK
import numpy as np

infile = open('cluster.txt', 'r')
n = int(infile.readline())
x = np.array([0.0]*n)
y = np.array([0.0]*n)
z = np.array([0.0]*n)    # it is a 2d structure, so left the z being 0.0
value = np.zeros(n)
for i in xrange(n) :
    line = infile.readline()
    text = line.split()
    x[i] = float(text[0])
    y[i] = float(text[1])
    value[i] = int(text[2])

x = x - np.mean(x)
y = y - np.mean(y)

pointsToVTK('dlca', x, y, z, data = {"value" : value})
infile.close()
