import sys
from props import *

for line in sys.stdin:
    x,y = line.split()
    if y == "None": y = None
    sim2D = similarity(x, y)
    try:
        if sim2D > .2:
            print x, y, sim2D#, bandgap(y)
            #print sim2D
            #print "\n"
    except Exception as e:
        #print 0.0
        #print x, y, 0.0
        #print "\n"
        pass
        