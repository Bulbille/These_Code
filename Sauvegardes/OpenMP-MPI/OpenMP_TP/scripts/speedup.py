#!/usr/bin/python
# encoding=utf8

import sys
import re
import collections
import os

try:
    import matplotlib.pyplot as plt
    canPlot = True
except ImportError:
    canPlot = False

numThreads = [0,1,2,4,6,8]

def parse_res(filename):
    try:
        D = collections.defaultdict(dict)
        n = collections.defaultdict(int)

        f = open(filename,'r')
        for line in iter(f):

            # Populate 'D'
            try:
                key = re.search(' *(.+?) *:', line).group(1)
                val = re.search(':  (.*)', line).group(1)
                #print "-%s-" % key
                #print "-%s-" % val
                D[key][numThreads[n[key]]] = val
                n[key] += 1

            except AttributeError:
                pass

        f.close()

        #print sorted(D, key=D.get)

        D2 = {}
        D2["+CPU (%)"] = {}
        D2["Speedup"] = {}

        for k in D:
        #    print k, ': ', D[k]

            DD=D[k];

            if k == "Temps CPU":

                refTime=float(DD[0].replace(' sec.' ,''));
                for kk in DD:
                    x = float(DD[kk].replace(' sec.', ''))
                    per = (x-refTime)*100/refTime
                    D2["+CPU (%)"][kk] = '+%.2f %%' % per

            elif k == "Temps elapsed":

                refTime=float(DD[0].replace(' sec.' ,''));
                for kk in DD:
                    x = float(DD[kk].replace(' sec.', ''))
                    per = refTime/x
                    D2["Speedup"][kk] = '%.2f' % per

        #   else:
        #       val = DD.values()
        #       r = (val[1:] == val[:-1])
        #       print k, ": ", r

        return D, D2
    except IOError:
        return None, None

if len(sys.argv) != 2:
    print 'Usage: ./speedup.py <omp_tp.res>'
    sys.exit(1)

filename=sys.argv[1]

print  "== ", filename, " =="

D, D2 = parse_res(filename)
if (D is None or D2 is None):
    print "Erreur : Le fichier", filename, "ne peut pas être ouvert."
    sys.exit(2)


filenameRef = os.path.split(filename.replace('res', 'ref'))
filenameRef = os.path.join(filenameRef[0], '../parametres', filenameRef[1])
Dref, D2ref = parse_res(filenameRef)

#print
print "Temps CPU seq.", D["Temps CPU"][0]
print "Speedup:      ", D2["Speedup"]
#print "Surcout CPU:  ", D2["+CPU (%)"]

if canPlot:
    x = D2["Speedup"].keys()[1:]
    plotPerfectSpeedup, = plt.plot(x, x, color='k', linestyle='--')
    plotSpeedup, = plt.plot(x, D2["Speedup"].values()[1:], color='b', marker='o')
    if (Dref is None or D2ref is None):
        plotSpeedupRef = None
    else:
        plotSpeedupRef, = plt.plot(D2ref["Speedup"].keys()[1:], D2ref["Speedup"].values()[1:], color='r', linestyle='--', marker='x')
    plt.title(u"Accélération " + filename)
    plt.axis([1,8,1,8])
    plt.xticks(x)
    plt.xlabel("Nombre de threads")
    plt.ylabel(u"Accélération")
    if (plotSpeedupRef is not None):
        plt.legend([plotSpeedup, plotSpeedupRef, plotPerfectSpeedup], [u"Accélération OpenMP", u"Accélération OpenMP de référence", u"Accélération idéale"], loc=2)
    else:
        plt.legend([plotSpeedup, plotPerfectSpeedup], [u"Accélération OpenMP", u"Accélération idéale"], loc=2)
    plt.grid(True)
    plt.savefig(filename.replace('res', 'png'))
    plt.show()
else:
    print "Le module matplotlib ne semble pas disponible, la courbe d'acceleration n'a donc pas pu etre generee."
    print "Executer la commande 'module load python/2.7.10' puis relancer le script devrait resoudre le probleme."

print  "=="
