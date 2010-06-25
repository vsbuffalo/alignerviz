"""
aviz - AlignerVIZ
Author: Vince Buffalo <vsbuffaloAAAA@gmail.com> (with the poly-A tail removed)
Copyright (c) 2010 The Regents of University of California, Davis.
All rights reserved.
Based on the paper by Vlachos, Taneri, Keogh, and Yu, 2007
"""

USAGE = """%prog [-p graph_out.pdf] file.fasta
Fastest on a few relatively short (>1000bp) sequences.
"""

from optparse import OptionParser, OptionGroup
import os.path

import_msg = "Package '%s' was not found; please install it."
try:
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except ImportError:
    print import_msg % 'matplotlib'
try:
    import numpy as np
except ImportError:
    print import_msg % 'numpy'

from operator import itemgetter

V = {'A': (0, 1), 
     'T': (1, 0),
     'C': (0, -1),
     'G': (-1, 0)}

def basepair_vector(base):
    return V[base]

def trajectorize(seq):
    trajectory = []
    for base in seq:
        trajectory.append(basepair_vector(base))
    return trajectory

def render(trajectory):
    x, y, z = [0], [0], [0]
    for i, points in enumerate(trajectory):
        tx, ty = points
        x.append(x[-1] + tx)
        y.append(y[-1] + ty)
        z.append(i + 1)  # increment is a visual fix for first base

    points = (x, y, z)
    return points

def plot(seq_file, alpha, black, lwidth=1, show=True):
    mpl.rcParams['legend.fontsize'] = 10
    seqs = dict()
    file_obj = open(seq_file)
    i = 0
    for line in file_obj:
        if line.startswith('>'):
            # grab header, use as key
            line = line.strip()
            header = line[1:]

            # grab sequence (next in iterator) and trajectorizesequence 
            line = file_obj.next()
            line = line.strip()
            t = trajectorize(line)
            seqs[header] = t
        i += 1

    fig = plt.figure()
    ax = Axes3D(fig)

    for seq in seqs:
        x, y, z = render(seqs[seq])
        if black:
            ax.plot(x, y, z, label=seq, alpha=alpha, color='black')
        else:
            ax.plot(x, y, z, label=seq, alpha=alpha, linewidth=lwidth)

    ax.legend()
    if show:
        plt.show()
    return (plt, seqs)
    
if __name__ == "__main__":
    ## option parsing
    parser = OptionParser(usage=USAGE)
    parser.add_option("-p", "--plotpdf", dest="plotpdf", metavar="FILE", 
                      help="Output graph to PDF file specified (default: off)",
                      default=None)
    parser.add_option("-a", "--alpha", dest="alpha", default=1,
                      help="Set alpha (default: 1)")
    parser.add_option("-b", "--black", dest="black", default=False, action="store_true",
                      help="Force black color (default: off)")
    parser.add_option("-l", "--linewidth", dest="linewidth", default=1, 
                      help="Line width (default: 1)")
    (options, args) = parser.parse_args()

    if os.path.exists(args[0]):
        show = True if options.plotpdf is None else False
        plt, seqs = plot(args[0], black=options.black, lwidth=int(options.linewidth),
                         alpha=float(options.alpha), show=show)
        if options.plotpdf is not None:
            plt.savefig(options.plotpdf, format='pdf', facecolor='black')
    else:
        parser.error("file specified does not exist")
    
