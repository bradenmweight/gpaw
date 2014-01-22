#!/usr/bin/env python

from optparse import OptionParser
import pylab as pl
#from glob import glob

p = OptionParser(usage='%prog [OPTION] FILE...',
                 description='plot timings from gpaw parallel timer.  '
                 'Typically one would run %prog timings.*.txt to plot '
                 'timings on all cores.  (Note: The plotting code is '
                 'rather hacky and ugly at the moment.)')
opts, fnames = p.parse_args()

fnames.sort()

class Entry:
    def __init__(self, name, t1, parent=None, childnumber=None):
        self.name = name
        self.t1 = t1
        self.parent = parent
        self.children = []
        self.childnumber = childnumber
        if parent is None:
            self.level = -1
        else:
            self.level = parent.level + 1

    def stop(self, t2):
        self.t2 = t2
        self.dt = self.t2 - self.t1

    def subentry(self, name, time):
        subentry = Entry(name, time, self, len(self.children))
        self.children.append(subentry)
        return subentry

    def iterate(self):
        for child1 in self.children:
            yield child1
            for child2 in child1.iterate():
                yield child2

    def normalize(self, start, stop):
        offset = start
        scale = 1.0 / (stop - start)
        self.transform(offset, scale)

    def transform(self, offset, scale):
        self.t1 = scale * (self.t1 - offset)
        self.t2 = scale * (self.t2 - offset)
        self.dt = self.t2 - self.t1
        for child in self.children:
            child.transform(offset, scale)

class EntryCollection:
    def __init__(self, entries):
        tstart = min([entry.t1 for entry in entries])
        tstop = max([entry.t2 for entry in entries])
        
        for entry in entries:
            assert entry.level == -1
            entry.normalize(tstart, tstop)
        self.entries = entries
        self.totals = self.get_totals()
    
    def get_totals(self):
        totals = {}
        for entry in self.entries:
            for child in entry.iterate():
                if child.name not in totals:
                    totals[child.name] = 0.0
                totals[child.name] += child.dt / len(self.entries)
        return totals

def get_timings(fname):
    root = Entry('root', 0.0)
    head = root
    for line in open(fname):
        line = line.strip()
        part1, part2, action = line.rsplit(' ', 2)
        tokens = part1.split(' ', 3)
        t = float(tokens[2])
        name = tokens[3]

        if action == 'started':
            head = head.subentry(name, t)
        else:
            assert action == 'stopped'
            assert name == head.name
            head.stop(t)
            head = head.parent
    assert head == root
    root.t1 = root.children[0].t1
    root.stop(root.children[-1].t2)
    return root

alltimings = []
metadata = None
    
for rank, fname in enumerate(fnames):
    if fname == 'timings.metadata.txt':
        assert metadata is None
        metadata = [line.strip() for line in open(fname)]
    else:
        alltimings.append(get_timings(fname))

if len(alltimings) == 0:
    p.error('no timings found')

if metadata is None:
    metadata = range(len(names))

entries = EntryCollection(alltimings)

ordered_names = []

fig = pl.figure(figsize=(12, 6))
fig2 = pl.figure()
ax = fig.add_subplot(111)
#nameax = fig.add_subplot(122)
nameax = fig2.add_subplot(111)
nameax.set_yticks([])
nameax.set_xticks([])
#nameax = [fig2.add_subplot(2, 2, i + 1) for i in range(4)]

styles = {}


for rank, rootnode in enumerate(alltimings):
    thecolors = 'bgrcmy'
    thehatches = ['', '//', r'\\', '.', 'o', '*']
    def getstyle(i):
        return thecolors[i % len(thecolors)], thehatches[i // len(thecolors)]

    nstyles_used = 0

    for child in rootnode.iterate():
        if child.name not in styles:
            if entries.totals[child.name] < 0.01:
                color, hatch = ('k', '')
            else:
                color, hatch = getstyle(nstyles_used)
                nstyles_used += 1
                ordered_names.append(child.name)
                #nameax[child.level].plot([], [],
                #                        color=color,
                #                        #hatch=hatch,
                #                        marker='s', ls='',
                #                        label=child.name)
            styles[child.name] = color, hatch
        
        color, hatch = styles[child.name]
        centerx = child.t1 + child.dt * 0.5
        centery = child.level + 0.5
        x = [child.t1, child.t2, child.t2, child.t1]
        y = [child.level, child.level, child.level + 1, child.level + 1]
        # hardcoded to max 5.  Graphs will overlap if larger.
        compression = 5.0
        y = [y1 / compression + rank - 0.25 for y1 in y]
        ax.fill(x, y, color=color,
                label='__nolegend__')
        ax.fill(x, y, color='None', hatch=hatch, edgecolor='k',
                alpha=0.7,
                label='__nolegend__')

ax.set_ylabel('rank')
ax.set_yticks(range(len(metadata)))
ax.set_yticklabels([txt.replace('=', '') for txt in metadata])

ax.set_xlabel('normalized time')

for i, name in enumerate(ordered_names):
    color, hatch = styles[name]
    x = [0, 0.8, 0.8, 0]
    y = [i, i, i + 0.8, i + .8]
    nameax.fill(x, y, color=color,
                label='__nolegend__')
    nameax.fill(x, y, color='None', hatch=hatch, edgecolor='k',
                alpha=0.7,
                label='__nolegend__')
    nameax.text(1.0, i + 0.4, name, va='center', ha='left', color='k')

nameax.axis(xmin=-1, xmax=len(ordered_names))
    
pl.show()
