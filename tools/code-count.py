#!/usr/bin/python

import os
import datetime
import pylab as pl


def count(dir, pattern):
    if not os.path.isdir(dir):
        return 0
    p = os.popen('wc -l `find %s -name %s` | tail -1' % (dir, pattern), 'r')
    return int(p.read().split()[0])


def polygon(x, y1, y2, *args, **kwargs):
    x = pl.concatenate((x, x[::-1]))
    y = pl.concatenate((y1, y2[::-1]))
    pl.fill(x, y, *args, **kwargs)


def plot_count(fname):
    # Load data
    date, libxc, c, code, test, doc = pl.load(fname, unpack=True)
    zero = pl.zeros_like(date)

    fig = pl.figure(1, figsize=(7.5, 5))
    ax = fig.add_subplot(111)
    polygon(date, c + code + test, c + code + test + doc,
             facecolor='m', label='Documentation')
    polygon(date, c + code, c + code + test,
             facecolor='y', label='Tests')
    polygon(date, c, c + code,
            facecolor='g', label='Python-code')
    polygon(date, zero, c,
            facecolor='r', label='C-code')

    months = pl.MonthLocator()
    months3 = pl.MonthLocator(interval=3)
    month_year_fmt = pl.DateFormatter("%b '%y")

    ax.xaxis.set_major_locator(months3)
    ax.xaxis.set_minor_locator(months)
    ax.xaxis.set_major_formatter(month_year_fmt)
    labels = ax.get_xticklabels()
    pl.setp(labels, rotation=30)
    pl.axis('tight')
    pl.legend(loc='upper left')
    pl.title('Number of lines')
    pl.savefig(fname.split('.')[0] + '.png')


if __name__ == '__main__':
    # Milestones:
    # Rev    1: 10/19/05 16:35:46  start revision log
    # Rev  383: 08/28/06 08:21:40  gridpaw -> gpaw
    # Rev  887: 07/11/07 10:33:37  libxc introduced
    # Rev 2050: 06/20/08 09:54:24  /doc in svn
    # Rev 5020: 09/28/09 14:43:03  /test starts moving to /gpaw/test
    
    # Check if stat file already exists, and stat from last checked day
    if os.path.isfile('stat.dat'):
        date1 = int(os.popen('tail -1 stat.dat','r').read().split()[0]) + 1
        stat = open('stat.dat', 'a')
    else:
        date1 = datetime.date(2005, 10, 20).toordinal()
        stat = open('stat.dat', 'w')
    
    date2 = datetime.date.today().toordinal()
    for datenum in range(date1, date2 + 1):
        datestr = datetime.date.fromordinal(datenum).isoformat()
        print datestr

        # Checkout of relevant gpaw folders
        svn = ('svn export --revision {%s} '
               'https://svn.fysik.dtu.dk/projects/gpaw/trunk' % datestr)
        e = os.system(svn + ' temp-gpaw > /dev/null')

        if e != 0:
            os.system('rm -rf temp-gpaw')
            continue

        # Remove gui:
        os.system('rm -rf temp-gpaw/gpaw/gui')

        libxc = count('temp-gpaw/c/libxc', '\*.[ch]')
        ch = count('temp-gpaw/c', '\*.[ch]') - libxc
        py = count('temp-gpaw/gridpaw', '\*.py')
        py += count('temp-gpaw/gpaw', '\*.py')
        test = count('temp-gpaw/test', '\*.py')
        test += count('temp-gpaw/gpaw/test', '\*.py')
        doc = count('temp-gpaw/doc', '\*.py -or -name \*.rst')

        # Clean up
        os.system('rm -rf temp-gpaw')

        # Dump data to stat file
        print >> stat, datenum, libxc, ch, py, test, doc

    stat.close()
    plot_count('stat.dat')

