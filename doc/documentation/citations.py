# creates: citations.png
import os
import datetime

import matplotlib.pyplot as plt


months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def f(filename):
    papers = {}
    lines = open(filename).readlines()
    n = 0
    dois = set()
    while n < len(lines):
        line = lines[n]
        tag = line[:2]
        if tag == 'TI':
            ntitle = n
            y = None
            m = 1
            d = 15
        elif tag == 'SO':
            title = ' '.join(lines[i][3:-1] for i in range(ntitle, n))
        elif tag == 'DI':
            doi = line[3:-1]
        elif tag == 'PY':
            y = int(line.split()[1])
        elif tag == 'PD':
            for w in line.split()[1:]:
                if w[0].isdigit():
                    w = int(w)
                    if w < 100:
                        d = w
                    else:
                        y = w
                else:
                    if '-' in w:
                        w = w.split('-')[-1]
                    m = months.index(w) + 1
        elif tag == '\n' and y is not None:
            date = datetime.date(y, m, d)
            if doi not in dois:
                dois.add(doi)
                papers[doi] = (date, title)
        n += 1

    return papers


# The papers here are:
label_bib = {
    'gpaw1':
        'Mortensen et al., Phys. Rev. B (2005)',
        # https://doi.org/10.1103/PhysRevB.71.035109
    'gpaw2':
        'Enkovaara et al., J. Phys.: Condens. Matter (2010)',
        # https://doi.org/10.1088/0953-8984/22/25/253202
    'lcao':
        'Larsen et al., Phys. Rev. B (2009)',
        # https://doi.org/10.1103/PhysRevB.80.195112
    'tddft':
        'Walter et al., J. Chem. Phys. (2008)',
        # https://doi.org/10.1063/1.2943138
    'response':
        'Yan et al., Phys. Rev. B (2011)',
        # https://doi.org/10.1103/PhysRevB.83.245122
}

plt.figure(figsize=(8, 4))
total = {}
# for bib in ['gpaw1', 'tddft', 'lcao', 'gpaw2', 'response']:
for bib in ['gpaw1', 'gpaw2']:
    papers = {}
    for line in open(bib + '.txt'):
        date, doi, title = line.split(' ', 2)
        papers[doi] = (datetime.date(*[int(x) for x in date.split('-')]),
                       title.strip())
    if os.path.isfile(bib + '.bib'):
        papers.update(f(bib + '.bib'))
    papers = sorted((papers[doi][0], doi, papers[doi][1]) for doi in papers)
    plt.plot([paper[0] for paper in papers], range(1, len(papers) + 1),
             '-o', label=label_bib[bib])
    with open(bib + '.txt', 'w') as fd:
        for date, doi, title in papers:
            fd.write('%d-%02d-%02d %s %s\n' % (date.year, date.month, date.day,
                                               doi, title))
            # assert '"' not in title, title
            total[doi] = (date, title)
    x = dict([(p[1], 0) for p in papers])
    print((bib, len(papers), len(x), len(total)))


allpapers = sorted((paper[0], doi, paper[1]) for doi, paper in total.items())
plt.plot([paper[0] for paper in allpapers], range(1, len(allpapers) + 1),
         '-o', label='Total')

if 0:
    with open('citations.csv', 'w') as fd:
        n = len(allpapers)
        for date, doi, title in allpapers[::-1]:
            fd.write('%d,":doi:`%s <%s>`"\n' % (n, title, doi))
            n -= 1

plt.xlabel('date')
plt.ylabel('number of citations')
plt.legend(loc='upper left')
plt.savefig('citations.png')
