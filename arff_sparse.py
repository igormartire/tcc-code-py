import sqlite3
import os
from time import time
from itertools import chain, combinations

ORGANISMS = ['ce', 'dm', 'mm', 'sc']
CATEGORIES = ['bp', 'cc', 'mf']
CATEGORIES_SUBSETS = list(chain(*map(lambda i: combinations(CATEGORIES, i),
                                     range(1, len(CATEGORIES) + 1))))
GO_THRESHOLD = PPI_THRESHOLD = 3

OUTPUT_DIR = 'output_sparse'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

conn = sqlite3.connect('db.sqlite')
c = conn.cursor()

for org in ORGANISMS:
    print 'Preparing to create ARFF files for organism %s.' % org

    print '\tGetting the genes for this organism.'
    start = time()
    q = c.execute("SELECT DISTINCT gene_id "
                  "FROM gene_class WHERE org = ?", (org,))
    gene_list = [str(row[0]) for row in q]
    print '\tFinished in %s seconds.' % int(time() - start)

    print '\tGetting the class value for each gene.'
    start = time()
    class_value = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, class "
                  "FROM gene_class WHERE org = ?", (org,))
    for row in q:
        class_value[row[0]] = row[1]
    print '\tFinished in %s seconds.' % int(time() - start)

    print '\tGetting the proteins with which the genes interact.'
    start = time()
    q = c.execute("SELECT prot_id FROM gene_prot_%s "
                  "GROUP BY prot_id HAVING COUNT(DISTINCT gene_id) >= ?" % org,
                  (PPI_THRESHOLD,))
    prot_list = [str(row[0]) for row in q]
    print '\tFinished in %s seconds.' % int(time() - start)

    print '\tGetting the scores of each protein-protein interaction.'
    start = time()
    ppis = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, prot_id, score FROM gene_prot_%s" % org)
    for row in q:
        ppis[row[0]][row[1]] = '0.' + row[2]
    print '\tFinished in %s seconds.' % int(time() - start)

    print '\tPreparing all the scores to be written.'
    start = time()
    idx = {}
    for i, p in enumerate(prot_list):
        idx[p] = i
    ppi_values = {gene: [] for gene in gene_list}
    for gene in gene_list:
        for p in ppis[gene].keys():
            try:
                ppi_values[gene].append('%d %s' % (idx[p], ppis[gene][p]))
            except KeyError:
                pass  # gene interacts with protein eliminated by the threshold
        ppi_values[gene].sort(key=(lambda x: int(x.split(' ')[0])))
    print '\tFinished in %s seconds.' % int(time() - start)

    n_prot = len(prot_list)

    for cats in CATEGORIES_SUBSETS:
        print '\tPreparing to create ARFF files for category(ies) %s.' \
            % str(cats)

        q = c.execute("SELECT go_id FROM gene_go WHERE org = ? AND "
                      "(cat = '" + ("' OR cat = '".join(cats)) + "') "
                      "GROUP BY go_id HAVING COUNT(DISTINCT gene_id) >= ?",
                      (org, GO_THRESHOLD))
        go_list = [str(row[0]) for row in q]

        gos = {gene: [] for gene in gene_list}
        q = c.execute("SELECT DISTINCT gene_id, go_id "
                      "FROM gene_go WHERE org = ?", (org,))
        for row in q:
            gos[row[0]].append(row[1])

        idx = {}
        for i, go in enumerate(go_list):
            idx[go] = i + n_prot
        go_values = {gene: [] for gene in gene_list}
        for gene in gene_list:
            for go in gos[gene]:
                try:
                    go_values[gene].append('%d 1' % idx[go])
                except KeyError:
                    pass  # gene has go term eliminated by the threshold
            go_values[gene].sort(key=(lambda x: int(x.split(' ')[0])))

        n_go = len(go_list)

        file_name = OUTPUT_DIR + '/' + org + '-' + '+'.join(cats) + '.arff'
        print '\tCreating ARFF file: %s.' % file_name
        with open(file_name, 'w') as f:
            f.write('@relation GO_PPI_CE\n')
            f.write('\n')
            # f.write('@attribute gene string\n')
            for prot in prot_list:
                f.write('@attribute ' + prot + ' numeric\n')
            for go in go_list:
                f.write('@attribute ' + go + ' {0,1}\n')
            f.write('@attribute class {a,p}\n')
            f.write('\n')
            f.write('@data\n')
            for gene in gene_list:
                f.write('{'  # gene + ','
                        + ','.join(ppi_values[gene] + go_values[gene])
                        + ',%d %s' % (n_prot + n_go, class_value[gene])
                        + '}\n')
