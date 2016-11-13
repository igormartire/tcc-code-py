import sqlite3
import os
from itertools import chain, combinations

ORGANISMS = ['ce', 'dm', 'mm', 'sc']
CATEGORIES = ['bp', 'cc', 'mf']
CATEGORIES_SUBSETS = list(chain(*map(lambda i: combinations(CATEGORIES, i),
                                     range(1, len(CATEGORIES) + 1))))
GO_THRESHOLD = PPI_THRESHOLD = 3
SCORE_THRESHOLD = 400

OUTPUT_DIR = 'output_binary_score'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

conn = sqlite3.connect('db.sqlite')
c = conn.cursor()

for org in ORGANISMS:
    print 'Preparing to create ARFF files for organism %s.' % org

    q = c.execute("SELECT DISTINCT gene_id "
                  "FROM gene_class WHERE org = ?", (org,))
    gene_list = [str(row[0]) for row in q]

    class_value = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, class "
                  "FROM gene_class WHERE org = ?", (org,))
    for row in q:
        class_value[row[0]] = row[1]

    q = c.execute("SELECT prot_id FROM gene_prot_%s "
                  "GROUP BY prot_id HAVING COUNT(DISTINCT gene_id) >= ?" % org,
                  (PPI_THRESHOLD,))
    prot_list = [str(row[0]) for row in q]

    ppis = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, prot_id, score FROM gene_prot_%s" % org)
    for row in q:
        ppis[row[0]][row[1]] = '0.' + row[2]

    ppi_values = {gene: [] for gene in gene_list}
    for gene in gene_list:
        for p in prot_list:
            if p in ppis[gene].keys():
                score = int(ppis[gene][p])
                ppi_values[gene].append('1' if score >= SCORE_THRESHOLD
                                        else '0')
            else:
                ppi_values[gene].append('0')

    for cats in CATEGORIES_SUBSETS:
        print '\tPreparing to create ARFF files for category(ies) %s.' \
            % str(cats)

        q = c.execute("SELECT go_id FROM gene_go WHERE org = ? AND "
                      "(cat = '" + ("' OR cat = '".join(cats)) + "') "
                      "GROUP BY go_id HAVING COUNT(DISTINCT gene_id) >= ?",
                      (org, GO_THRESHOLD))
        go_list = [str(row[0]) for row in q]

        gos = {gene: [] for gene in gene_list}
        q = c.execute("SELECT gene_id, go_id "
                      "FROM gene_go WHERE org = ?", (org,))
        for row in q:
            gos[row[0]].append(row[1])

        go_values = {gene: [] for gene in gene_list}
        for gene in gene_list:
            for go in go_list:
                go_values[gene].append('1' if go in gos[gene] else '0')

        file_name = 'output/' + org + '-' + '+'.join(cats) + '.arff'
        print '\tCreating ARFF file: %s.' % file_name
        with open(file_name, 'w') as f:
            f.write('@relation GO_PPI_CE\n')
            f.write('\n')
            # f.write('@attribute gene string\n')
            for go in go_list:
                f.write('@attribute ' + go + ' {0,1}\n')
            for prot in prot_list:
                f.write('@attribute ' + prot + ' numeric\n')
            f.write('@attribute class {a,p}\n')
            f.write('\n')
            f.write('@data\n')
            for gene in gene_list:
                f.write(''  # gene + ','
                        + ','.join(go_values[gene] + ppi_values[gene])
                        + ',%s\n' % class_value[gene])
