import sqlite3
import os
from time import time
from itertools import chain, combinations
import numpy as np

ORGANISMS = ['ce', 'dm', 'mm', 'sc']
CATEGORIES = ['bp', 'cc', 'mf']
CATEGORIES_SUBSETS = list(chain(*map(lambda i: combinations(CATEGORIES, i),
                                     range(1, len(CATEGORIES) + 1))))
GO_THRESHOLD = PPI_THRESHOLD = 3
SCORE_THRESHOLD = 0.500

NUM_SAMPLES = 11

OUTPUT_DIR = 'output'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

conn = sqlite3.connect('db.sqlite')
c = conn.cursor()

for org in ORGANISMS:
    print 'Organism %s.' % org

    # print '\tGetting the genes for this organism.'
    # start = time()
    q = c.execute("SELECT DISTINCT gene_id "
                  "FROM gene_class WHERE org = ?", (org,))
    gene_list = [str(row[0]) for row in q]
    # print '\tFinished in %s seconds.' % int(time() - start)

    # print '\tGetting the class value for each gene.'
    # start = time()
    class_value = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, class "
                  "FROM gene_class WHERE org = ?", (org,))
    for row in q:
        class_value[row[0]] = 0 if row[1] == 'a' else 1
    # print '\tFinished in %s seconds.' % int(time() - start)

    # print '\tGetting the proteins with which the genes interact.'
    # start = time()
    q = c.execute("SELECT prot_id FROM gene_prot_%s "
                  "GROUP BY prot_id HAVING COUNT(DISTINCT gene_id) >= ?" % org,
                  (PPI_THRESHOLD,))
    prot_list = [str(row[0]) for row in q]
    # print '\tFinished in %s seconds.' % int(time() - start)

    # print '\tGetting the scores of each protein-protein interaction.'
    # start = time()
    ppis = {gene: {} for gene in gene_list}
    q = c.execute("SELECT gene_id, prot_id, score FROM gene_prot_%s" % org)
    for row in q:
        ppis[row[0]][row[1]] = float(row[2])/1000
    # print '\tFinished in %s seconds.' % int(time() - start)

    print '\tPreparing all the scores to be written.'
    start = time()
    ppi_values = {gene: [] for gene in gene_list}
    ppi_values_sample = [{gene: [] for gene in gene_list}
                         for i in range(NUM_SAMPLES)]
    for gene in gene_list:
        for p in prot_list:
            try:
                ppi_values[gene].append(ppis[gene][p])
            except KeyError:
                ppi_values[gene].append(0)
    print '\tFinished in %s seconds.' % int(time() - start)

    for cats in CATEGORIES_SUBSETS:
        # print '\tPreparing to create ARFF files for category(ies) %s.' \
        #     % str(cats)

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

        go_values = {gene: [] for gene in gene_list}
        for gene in gene_list:
            for go in go_list:
                go_values[gene].append(1 if go in gos[gene] else 0)

        dataset_name = org + '-' + '+'.join(cats)

        print '\tPreparing data.'
        start = time()
        data = np.array([go_values[gene] + ppi_values[gene] for gene in gene_list])
        data_binary = np.copy(data)
        data_binary[data_binary >= SCORE_THRESHOLD] = 1
        data_binary[data_binary < SCORE_THRESHOLD] = 0
        target = np.array([class_value[gene] for gene in gene_list])
        print '\tFinished in %s seconds.' % int(time() - start)
        start = time()
        from sklearn.model_selection import cross_val_score
        from sklearn.svm import SVC
        clf = SVC(kernel='linear', C=1)
        # scores = cross_val_score(clf, data, target, cv=5)
        # print("\tAccuracy SVC (%s):\t%0.2f (+/- %0.2f)"
        #       % (dataset_name, scores.mean(), scores.std() * 2))
        # scores = cross_val_score(clf, data_binary, target, cv=5)
        # print("\tAccuracy SVC b(%s):\t%0.2f (+/- %0.2f)"
        #       % (dataset_name, scores.mean(), scores.std() * 2))
        # from sklearn.naive_bayes import GaussianNB
        # gnb = GaussianNB()
        # scores = cross_val_score(gnb, data, target, cv=5)
        # print("\tAccuracy GNB (%s):\t%0.2f (+/- %0.2f)"
        #       % (dataset_name, scores.mean(), scores.std() * 2))
        # scores = cross_val_score(gnb, data_binary, target, cv=5)
        # print("\tAccuracy GNB b(%s):\t%0.2f (+/- %0.2f)"
        #       % (dataset_name, scores.mean(), scores.std() * 2))
        from scipy import stats
        from sklearn import metrics
        from sklearn.model_selection import cross_val_predict
        print '\tEnsemble method.'
        start = time()
        predicted = []
        for i in range(NUM_SAMPLES):
            sample = np.random.rand(*data.shape)
            sample[sample >= (1-data)] = 1
            sample[sample < (1-data)] = 0
            pred = cross_val_predict(clf, sample, target, cv=5)
            print '\t%s' % pred[:10]
            predicted.append(pred)
        m = stats.mode(np.array(predicted))
        score = metrics.accuracy_score(target, m[0][0])
        print("\tAccuracy ESB (%s):\t%0.2f" % (dataset_name, score))
        print '\tFinished in %s seconds.' % int(time() - start)


# >>> from sklearn.model_selection import cross_val_predict
# >>> predicted = cross_val_predict(clf, iris.data, iris.target, cv=10)
# >>> metrics.accuracy_score(iris.target, predicted)
# 0.966...
