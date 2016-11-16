import sqlite3
from itertools import chain, combinations
import numpy as np
from multiprocessing import Pool
from sklearn.svm import SVC
from scipy import stats
from sklearn import metrics
from sklearn.model_selection import cross_val_predict
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier

ORGANISMS = ['ce', 'dm', 'mm', 'sc']
CATEGORIES = ['bp', 'cc', 'mf']
CATEGORIES_SUBSETS = list(chain(*map(lambda i: combinations(CATEGORIES, i),
                                     range(1, len(CATEGORIES) + 1))))
GO_THRESHOLD = PPI_THRESHOLD = 3
SCORE_THRESHOLD = 0.500

NUM_SAMPLES = 11
NUM_PROCESSES = 11

K_FOLD = 10

METRICS_FILE_HEADER = [
    'Organism', 'Categories', 'GO Threshold', 'PPI Threshold', 'Score Type',
    'Classifier',
    'Accuracy', '(TPR*TNR)^1/2',
    'Precision A', 'Recall A', 'F1-score A', 'Support A',
    'Precision P', 'Recall P', 'F1-score P', 'Support P',
    'TN', 'FP', 'FN', 'TP', 'TPR', 'TNR'
]


def main(metrics_file):
    conn = sqlite3.connect('db.sqlite')
    c = conn.cursor()

    for org in ORGANISMS:
        print 'Organism %s.' % org

        q = c.execute("SELECT DISTINCT gene_id "
                      "FROM gene_class WHERE org = ?", (org,))
        gene_list = [str(row[0]) for row in q]

        class_value = {gene: {} for gene in gene_list}
        q = c.execute("SELECT gene_id, class "
                      "FROM gene_class WHERE org = ?", (org,))
        for row in q:
            class_value[row[0]] = row[1]

        q = c.execute("SELECT prot_id FROM gene_prot_%s "
                      "GROUP BY prot_id HAVING COUNT(DISTINCT gene_id) >= ?"
                      % org, (PPI_THRESHOLD,))
        prot_list = [str(row[0]) for row in q]

        ppis = {gene: {} for gene in gene_list}
        q = c.execute("SELECT gene_id, prot_id, score FROM gene_prot_%s" % org)
        for row in q:
            ppis[row[0]][row[1]] = float(row[2])/1000

        ppi_values = {gene: [] for gene in gene_list}
        for gene in gene_list:
            for p in prot_list:
                try:
                    ppi_values[gene].append(ppis[gene][p])
                except KeyError:
                    ppi_values[gene].append(0)

        # for cats in CATEGORIES_SUBSETS:
        for cats in [('bp', 'cc', 'mf')]:  # only 1 - all the cats
            q = c.execute("SELECT go_id FROM gene_go WHERE org = ? "
                          "AND (cat = '" + ("' OR cat = '".join(cats)) + "') "
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

            data = np.array([go_values[gene] + ppi_values[gene]
                             for gene in gene_list])
            data_binary = np.copy(data)
            data_binary[data_binary >= SCORE_THRESHOLD] = 1
            data_binary[data_binary < SCORE_THRESHOLD] = 0
            target = np.array([class_value[gene] for gene in gene_list])

            dataset_name = org + '-' + '+'.join(cats)
            print dataset_name

            dataset_attrs = [org, '+'.join(cats), GO_THRESHOLD, PPI_THRESHOLD]

            print 'SVC'
            clf = SVC(kernel='linear', C=1)
            attr_list = dataset_attrs + ['real', 'SVC linear C=1']
            pred = cross_val_predict(clf, data, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            attr_list = dataset_attrs + ['bin', 'SVC linear C=1']
            pred = cross_val_predict(clf, data_binary, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            print 'GNB'
            clf = GaussianNB()
            attr_list = dataset_attrs + ['real', 'GaussianNB']
            pred = cross_val_predict(clf, data, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            attr_list = dataset_attrs + ['bin', 'GaussianNB']
            pred = cross_val_predict(clf, data_binary, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            print 'KNN'
            clf = KNeighborsClassifier(n_neighbors=5, n_jobs=-1)
            attr_list = dataset_attrs + ['real', 'KNN 5']
            pred = cross_val_predict(clf, data, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            attr_list = dataset_attrs + ['bin', 'KNN 5']
            pred = cross_val_predict(clf, data_binary, target, cv=K_FOLD)
            writeMetrics(metrics_file, target, pred, attr_list)

            print 'ESB'
            attr_list = dataset_attrs + ['real',
                                         'Ensemble S=%d SVC linear C=1'
                                         % NUM_SAMPLES]
            pool = Pool(NUM_PROCESSES)
            predictions = pool.map(predictstar,
                                   ((data, target)
                                    for i in xrange(NUM_SAMPLES)))
            pred = stats.mode(np.array(predictions))[0][0]
            writeMetrics(metrics_file, target, pred, dataset_attrs)


def writeMetrics(f, target, prediction, attr_list):
    accuracy = metrics.accuracy_score(target, prediction)
    # print "\tAccuracy: %0.2f" % accuracy
    # print metrics.classification_report(target, prediction)
    cm = metrics.confusion_matrix(target, prediction)
    TN, FP, FN, TP = cm[0][0], cm[0][1], cm[1][0], cm[1][1]
    TPR, TNR = TP/float(TP+FN), TN/float(TN+FP)
    geometric_mean = ((TPR * TNR) ** (1.0/2.0))
    precision, recall, f1_score, support = (
        metrics.precision_recall_fscore_support(target, prediction))
    attr_list += [
        accuracy, geometric_mean,
        precision[0], recall[0], f1_score[0], support[0],
        precision[1], recall[1], f1_score[1], support[1],
        TN, FP, FN, TP, TPR, TNR
    ]
    f.write(','.join(map(converter, attr_list)) + '\n')


def converter(x):
    if isinstance(x, int):
        return '%d' % x
    elif isinstance(x, float):
        return '%0.2f' % x
    else:
        return x


def predictstar(args):
    return predict(*args)


def predict(data, target):
    clf = SVC(kernel='linear', C=1)
    np.random.seed()
    sample = np.random.rand(*data.shape)
    sample[sample >= (1-data)] = 1
    sample[sample < (1-data)] = 0
    prediction = cross_val_predict(clf, sample, target, cv=K_FOLD)
    return prediction

if __name__ == '__main__':
    with open('metrics_output.csv', 'w') as f:
        f.write(','.join(METRICS_FILE_HEADER) + '\n')
        # for GO_THRESHOLD in [3, 10, 50]:
        #     for PPI_THRESHOLD in [3, 10, 50]:
        main(f)
