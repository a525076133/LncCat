import warnings

warnings.filterwarnings("ignore")
import sys
sys.path.append("/")
from multiprocessing import Pool, cpu_count
import numpy as np

from encode import pip_fs, read_fa

num = 10000000000000000
top=3

path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Human/lnc.test.human.B.fa"
path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Human/pct.test.human.B.fa"

# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Mouse/lnc.test.Mouse.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Mouse/pct.test.Mouse.fa"

# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Wheat/lnc.test.Wheat.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Wheat/pct.test.Wheat.fa"

# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/lnc.test.zebrafish.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/pct.test.zebrafish.fa"

# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Chicken/lnc.test.chicken.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Chicken/pct.test.chicken.fa"


# seqs_p_test = []
# tmp_list = list(read_fa(path_p_test).values())[0:num]
# seqs_p_test.extend(tmp_list)
#
# seqs_n_test = []
# tmp_list = list(read_fa(path_n_test).values())[0:num]
# seqs_n_test.extend(tmp_list)
#
# pool = Pool(cpu_count())
# p_data_test = pool.map(pip_fs, seqs_p_test)
# n_data_test = pool.map(pip_fs, seqs_n_test)
#
# np.save("human_lnc_data.npy",p_data_test)
# np.save("human_pct_data.npy",n_data_test)

p_data_test = np.load("human_lnc_data.npy")
n_data_test = np.load("human_pct_data.npy")




data_test = []
data_test.extend(p_data_test)
data_test.extend(n_data_test)
label_test = [1] * len(p_data_test) + [0] * len(n_data_test)

test_data = np.array(data_test)
test_label = np.array(label_test)

print("#######################")
print("#######################")
print("#######################")
print("#######################")
print("test data", test_data.shape)
print("test label", test_label.shape)

from catboost import CatBoostClassifier
clf = CatBoostClassifier()
from catboost import CatBoostClassifier
clf = CatBoostClassifier(
                         learning_rate = 0.01,
                         iterations=1000000000000000,
                         # iterations=1,
                         # n_estimators=3000,
                         depth = 10,
                         loss_function="Logloss",
                         early_stopping_rounds=5000,
                         # eval_metric="AUC",
                         # eval_metric="MCC",
                         eval_metric="Accuracy",
                         verbose=1000,
                         task_type="GPU",
                         devices='1:1'
                         # task_type="CPU",
                         )

clf.fit(X=train_data,y=train_label,
        eval_set=(val_data,val_label))

# clf.load_model("human_model")

import matplotlib.pyplot as plt
res = list(clf.feature_importances_)
fig = plt.figure()
plt.bar(x=list(range(len(res))),height=res)
plt.savefig("important.png")

res = list(clf.feature_importances_)
f = open("import.txt","w")
for i,x in enumerate(res):
    tmp = str(i)+","+str(x)+"\n"
    f.write(tmp)
f.close()

def met(model, val_data, val_label):
    pred_res = model.predict_proba(X=val_data)[:,1]
    # pred_res = pred_res[:,1]
    t = 0.5
    pred_label = []
    for x in pred_res:
        if x < t:
            pred_label.append(0)
        else:
            pred_label.append(1)

    from sklearn import metrics
    acc = metrics.accuracy_score(y_pred=pred_label, y_true=val_label)
    ap = metrics.average_precision_score(y_score=pred_res, y_true=val_label)
    mcc = metrics.matthews_corrcoef(y_pred=pred_label, y_true=val_label)
    f1 = metrics.f1_score(y_pred=pred_label, y_true=val_label)
    recall = metrics.recall_score(y_pred=pred_label, y_true=val_label)
    precision = metrics.precision_score(y_pred=pred_label, y_true=val_label)
    auc = metrics.roc_auc_score(y_true=val_label, y_score=pred_res)
    kappa = metrics.cohen_kappa_score(val_label, pred_label)
    tn, fp, fn, tp = metrics.confusion_matrix(y_true=val_label, y_pred=pred_label).ravel()
    specificity = tn / (tn + fp)

    print("acc",acc)
    print("ap",ap)
    print("mcc",mcc)
    print("f1",f1)
    print("recall",recall)
    print("precision",precision)
    print("auc",auc)
    print("kappa",kappa)
    print("specificity",specificity)
    print("tn",tn)
    print("tp",tp)
    print("fn",fn)
    print("fp",fp)
met(clf, test_data, test_label)