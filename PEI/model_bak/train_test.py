#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: train_test.py
# @time: 2020/3/11 18:16

from time import time
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import roc_auc_score, f1_score, \
    recall_score, precision_score
from sklearn.preprocessing import StandardScaler


def train_model(clf, x_train, y_train, class_cv):
    x_train.index = list(range(x_train.shape[0]))
    y_train.index = list(range(x_train.shape[0]))
    preds_train = np.zeros((x_train.shape[0], 2), dtype=np.float)
    f1_train, f1_test = [], []
    precision_train, precision_test = [], []
    recall_train, recall_test = [], []
    auc_roc_train, auc_roc_test = [], []
    feature_importance = []

    i = 1
    for idx_train, idx_test in class_cv.split(x_train, y_train):
        sub_x_train, sub_y_train = x_train.loc[idx_train, :], \
            y_train.loc[idx_train]
        sub_x_test, sub_y_test = x_train.loc[idx_test, :], \
            y_train.loc[idx_test]

        clf.fit(sub_x_train, sub_y_train, eval_set=[(sub_x_test, sub_y_test)],
                early_stopping_rounds=500, verbose=False)

        pred_y_train = clf.predict_proba(sub_x_train)
        pred_y_test = clf.predict_proba(sub_x_test)
        feature_importance.append(clf.feature_importances_)
        print(clf.feature_importances_)
        sub_pred_y_train = np.argmax(pred_y_train, axis=1)
        sub_pred_y_test = np.argmax(pred_y_test, axis=1)
        f1_train.append(f1_score(sub_y_train, sub_pred_y_train))
        f1_test.append(f1_score(sub_y_test, sub_pred_y_test))
        precision_train.append(precision_score(sub_y_train, sub_pred_y_train))
        precision_test.append(precision_score(sub_y_test, sub_pred_y_test))
        recall_train.append(recall_score(sub_y_train, sub_pred_y_train))
        recall_test.append(recall_score(sub_y_test, sub_pred_y_test))
        auc_roc_train.append(roc_auc_score(
            sub_y_train, pred_y_train[:, 1] - pred_y_train[:, 0]))
        auc_roc_test.append(roc_auc_score(
            sub_y_test, pred_y_test[:, 1] - pred_y_test[:, 0]))

        preds_train[idx_test] = clf.predict_proba(sub_x_test)[:]

        print('{0}|F1 Score: Train {1:0.7f} Validation {2:0.7f}'.format(
            i, f1_train[-1], f1_test[-1]))
        print('-' * 50)
        i += 1

    print('F1 Score: Mean_Train{0:0.5f}\nMean_Test{1:0.5f}\n'.format(
        np.mean(f1_train), np.mean(f1_test)))
    print('-' * 50)
    print('Precision: Mean_Train{0:0.5f}\nMean_Test{1:0.5f}\n'.format(
        np.mean(precision_train), np.mean(precision_test)))
    print('-' * 50)
    print('Recall: Mean_Train{0:0.5f}\nMean_Test{1:0.5f}\n'.format(
        np.mean(recall_train), np.mean(recall_test)))
    print('-' * 50)
    print('AUC of ROC: Mean_Train{0:0.5f}\nMean_Test{1:0.5f}\n'.format(
        np.mean(auc_roc_train), np.mean(auc_roc_test)))
    print('-' * 50)
    array_importance = np.array(feature_importance)
    mean_importance = np.mean(array_importance, axis=0)
    normal_importance = mean_importance / np.sum(mean_importance)
    print(mean_importance)
    print(normal_importance)

    return


if __name__ == '__main__':
    time_start = time()
    # parameters
    n_fold = 10
    skf = StratifiedKFold(n_splits=n_fold, shuffle=True)

    # dataset
    file_input = \
        '/local/zy/PEI/mid_data/cell_line/model_input/GM12878/input_file.txt'
    df_input = pd.read_csv(file_input, sep='\t')
    df_label = df_input['label']
    df_features = (df_input.iloc[:, 3:]).drop('label', axis=1)
    # df_features = df_features.loc[
    #               :, ['score_dhs_enhancer',
    #                   'score_h3k27ac_enhancer', 'pval_h3k27ac_enhancer',
    #                   'distance', 'score_dhs_promoter',
    #                   'score_h3k4me3_promoter', 'pval_h3k4me3_promoter',
    #                   'gene_expression',
    #                   'score_dhs_insulator', 'score_ctcf_insulator',
    #                   'pcHi-C_ng2019', '3DIV', 'Thurman']]
    # df_features = df_features.loc[
    #               :, ['score_h3k27ac_enhancer',
    #                   'distance', 'score_h3k4me3_promoter',
    #                   'gene_expression', 'score_ctcf_insulator',
    #                   'pcHi-C_ng2019', '3DIV', 'Thurman']]
    df_features = df_features.loc[
                  :, ['score_h3k27ac_enhancer',
                      'distance', 'score_h3k4me3_promoter',
                      'gene_expression', 'score_ctcf_insulator',
                      'pcHi-C_ng2019', '3DIV', 'Thurman']]
    scale_scores = StandardScaler().fit_transform(df_features)
    df_x_train, df_x_test, df_y_train, df_y_test = train_test_split(
        df_features, df_label, test_size=0.1, random_state=331)
    scaler = StandardScaler()
    X_train = scaler.fit_transform(df_x_train)
    X_test = scaler.transform(df_x_test)

    # lightGB
    params = {
        'objective': 'binary',
        'learning_rate': 0.01,
        'min_child_samples': 15,
        'max_depth': 6,
        'num_leaves': 75,
        'lambda_l1': 1,
        'boosting': 'gbdt',
        # 'boosting': 'dart',
        'n_estimators': 2000,
        'metric': 'binary',
        'feature_fraction': .75,
        'bagging_fraction': .85,
        'bagging_freq': 3,
        'seed': 99,
        'num_threads': 40,
        'verbose': -1,
    }
    # params = {
    #     'boosting': 'gbdt',
    #     'objective': 'binary',
    #     'max_depth': 5,
    #     'num_leaves': 36,
    # }
    model_lgb = lgb.LGBMClassifier(**params)
    train_model(model_lgb, df_x_train, df_y_train, skf)

    time_end = time()
    print(time_end - time_start)
