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
    recall_score, precision_score, auc, precision_recall_curve
from sklearn.preprocessing import StandardScaler


def train_model(clf, x_train, y_train, class_cv):
    x_train.index = list(range(x_train.shape[0]))
    y_train.index = list(range(x_train.shape[0]))
    probs_test = np.zeros((x_train.shape[0], 2), dtype=np.float)
    preds_test = np.zeros(x_train.shape[0], dtype=np.int)
    f1_train, f1_test = [], []
    precision_train, precision_test = [], []
    recall_train, recall_test = [], []
    auc_roc_train, auc_roc_test = [], []
    auc_prc_train, auc_prc_test = [], []
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
        _precision, _recall, _ = precision_recall_curve(
            sub_y_train, pred_y_train[:, 1] - pred_y_train[:, 0])
        auc_prc_train.append(auc(_recall, _precision))
        _precision, _recall, _ = precision_recall_curve(
            sub_y_test, pred_y_test[:, 1] - pred_y_test[:, 0])
        auc_prc_test.append(auc(_recall, _precision))

        probs_test[idx_test] = pred_y_test
        preds_test[idx_test] = sub_pred_y_test

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
    print('AUC of PRC: Mean_Train{0:0.5f}\nMean_Test{1:0.5f}\n'.format(
        np.mean(auc_prc_train), np.mean(auc_prc_test)))
    print('-' * 50)
    array_importance = np.array(feature_importance)
    mean_importance = np.mean(array_importance, axis=0)
    normal_importance = mean_importance / np.sum(mean_importance)
    print(mean_importance)
    print(normal_importance)

    df_prob = pd.DataFrame(probs_test, columns=['prob_0', 'prob_1'])
    df_pred = pd.DataFrame(preds_test, columns=['pred'])
    df_out = pd.concat([df_pred, df_prob], axis=1)

    return df_out


def predict_metric(file_input_predict, model_in):
    df_input_predict = pd.read_csv(file_input_predict, sep='\t')
    df_y_predict = df_input_predict['label']
    df_x_predict = df_input_predict.loc[:, use_col]
    pred_y_predict = model_in.predict_proba(df_x_predict)
    bool_pred_y_predict = np.argmax(pred_y_predict, axis=1)
    f1_predict = f1_score(df_y_predict, bool_pred_y_predict)
    precision_predict = \
        precision_score(df_y_predict, bool_pred_y_predict)
    recall_predict = recall_score(df_y_predict, bool_pred_y_predict)
    auc_roc_predict = roc_auc_score(
        df_y_predict, pred_y_predict[:, 1] - pred_y_predict[:, 0])
    _precision, _recall, _ = precision_recall_curve(
        df_y_predict, pred_y_predict[:, 1] - pred_y_predict[:, 0])
    auc_prc_predict = auc(_recall, _precision)

    print('F1 Score: ', f1_predict)
    print('Precision: ', precision_predict)
    print('Recall: ', recall_predict)
    print('AUC of ROC: ', auc_roc_predict)
    print('AUC of PRC: ', auc_prc_predict)

    return


if __name__ == '__main__':
    time_start = time()
    # parameters
    n_fold = 10
    skf = StratifiedKFold(n_splits=n_fold, shuffle=True)
    path_root = '/local/zy/PEI'

    # dataset
    path_label = \
        path_root + '/mid_data/training_label/label_interactions_V1'
    file_input = path_label + '/training_set.txt'
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

    # df_features = df_features.loc[
    #               :, ['DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
    #                   'H3K4me3_H3K27ac', 'score_ctcf_insulator']]

    # use_col = ['DHS_DHS', 'H3K4me3_DHS', 'DHS_H3K27ac',
    #            'H3K4me3_H3K27ac', 'score_ctcf_insulator']
    use_col = ['score_ctcf_insulator']
    # scale_scores = StandardScaler().fit_transform(df_features)
    df_x_train_pre, df_x_test_pre, df_y_train, df_y_test = train_test_split(
        df_features, df_label, test_size=0.1, random_state=331)
    df_x_train = df_x_train_pre.loc[:, use_col]
    df_x_test = df_x_test_pre.loc[:, use_col]
    # scaler = StandardScaler()
    # X_train = scaler.fit_transform(df_x_train)
    # X_test = scaler.transform(df_x_test)

    # lightGB
    params = {
        'objective': 'binary',
        'learning_rate': 0.01,
        'min_child_samples': 15,
        'max_depth': 5,
        'num_leaves': 35,
        'lambda_l1': 1,
        'boosting': 'gbdt',
        # 'boosting': 'dart',
        'n_estimators': 2000,
        'metric': 'binary',
        'feature_fraction': .75,
        'bagging_fraction': .85,
        'bagging_freq': 3,
        'seed': 99,
        'num_threads': 20,
        'verbose': -1,
    }
    # params = {
    #     'boosting': 'gbdt',
    #     'objective': 'binary',
    #     'max_depth': 5,
    #     'num_leaves': 36,
    # }
    model_lgb = lgb.LGBMClassifier(**params)
    # train and cross validation
    res_pred = train_model(model_lgb, df_x_train, df_y_train, skf)
    df_res_train = pd.concat([df_x_train, df_y_train, res_pred], axis=1)

    # train model
    model_all = lgb.LGBMClassifier(**params)
    model_all.fit(df_x_train, df_y_train, eval_set=[(df_x_test, df_y_test)],
                  early_stopping_rounds=500, verbose=False)

    # test
    pred_y = model_all.predict_proba(df_x_test)
    bool_pred_y = np.argmax(pred_y, axis=1)
    f1_test_independent = f1_score(df_y_test, bool_pred_y)
    precision_test_independent = precision_score(df_y_test, bool_pred_y)
    recall_test_independent = recall_score(df_y_test, bool_pred_y)
    auc_roc_test_independent = roc_auc_score(
        df_y_test, pred_y[:, 1] - pred_y[:, 0])
    precision, recall, _thresholds = precision_recall_curve(
        df_y_test, pred_y[:, 1] - pred_y[:, 0])
    auc_prc_test_independent = auc(recall, precision)
    print('F1 Score of test set: ', f1_test_independent)
    print('Precision of test set: ', precision_test_independent)
    print('Recall of test set: ', recall_test_independent)
    print('AUC of ROC of test set: ', auc_roc_test_independent)
    print('AUC of PRC of test set: ', auc_prc_test_independent)
    df_prob_test = pd.DataFrame(pred_y, columns=['prob_0', 'prob_1'])
    df_pred_test = pd.DataFrame(bool_pred_y, columns=['pred'])
    df_res_test = pd.concat(
        [df_x_test, df_y_test, df_pred_test, df_prob_test], axis=1)
    df_res = pd.concat([df_res_train, df_res_test])
    # file_res = path_label + '/result_presiction.txt'
    # df_res.to_csv(file_res, sep='\t', index=None)

    # prediction
    file_input_h1 = path_label + '/H1/test_set.txt'
    print('H1')
    predict_metric(file_input_h1, model_all)

    file_input_IMR90 = \
        path_label + '/IMR-90/test_set.txt'
    print('IMR-90')
    predict_metric(file_input_IMR90, model_all)

    time_end = time()
    print(time_end - time_start)
