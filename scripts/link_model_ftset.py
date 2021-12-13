import pandas as pd
import json

#####
f_imp = 'src/classifier_metrics_20210821/feature_importance.tsv'
df_imp_raw = pd.read_csv(f_imp, sep='\t')
####

team_options= ['aklimate', 'CF', 'jadbio', 'subSCOPE', 'skgrid']
team_imp_options = ['aklimate', 'cloudforest', 'jadbio', 'subscope', 'skgrid']
cancer_list = ['ACC', 'BLCA', 'BRCA', 'CESC', 'COADREAD', 'ESCC', 'GEA', 'HNSC', 'KIRCKICH', 'KIRP', 'LGGGBM', 'LIHCCHOL', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'SARC', 'SKCM', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UVM']


def skgrid_get_more_model_info(cancer):
    '''skgrid in normal pipeline doesnt have enough info to pinpoint one ft selection
    method and classification so need to pull classification info

    Reason: skgrid needs info on classifier to match best model
    '''
    # Open file
    f_pred = 'src/classifier_metrics_20210821/top_performing_models_lte_100_features.tsv'
    df_pred = pd.read_csv(f_pred, sep='\t')

    # Get classifcation info
    skgrid_s1 = df_pred[df_pred['feature_list_method']=='skgrid']
    skgrid_s1 = skgrid_s1[skgrid_s1['cohort']==cancer].reset_index(drop=True)
    if skgrid_s1.shape[0]==1:
        selected_skgrid_model = skgrid_s1['model'][0]
        return(
            '## SKGRID ONLY. featureID and model info\n{}\n{}'.format(skgrid_s1['featureID'][0], selected_skgrid_model),
            selected_skgrid_model
        )
    else:
        return(
            'MULTIPLE TIED PERFORMING MODELS (N={})'.format(skgrid_s1.shape[0]),
            list(skgrid_s1['model'])
        )

skgrid_function_dict = {
    'ExtraTrees' : 'ExtraTreesClassifier',
    'LogisticRegression' : 'LogisticRegression', #nochange
    'RandomForest' : 'RandomForestClassifier',
    'AdaBoost' : 'AdaBoostClassifier',
    'BernoulliNB' : 'BernoulliNB', #nochange
    'DecisionTree' : 'DecisionTreeClassifier',
    'GaussianNB' : 'GaussianNB', #nochange
    'GaussianProcess' : 'GaussianProcessClassifier',
    'KNeighbors' : 'KNeighborsClassifier',
    'SGD' : 'SGDClassifier',
    'SVC' : 'SVC', #nochange
}

# If ties then resolve from past manual picking
previous_selected = {
    tuple(['skgrid_RandomForest(criterion=entropy,n_estimators=150)|skgrid_BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=150)|skgrid_BRCA.tsv_skgrid_BRCA_fbedeBIC_perplatformGEXP|BRCA.tsv_skgrid_BRCA_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=200)|skgrid_BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=200)|skgrid_BRCA.tsv_skgrid_BRCA_fbedeBIC_perplatformGEXP|BRCA.tsv_skgrid_BRCA_fbedeBIC_perplatformGEXP|2021-01-13|c'])
    : 'skgrid_RandomForest(criterion=entropy,n_estimators=150)|skgrid_BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|BRCA.tsv_skgrid_BRCA_fbedeBIC_combined|2021-01-13|c',

    tuple(['skgrid_LogisticRegression(C=100,max_iter=500,solver=newton-cg)|skgrid_BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=newton-cg)|skgrid_BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=10,max_iter=500,solver=newton-cg)|skgrid_BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|BLCA.tsv_skgrid_BLCA_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_LogisticRegression(C=0.01,max_iter=500,solver=newton-cg)|skgrid_HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=0.1,max_iter=500,solver=newton-cg)|skgrid_HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=0.01,max_iter=500,solver=newton-cg)|skgrid_HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|HNSC.tsv_skgrid_HNSC_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_SGD(alpha=0.001,loss=modified_huber,penalty=l1)|skgrid_KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_combined|KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_combined|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=modified_huber,penalty=l1)|skgrid_KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_perplatformGEXP|KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_perplatformGEXP|2021-01-13|c'])
    : 'skgrid_SGD(alpha=0.001,loss=modified_huber,penalty=l1)|skgrid_KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_combined|KIRCKICH.tsv_skgrid_KIRCKICH_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=newton-cg)|skgrid_LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|LGGGBM.tsv_skgrid_LGGGBM_rfe15_perplatformMETH|2021-01-13|c',
    tuple(['skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|PAAD.tsv_skgrid_PAAD_fbedeBIC_perplatformGEXP|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|PAAD.tsv_skgrid_PAAD_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_SVC(C=1.5,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_combined|SARC.tsv_skgrid_SARC_fbedeBIC_combined|2021-01-13|c',
 'skgrid_SVC(C=1.5,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_combined|SARC.tsv_skgrid_SARC_fbedeBIC_combined|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_combined|SARC.tsv_skgrid_SARC_fbedeBIC_combined|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|SARC.tsv_skgrid_SARC_fbedeBIC_perplatformGEXP|2021-01-13|c'])
    : 'skgrid_SVC(C=1,kernel=linear)|skgrid_SARC.tsv_skgrid_SARC_fbedeBIC_combined|SARC.tsv_skgrid_SARC_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_DecisionTree(criterion=entropy,max_depth=4,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=5,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=6,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=7,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=8,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=4,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=5,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=6,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=7,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=8,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianProcess(kernel=sklearn.gaussian_process.kernels.RationalQuadratic,max_iter_predict=1000)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c']) :
    'skgrid_DecisionTree(criterion=gini,max_depth=4,min_samples_split=3)|skgrid_SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|SKCM.tsv_skgrid_SKCM_fbedeBIC_perplatformMUTA|2021-01-13|c',
    tuple(['skgrid_LogisticRegression(C=0.01,max_iter=500,solver=lbfgs)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=0.01,max_iter=500,solver=newton-cg)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=lbfgs)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=newton-cg)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=0.01,max_iter=500,solver=lbfgs)|skgrid_TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|TGCT.tsv_skgrid_TGCT_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_GaussianNB(var_smoothing=1e-06)|skgrid_THYM.tsv_skgrid_THYM_fbedeBIC_perplatformMIR|THYM.tsv_skgrid_THYM_fbedeBIC_perplatformMIR|2021-01-13|c',
 'skgrid_SVC(C=0.2,kernel=linear)|skgrid_THYM.tsv_skgrid_THYM_fbedeBIC_combined|THYM.tsv_skgrid_THYM_fbedeBIC_combined|2021-01-13|c'])
    : 'skgrid_GaussianNB(var_smoothing=1e-06)|skgrid_THYM.tsv_skgrid_THYM_fbedeBIC_perplatformMIR|THYM.tsv_skgrid_THYM_fbedeBIC_perplatformMIR|2021-01-13|c',
    tuple(['skgrid_ExtraTrees(criterion=gini,n_estimators=128)|skgrid_UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=96)|skgrid_UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|2021-01-13|c'])
    : 'skgrid_ExtraTrees(criterion=gini,n_estimators=96)|skgrid_UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|UCEC.tsv_skgrid_UCEC_fbedeBIC_perplatformALL|2021-01-13|c',
    tuple(['skgrid_RandomForest(criterion=entropy,n_estimators=200)|skgrid_UVM.tsv_skgrid_UVM_fbedeBIC_combined|UVM.tsv_skgrid_UVM_fbedeBIC_combined|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=200)|skgrid_UVM.tsv_skgrid_UVM_fbedeBIC_perplatformGEXP|UVM.tsv_skgrid_UVM_fbedeBIC_perplatformGEXP|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=60)|skgrid_UVM.tsv_skgrid_UVM_fbedeBIC_combined|UVM.tsv_skgrid_UVM_fbedeBIC_combined|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=60)|skgrid_UVM.tsv_skgrid_UVM_fbedeBIC_perplatformGEXP|UVM.tsv_skgrid_UVM_fbedeBIC_perplatformGEXP|2021-01-13|c'])
    : 'skgrid_RandomForest(criterion=entropy,n_estimators=60)|skgrid_UVM.tsv_skgrid_UVM_fbedeBIC_combined|UVM.tsv_skgrid_UVM_fbedeBIC_combined|2021-01-13|c',
    tuple(['skgrid_AdaBoost(learning_rate=0.01,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.01,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.02,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.02,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.02,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.05,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.05,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.05,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.05,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.1,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.1,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.1,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.1,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=0.1,n_estimators=50)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=15)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=30)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1.5,n_estimators=50)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=15)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=30)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=500)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_AdaBoost(learning_rate=1,n_estimators=50)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.4)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.6)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.7)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.8)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=0.9)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_BernoulliNB(alpha=1.0)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=3,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=3,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=4,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=4,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=5,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=5,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=6,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=6,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=7,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=7,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=8,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=entropy,max_depth=8,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=3,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=3,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=4,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=4,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=5,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=5,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=6,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=6,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=7,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=7,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=8,min_samples_split=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_DecisionTree(criterion=gini,max_depth=8,min_samples_split=3)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=entropy,n_estimators=128)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=entropy,n_estimators=16)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=entropy,n_estimators=32)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=entropy,n_estimators=64)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=entropy,n_estimators=96)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=128)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=16)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=32)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=64)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_ExtraTrees(criterion=gini,n_estimators=96)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianNB(var_smoothing=0.001)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianNB(var_smoothing=0.01)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianNB(var_smoothing=0.1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianProcess(kernel=sklearn.gaussian_process.kernels.DotProduct,max_iter_predict=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianProcess(kernel=sklearn.gaussian_process.kernels.Matern,max_iter_predict=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianProcess(kernel=sklearn.gaussian_process.kernels.RationalQuadratic,max_iter_predict=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_GaussianProcess(kernel=sklearn.gaussian_process.kernels.RBF,max_iter_predict=1000)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=15,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=30,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=auto,leaf_size=45,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=15,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=30,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=ball_tree,leaf_size=45,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=1,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=1,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=1,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=2,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=2,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=2,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=15,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=1,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=1,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=1,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=2,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=2,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=2,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=30,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=1,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=1,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=1,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=2,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=2,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=2,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=brute,leaf_size=45,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=15,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=30,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=3,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=3,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=3,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=4,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=4,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=4,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=5,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=5,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=5,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=6,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=6,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=6,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=7,p=1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=7,p=2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_KNeighbors(algorithm=kd_tree,leaf_size=45,n_neighbors=7,p=5)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=lbfgs)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=liblinear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=100,max_iter=500,solver=newton-cg)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=lbfgs)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=liblinear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=liblinear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=newton-cg)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_LogisticRegression(C=10,max_iter=500,solver=newton-cg)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=120)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=150)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=60)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=entropy,n_estimators=80)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=100)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=120)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=150)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=200)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=60)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_RandomForest(criterion=gini,n_estimators=80)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.0001,loss=hinge,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.0001,loss=log,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.0001,loss=log,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=hinge,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=hinge,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=hinge,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=log,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=log,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=log,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=modified_huber,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=squared_hinge,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=squared_hinge,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.001,loss=squared_hinge,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=hinge,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=hinge,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=hinge,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=modified_huber,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=modified_huber,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=modified_huber,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=squared_hinge,penalty=elasticnet)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=squared_hinge,penalty=l1)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SGD(alpha=0.01,loss=squared_hinge,penalty=l2)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.2,kernel=linear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.2,kernel=poly)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.2,kernel=rbf)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.2,kernel=sigmoid)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.5,kernel=linear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.5,kernel=poly)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.5,kernel=rbf)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=0.5,kernel=sigmoid)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1.5,kernel=linear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1.5,kernel=poly)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1.5,kernel=rbf)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1.5,kernel=sigmoid)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=linear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=poly)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=rbf)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=1,kernel=sigmoid)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=linear)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=poly)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=rbf)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
 'skgrid_SVC(C=2,kernel=sigmoid)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c'])
    : 'skgrid_LogisticRegression(C=1.0,max_iter=500,solver=lbfgs)|skgrid_THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|THCA.tsv_skgrid_THCA_fbedeBIC_perplatformMUTA|2021-01-13|c',
}



for i_team in range(0, len(team_options)):

    # init conversion dictionary
    perform2imp = {}

    for cancer in cancer_list:
        print(cancer)
        # Select Team
        selected_team = team_options[i_team]
        selected_team_imp = team_imp_options[i_team]
        print('{} and {} selected from list'.format(selected_team, selected_team_imp))

        # Set up out file
        f_out = 'src/conversions/' + selected_team + '.json'

        # Find top model for a team
        f_top = 'data/figure_panel_a/best_models_{}.tsv'.format(cancer)
        df_top = pd.read_csv(f_top, sep='\t', index_col=0)
        top_models = list(df_top.columns)
        print('Best models options:\n')
        for t in top_models:
            print(t)

        # Create source dictionary of conversions
        # k:v == team_options to model prefix in top_models
        mini_conversion_prefix = {
            'aklimate' : 'AKLIMATE',
            'CF' : 'CF', #no change
            'jadbio' : 'jadbio', #no change
            'subSCOPE' : 'subSCOPE',
            'skgrid' :'skgrid' #nochange
        }


        # Select Team to work on
        team_prefix = mini_conversion_prefix[selected_team]
        for m in top_models:
            if m.startswith(team_prefix):
                selected_model = m
                print(selected_model, '\nwas assigned to\n', selected_team)
                exit

        # Subset for team and cancer
        df_imp = df_imp_raw[df_imp_raw['method']==selected_team_imp]
        keep = []
        for v in df_imp['feature_importance_ID']:
            if cancer in v:
                keep.append(v)
        df_imp= df_imp[df_imp['feature_importance_ID'].isin(keep)].reset_index(drop=True)
        print('{} models found'.format(df_imp.shape[0]))

        if selected_team != 'skgrid':
            # Find row that matches with selected_model description
            # src dictionary specfic to team key words in feature_importance_ID column
            substring_dict = {
                'aklimate' : { # no MIR importances reported
                    'CNVR_ONLY' : 'CNVR_ONLY', #nochange
                    'GEXP_ONLY' : 'GEXP_ONLY', #nochange
                    'METH_ONLY' : 'METH_ONLY', #nochange
                    'MULTI_DATA' : 'MULTI_DATA', #nochange
                },
                'CF' :{
                    'All' : 'All', #nochange
                    'CNVR' : 'CNVR', #nochange
                    'GEXP' : 'GEXP', #nochange
                    'METH' : 'METH', #nochange
                    'MIR' : 'MIR', #nochange
                    'MUTA' : 'MUTA', #nochange
                },
                'jadbio' : {
                    'CNVR' : 'CNVR', #nochange
                    'GEXP' : 'GEXP', #nochange
                    'METH' : 'METH', #nochange
                    'MIR' : 'MIR', #nochange
                    'MUTA' : 'MUTA', #nochange
                    'MULTIDATATYPE' : 'MULTIDATATYPE', #nochange
                },
                'subSCOPE' : {
                    'CNVR' : 'CNVR', #nochange
                    'GEXP' : 'GEXP', #nochange
                    'METH' : 'METH', #nochange
                    'MIR' : 'MIR', #nochange
                    'MUTA' : 'MUTA', #nochange
                    'ENSEMBLE' : 'ENSEMBLE', #nochange
                }
            }

            # src dictionary if no matches from substring_dict. these are the assumed values
            gap_substring_dict = {
                'aklimate' : 'MULTI_DATA'
            }


            # 1. Find substring present in selected model
            found = 'false'
            for potential_substring in substring_dict[selected_team].keys():
                if potential_substring in selected_model:
                    lookup_key = potential_substring
                    found = 'true'
                    exit
            if found == 'false': # if no hits from above
                lookup_key = gap_substring_dict[selected_team]
                print('uses this')
            # 2. Use that to find substring to use in df_imp
            df_lookup_key = substring_dict[selected_team][lookup_key]
            print(df_lookup_key)

            # 3. Find matching model and add to perform2imp
            for i in range(0, df_imp.shape[0]):
                if df_lookup_key in df_imp['feature_importance_ID'][i]:
                    df_model = df_imp.iloc[i,:]['feature_importance_ID']
                    print('at row {} found match of\n{}\n\tto\n{}'.format(i, selected_model, df_model))
                    perform2imp[selected_model]= df_imp.iloc[i,:]['feature_importance_ID']
                    exit

            # Output conversion keys - will overwrite old file with new one each loop
            with open(f_out, 'w') as out:
                out.write(json.dumps(perform2imp))
                out.write('\n')


        else: #for skgrid only
            info , selected_skgrid_model = skgrid_get_more_model_info(cancer)
            print(info)

            # If ties then resolve from past manual picking
            if type(selected_skgrid_model)== list:
                if tuple(selected_skgrid_model) in previous_selected:
                    selected_skgrid_model = previous_selected[tuple(selected_skgrid_model)]
                    print('updated to {}'.format(selected_skgrid_model))

            # Input skgrid selected_model object - redundant
            input_model = selected_skgrid_model

            # Remove skgrid prefix
            modified_model = '_'.join(input_model.strip().split('_')[1:])

            # Update ML function name if needed
            func_name = modified_model.split('(')[0]
            new_func_name = skgrid_function_dict[func_name]
            modified_model = new_func_name +'('+'('.join(modified_model.split('(')[1:])

            # Clean up ft selection and other misc info in model string
            modified_model = modified_model.strip().split('|')
            del modified_model[1] # remove second item

            new_ft_method = modified_model[1].split('.tsv_')[1] # remove cancer.tsv_ prefix
            final_model = '|'.join([modified_model[0], new_ft_method, modified_model[2], modified_model[3]])

            # Now look up model in df
            df_imp =df_imp[df_imp['feature_importance_ID']==final_model].reset_index(drop=True)
            assert df_imp.shape[0] == 1

            # 3. Find matching model and add to perform2imp
            perform2imp[selected_model]=df_imp['feature_importance_ID'][0]

            # Output conversion keys
            with open(f_out, 'w') as out:
                out.write(json.dumps(perform2imp))
                out.write('\n')
