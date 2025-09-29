from data import data_prepare,data_split,normal
import warnings
import joblib
import numpy as np
import pandas as pd
import re
from scipy.stats import spearmanr
from skopt import gp_minimize
from skopt.utils import use_named_args
from skopt.space import Real, Integer
from xgboost import XGBRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, mean_absolute_error,r2_score
from sklearn.exceptions import ConvergenceWarning

np.int = int
warnings.filterwarnings("ignore", category=FutureWarning, module='xgboost')
warnings.filterwarnings('ignore', category=ConvergenceWarning)

def data(data_csv="real_final_dataset.csv", tar='log_target', verbos=False, saveto="dataset4IC.csv"):
    df_all_features = data_prepare(data_csv,verbos=verbos, saveto=saveto)
    print(df_all_features.head())
    X_train, y_train, X_test, y_test = data_split(data_csv=saveto,ratio = 0.8,verbos = True)
    Xtrain, Ytrain, Xtest, Ytest = normal(X_train, y_train, X_test, y_test,tar=tar)
    return Xtrain, Ytrain, Xtest, Ytest

def _spearman_no_nan(y_true, y_pred):
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    if mask.sum() < 2:
        return np.nan, np.nan
    rho, p = spearmanr(y_true[mask], y_pred[mask])
    return float(rho), float(p)

def gbr(Xtrain, Ytrain, Xtest, Ytest, tar, n_job, call, save = True):
    warnings.filterwarnings("ignore", category=FutureWarning, module='xgboost')
    reg = XGBRegressor()
    space  = [Integer(1,500, name='n_estimators'),
            Integer(1, 10, name='max_depth'),
            Integer(1, 10, name='num_parallel_tree'),
            Integer(1, 10, name='min_child_weight'),
            Real(0.001,1,"log-uniform",name='learning_rate'),
            Real(0.01,1,name='subsample'),
            Real(0.001,10,"log-uniform",name='gamma'),
            Real(0, 1, name='alpha'),
            Real(2, 10, name='reg_alpha'),
            Real(10, 50, name='reg_lambda')
         ]
    @use_named_args(space)
    def objective(**params):
        reg.set_params(**params)
        result=-np.mean(cross_val_score(reg, Xtrain, Ytrain, cv=6, n_jobs=n_job,
                                            scoring="neg_mean_squared_error"))
        print(result)
        return result
    res_gp = gp_minimize(objective, space, n_calls=call)
    print("Best score=%.4f" % res_gp.fun)
    reg_opt = XGBRegressor(n_estimators=res_gp.x[0],
                            max_depth=res_gp.x[1],
                            num_parallel_tree=res_gp.x[2],
                            min_child_weight=res_gp.x[3],
                            learning_rate=res_gp.x[4],
                            subsample=res_gp.x[5],
                            gamma=res_gp.x[6],
                            alpha=res_gp.x[7],
                            reg_alpha=res_gp.x[8],
                            reg_lambda=res_gp.x[9]
                            )
    reg_opt.fit(Xtrain, Ytrain)
    r2_train = r2_score(Ytrain, reg_opt.predict(Xtrain))
    r2_test = r2_score(Ytest, reg_opt.predict(Xtest))
    mae_train = mean_absolute_error(Ytrain, reg_opt.predict(Xtrain))
    mae_test = mean_absolute_error(Ytest, reg_opt.predict(Xtest))
    rmse_train = np.sqrt(mean_squared_error(Ytrain,reg_opt.predict(Xtrain)))
    rmse_test = np.sqrt(mean_squared_error(Ytest, reg_opt.predict(Xtest)))
    spearman_train, spearman_p_train = _spearman_no_nan(Ytrain, reg_opt.predict(Xtrain))
    spearman_test,  spearman_p_test  = _spearman_no_nan(Ytest,  reg_opt.predict(Xtest))

    print("Training: r2,mae,rmse: ", r2_train,mae_train,rmse_train, spearman_train, spearman_p_train)
    print("Test: r2,mae,rmse: ", r2_test,mae_test,rmse_test, spearman_test,  spearman_p_test)

    if save:
        joblib.dump(reg_opt,"gbr_" + tar + ".pkl")
        result = pd.ExcelWriter("gbr_"+tar+".xlsx")
        df_result_train = pd.DataFrame({"Ytrain": Ytrain.values.reshape(-1),
                                        'Ytrain_pre': reg_opt.predict(Xtrain).ravel()})
        df_result_train.to_excel(result, index=False, sheet_name = "data_train")
        df_result_test = pd.DataFrame({"Ytest": Ytest.values.reshape(-1),
                                        'Ytest_pre': reg_opt.predict(Xtest).ravel()})
        df_result_test.to_excel(result, index=False, sheet_name = "data_test")
        feature_importance = reg_opt.feature_importances_
        feature_names = Xtrain.columns
        print("impact: ", feature_names, feature_importance)
        feature_importance_df = pd.DataFrame({"Feature": feature_names, "score": feature_importance})
        feature_importance_df.to_excel(result, index=False, sheet_name = "importance")
        score = pd.DataFrame({"r2_train": [r2_train],
                            "r2_test": [r2_test],
                            'mae_train': [mae_train],
                            'mae_test': [mae_test],
                            'rmse_train': [rmse_train],
                            'rmse_test': [rmse_test],
                            'spearman_train': [spearman_train],
                            'spearman_test':  [spearman_test],
                            'spearman_p_train': [spearman_p_train],
                            'spearman_p_test':  [spearman_p_test],
                            })
        score.to_excel(result, index=False, sheet_name = "score")
        result.close()
    