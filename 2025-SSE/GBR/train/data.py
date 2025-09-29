import os
import csv
import time
import joblib
import random
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

def data_prepare(data_csv,verbos=False, saveto: str="dataset4IC.csv")-> pd.DataFrame:
    df = pd.read_csv(data_csv)
    feature_cols = [c for c in df.columns if c != "log_target"]
    ordered_cols = feature_cols + ["log_target"]
    df_all_features = df.loc[:, ordered_cols].reset_index(drop=True)
    if verbos:
        print(df_all_features)
    if saveto:
        df_all_features.to_csv(saveto, index=True, index_label='Number')
    return df_all_features

def data_split(data_csv="dataset4IC.csv", ratio = 0.8, verbos = True):
    df = pd.read_csv(data_csv)
    N_materials = df.shape[0]
    N_features = df.shape[1] - 3
    diverse_ratio = ratio
    data_file_name = data_csv
    if os.path.exists("split.txt"):
        diverse_set=[]
        remaining_set=[]
        txt = open("split.txt",'r').read()
        s1=txt.find("[",0)
        s2=txt.find("]",s1)
        diverse_set=txt[s1+1:s2].split(", ")
        diverse_set=[int(i) for i in diverse_set]
        s3=txt.find("[",s2)
        s4=txt.find("]",s3)
        remaining_set=txt[s3+1:s4].split(", ")
        remaining_set=[int(i) for i in remaining_set]
    else:
        with open(data_file_name) as f:
            data_file = csv.reader(f)
            next(data_file)
            data = np.empty((N_materials, N_features))
            for i, d in enumerate(data_file):
                data[i] = np.asarray(d[2:2+N_features], dtype=np.float64)
        features_norm = []
        for j in range(N_features):
            col = col = data[:, j]
            col_norm = (col - np.min(col)) / (np.max(col) - np.min(col) + 1e-12)
            features_norm.append(col_norm.reshape(1, N_materials))
        
        x = np.concatenate(features_norm, axis=0)
        N_sample = int(N_materials * diverse_ratio)-1
        
        time.sleep(1)
        diverse_set = []
        remaining_set = list(range(N_materials))
        idx_init = random.sample(list(np.arange(N_materials)),1)[0]
        diverse_set.append(idx_init)
        remaining_set.remove(idx_init)
        N_diverse = 1
        while N_diverse <= N_sample:
            print("Selecting point ", N_diverse)
            min_d_to_diverse_set = np.zeros((N_materials-N_diverse,))
            for i in range(N_materials - N_diverse):
                d_from_each_diverse_pt = np.linalg.norm(x[:,diverse_set] - x[:,remaining_set[i]].reshape(N_features,1),axis=0)
                min_d_to_diverse_set[i] = np.min(d_from_each_diverse_pt)
            idx_select = remaining_set[np.argmax(min_d_to_diverse_set)]
            assert (len(remaining_set) == np.size(min_d_to_diverse_set))
            diverse_set.append(idx_select)
            remaining_set.remove(idx_select)
            N_diverse += 1
        if verbos:
            print("Total number of materials : ", data.shape[0])
            print("Number of features: ", N_features)
        with open("split.txt", "w") as f:
            f.write(str(diverse_set)+" "+str(remaining_set))
            
    all_cols = df.columns.tolist()
    num_col, comp_col = all_cols[0], all_cols[1]
    target_col = all_cols[-1]
    feature_cols = all_cols[2:-1]

    X = df[feature_cols].to_numpy(dtype=np.float64)
    y = df[[target_col]].to_numpy(dtype=np.float64)

    diverse_idx = list(map(int, diverse_set))
    remaining_idx = list(map(int, remaining_set))
    X_train, y_train = X[diverse_idx], y[diverse_idx]
    X_test,  y_test  = X[remaining_idx],  y[remaining_idx]
    df_Xtrain = pd.DataFrame(X_train, columns=feature_cols)
    df_Ytrain = pd.DataFrame(y_train, columns=[target_col])
    df_Xtest  = pd.DataFrame(X_test,  columns=feature_cols)
    df_Ytest  = pd.DataFrame(y_test,  columns=[target_col])
    
    return df_Xtrain, df_Ytrain, df_Xtest, df_Ytest


def normal(Xtrain, Ytrain, Xtest, Ytest, tar="log_target", scaler_path="scaler.gz"):
    x_cols = list(Xtrain.columns) if hasattr(Xtrain, "columns") else None
    idx_tr = Xtrain.index if hasattr(Xtrain, "index") else None
    idx_te = Xtest.index if hasattr(Xtest, "index") else None

    Xtr_np = Xtrain.values if isinstance(Xtrain, pd.DataFrame) else np.asarray(Xtrain, dtype=float)
    Xte_np = Xtest.values  if isinstance(Xtest,  pd.DataFrame) else np.asarray(Xtest,  dtype=float)

    if os.path.exists(scaler_path):
        scaler = joblib.load(scaler_path)
        Xtr_scaled = scaler.transform(Xtr_np)

    else:
        scaler = StandardScaler()
        Xtr_scaled = scaler.fit_transform(Xtr_np)
        joblib.dump(scaler, scaler_path)

    Xte_scaled = scaler.transform(Xte_np)

    if x_cols is not None:
        df_Xtrain = pd.DataFrame(Xtr_scaled, columns=x_cols, index=idx_tr)
        df_Xtest  = pd.DataFrame(Xte_scaled, columns=x_cols, index=idx_te)

    else:
        df_Xtrain = pd.DataFrame(Xtr_scaled)
        df_Xtest  = pd.DataFrame(Xte_scaled)


    if isinstance(Ytrain, pd.DataFrame):
        ytr_series = Ytrain[tar] if tar in Ytrain.columns else pd.Series(Ytrain.squeeze(), name=tar, index=Ytrain.index)

    else:
        ytr_series = pd.Series(np.asarray(Ytrain).squeeze(), name=tar)

    if isinstance(Ytest, pd.DataFrame):
        yte_series = Ytest[tar] if tar in Ytest.columns else pd.Series(Ytest.squeeze(), name=tar, index=Ytest.index)

    else:
        yte_series = pd.Series(n.asarray(Ytest).squeeze(), name=tar)()
    df_Ytrain = ytr_series.to_frame()
    df_Ytest  = yte_series.to_frame()

    return df_Xtrain, df_Ytrain, df_Xtest, df_Ytest

