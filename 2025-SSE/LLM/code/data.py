import preprocess
import pandas as pd
import numpy as np
import re
import json

cif_directory = './cif'
feature_directory = './feature'

train_list= f'{cif_directory}/train_list.txt'
test_list=f'{cif_directory}/test_list.txt'

train_cif_list=preprocess.txt2list(train_list)
test_cif_list=preprocess.txt2list(test_list)

train_disorder_list = preprocess.classify_disorders_from_cif(train_cif_list, cif_directory)
test_disorder_list = preprocess.classify_disorders_from_cif(test_cif_list, cif_directory)

train_formula_list = preprocess.formula_from_cif(train_cif_list, cif_directory)
test_formula_list = preprocess.formula_from_cif(test_cif_list, cif_directory)

train_symmetry_list = preprocess.symmetry_from_cif(train_cif_list, cif_directory)
test_symmetry_list = preprocess.symmetry_from_cif(test_cif_list, cif_directory)

feature = pd.read_csv(f'{feature_directory}/feature.csv')

train_conductivity_list = preprocess.conductivity_from_feature(train_cif_list, feature)
test_conductivity_list = preprocess.conductivity_from_feature(test_cif_list, feature)

train_formula_symmetry_list = preprocess.formula_symmetry(train_formula_list, train_symmetry_list)
test_formula_symmetry_list = preprocess.formula_symmetry(test_formula_list, test_symmetry_list)

train_formula_disorder_list = preprocess.formula_disorder(train_formula_list, train_disorder_list)
test_formula_disorder_list = preprocess.formula_disorder(test_formula_list, test_disorder_list)

train_formula_symmetry_disorder_list = preprocess.formula_symmetry_disorder(train_formula_list, train_symmetry_list, train_disorder_list)
test_formula_symmetry_disorder_list = preprocess.formula_symmetry_disorder(test_formula_list, test_symmetry_list, test_disorder_list)

###dataset
train = preprocess.get_Dataset(train_cif_list, train_formula_list, train_symmetry_list,train_disorder_list, train_formula_symmetry_list,train_formula_disorder_list,train_formula_symmetry_disorder_list, train_conductivity_list)

test = preprocess.get_Dataset(test_cif_list, test_formula_list, test_symmetry_list, test_disorder_list,test_formula_symmetry_list, test_formula_disorder_list, test_formula_symmetry_disorder_list, test_conductivity_list)
print(train)
print(test)
preprocess.get_alpaca_jsonl(train, './data/train.jsonl')
preprocess.get_alpaca_jsonl(test, './data/test.jsonl')
