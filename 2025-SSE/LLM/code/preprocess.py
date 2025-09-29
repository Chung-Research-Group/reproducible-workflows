import pandas as pd
import numpy as np
import re
import json
# -*- coding: utf-8 -*-
def txt2list(txt):
  cif_list=[]
  with open(txt, 'r') as f:
    for line in f:
      cif_list.append(line.strip())
  return cif_list
  
# Only English comments below

import os
import re
from collections import defaultdict
  
   
def classify_disorders_from_cif(cif_file_list, cif_directory):
  atom_list = [
    'H','Li','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ca',
    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','Se','Br',
    'Sr','Y','Zr','Nb','Mo','Ag','In','Sn','Sb','Te','I','Ba','Hf','Ta',
    'W','Pb','Bi','La','Ce','Pr','Nd','Sm','Gd','Er','Lu']
  disorder_list = []

  for name in cif_file_list:
      path = os.path.join(cif_directory, name)
      try:
          with open(path, "r") as f:
              lines = [ln.strip() for ln in f]
      except FileNotFoundError:
          disorder_list.append("Error: CIF file not found.")
          continue

      # initialize ONCE per file (not per line)
      has_li_disorder = False
      has_other_disorder = False
      in_disorder_section = False

      for line in lines:
          if not line or line.startswith("#"):
              # blank or comment
              if in_disorder_section:
                  in_disorder_section = False
              continue

          # start of disorder section (very naive: exactly one header line)
          if line.startswith("_atom_site_occupancy"):
              in_disorder_section = True
              continue

          # if another header/loop/comment appears, end the section
          if in_disorder_section and (line.startswith("_") or line.lower().startswith("loop_")):
              in_disorder_section = False
              # fall through to process this new header/loop in next iterations
              continue

          # parse data lines only while inside the section
          if in_disorder_section:
              parts = line.split()
              if len(parts) >= 2:
                  atom_symbol = parts[0]
                  occ_str = parts[-1]
                  try:
                      occ = float(occ_str)
                      if occ != 1.0:
                          if atom_symbol == "Li":
                              has_li_disorder = True
                          else:
                              has_other_disorder = True
                  except ValueError:
                      # not a numeric disorder, skip
                      pass

      # decide label per file
      if has_li_disorder and has_other_disorder:
          disorder_list.append("Li cation and Other atom disorder")
          #disorder_list.append("0")
      elif has_li_disorder:
          disorder_list.append("Li cation disorder")
          #disorder_list.append("1")
      elif has_other_disorder:
          disorder_list.append("Other atom disorder")
          #disorder_list.append("2")
      else:
          disorder_list.append("No disorder")
          #disorder_list.append("3")
  return disorder_list

def formula_from_cif(cif_file_list, cif_directory):
  formula_list=[]
  for name in cif_file_list:
    name = f'{cif_directory}/{name}'
    with open(name, 'r') as f:
      lines = f.readlines()
      for line in lines:
        if line.startswith("_chemical_formula_structural"):
          formula = line.strip().split()[-1]
          formula_list.append(formula)
  return formula_list  
                
def symmetry_from_cif(cif_file_list, cif_directory):
  symmetry_list=[]
  for name in cif_file_list:
    name = f'{cif_directory}/{name}'
    with open(name, 'r') as f:  
      lines = f.readlines()
      for line in lines:
        if line.startswith("_symmetry_space_group_name_H-M"):
          formula = line.strip().split(maxsplit=1)
          formula = formula[1].strip("'\"")
          symmetry_list.append(formula)
        elif line.startswith("_space_group_name_H-M_alt"):
          formula = line.strip().split(maxsplit=1)
          formula = formula[1].strip("'\"")
          symmetry_list.append(formula)
  return symmetry_list


def conductivity_from_feature(cif_file_list, feature_dataset):
  conductivity_list=[]
  for name in cif_file_list:
    cif = name.replace('.cif', '')
    conductivity = feature_dataset.loc[feature_dataset['file_name'] == cif, 'log_Conductivity']
    if not conductivity.empty:
      val = float(conductivity.values[0])
      conductivity_list.append(f"{val:.2f}")
    else :
      print('condcutivity is nan')
  return conductivity_list 

def formula_symmetry(formula_list, symmetry_list):
  formula_symmetry_list=[]
  for i in range(len(formula_list)):
    formula = formula_list[i]
    symmetry = symmetry_list[i]
    incor = str(f'<formula: {formula}> <symmetry: {symmetry}>')
    formula_symmetry_list.append(incor)   
  return formula_symmetry_list  

def formula_disorder(formula_list, disorder_list):
  formula_disorder_list=[]
  for i in range(len(formula_list)):
    formula = formula_list[i]
    disorder = disorder_list[i]
    incor = str(f'<formula: {formula}> <disordered atom: {disorder}>')
    formula_disorder_list.append(incor)   
  return formula_disorder_list

def formula_symmetry_disorder(formula_list, symmetry_list, disorder_list):
  formula_symmetry_disorder_list=[]
  for i in range(len(formula_list)):
    formula = formula_list[i]
    symmetry = symmetry_list[i]
    disorder = disorder_list[i]
    incor = str(f'<formula: {formula}> <symmetry: {symmetry}> <disordered atom: {disorder}>')
    formula_symmetry_disorder_list.append(incor)   
  return formula_symmetry_disorder_list

def get_Dataset(cif_file_name_list, formula_list, symmetry_list, disorder_list, formula_symmetry_list, formula_disorder_list, formula_symmetry_disorder_list, conductivity_list):
  datset =pd.DataFrame({'cif':cif_file_name_list,
  'formula':formula_list,
  'symmetry':symmetry_list,
  'disorder': disorder_list,
  'formula_symmetry': formula_symmetry_list,
  'formula_disorder': formula_disorder_list,
  'formula_symmetry_disorder': formula_symmetry_disorder_list,
  'conductivity':conductivity_list})
  return datset 


def get_alpaca_jsonl(datset, filename):
  with open(filename, 'w', encoding='utf-8') as f:
    for _, row in datset.iterrows():
      cif = row['cif']
      formula = row['formula']
      symmetry = row['symmetry']
      disorder = row['disorder']
      formula_symmetry = row['formula_symmetry']
      formula_disorder = row['formula_disorder']
      formula_symmetry_disorder = row['formula_symmetry_disorder']
      conductivity = row['conductivity']
            
      data = {
      "instruction": (
      "You are a domain expert in materials science specializing in ionic transport. "
      "Your task is to predict the logarithm (base-10) of the ionic conductivity of a solid electrolyte. "
      "You will be given the material's chemical formula and property description. "
      "Output ONLY a single numeric value (the predicted log10(conductivity)) on one line. "
      "Do not print units, JSON, labels, extra spaces, or any additional text. "
      "Be as precise as possible (use sufficient decimal places)"), 
      "CIF": f"{cif}",
      "FORMULA": f"{formula}",
      "SYMMETRY": f"{symmetry}",
      "DISORDER": f"{disorder}",
      "FORMULA_SYMMETRY": f"{formula_symmetry}",
      "FORMULA_disorder": f"{formula_disorder}",
      "FORMULA_SYMMETRY_disorder": f"{formula_symmetry_disorder}",      
      "CONDUCTIVITY": f"{conductivity}",
      }
      f.write(json.dumps(data, ensure_ascii=False) + '\n')





