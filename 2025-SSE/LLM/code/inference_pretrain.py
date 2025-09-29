from unsloth import FastLanguageModel
from transformers import TextStreamer
from transformers import StoppingCriteriaList

from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig
from peft import PeftModel
import torch

base_id = "unsloth/mistral-7b-instruct-v0.3-bnb-4bit"
model_type = input("Enter model type: ")
adapter_path = f"./pretrain_model/{model_type}/best_model"

if model_type == 'mistral_formula':
  feature_name = 'FORMULA'
elif model_type == 'mistral_for_sym':
  feature_name = 'FORMULA_SYMMETRY'
elif model_type == 'mistral_for_dis':
  feature_name = 'FORMULA_disorder'
elif model_type == 'mistral_for_sym_dis':
  feature_name = 'FORMULA_SYMMETRY_disorder'
    
bnb = BitsAndBytesConfig(load_in_4bit=False, bnb_4bit_compute_dtype=torch.bfloat16)
tokenizer = AutoTokenizer.from_pretrained(base_id)
base = AutoModelForCausalLM.from_pretrained(
    base_id, quantization_config=bnb, device_map="auto", torch_dtype=torch.bfloat16
)
model = PeftModel.from_pretrained(base, adapter_path)

model = FastLanguageModel.for_inference(model)

text_streamer = TextStreamer(tokenizer)

alpaca_prompt = """
<s>[INST]
{}

{}
[/INST]
{}
"""
from datasets import load_dataset
import re
from typing import Union

def formatting_prompts_func(examples):
    instructions = examples['instruction']
    inputs = examples[feature_name]
    outputs = examples['CONDUCTIVITY']
    texts =[]
    for instruction, input_, output in zip(instructions, inputs, outputs):
        text = alpaca_prompt.format(instruction, input_, output) + EOS_TOKEN
        texts.append(text)
    return {"text": texts,}
    
EOS_TOKEN = tokenizer.eos_token
    
dataset = load_dataset("json", data_files={'validation': './data/test.jsonl'})
dataset = dataset.map(
    formatting_prompts_func,
    batched=True,
)

train_conductivity_value_pred_list_total = []
train_conductivity_value_val_list_total =[]

NUMBER = r'[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?'    
    

valid_conductivity_value_pred_list_total = []
valid_conductivity_value_val_list_total =[]

for k in range(1):
  valid_conductivity_value_pred_list=[]
  valid_conductivity_value_val_list=[]

  for i in range(len(dataset['validation']['instruction'])):
    valid_instruction = dataset['validation']['instruction'][i]
    valid_input = dataset['validation'][feature_name][i]
    valid_output = dataset['validation']['CONDUCTIVITY'][i]
  
    prompt = alpaca_prompt.format(valid_instruction, valid_input,"")
    inputs = tokenizer([prompt], return_tensors="pt").to(model.device)
  
    generated_ids = model.generate(**inputs, streamer = text_streamer,max_new_tokens = 12, do_sample=True, temperature=0.7 , top_k=10)
    generated_text = tokenizer.decode(generated_ids[0], skip_special_tokens=True)
  
    match = re.findall(r"[-+]?\d*\.\d+|\d+", generated_text)
    valid_conductivity_value_pred = match[-1]
    valid_conductivity_value_pred_list.append(valid_conductivity_value_pred)
  
    valid_conductivity_value_val_list.append(valid_output)
  
  valid_conductivity_value_pred_list_total.append(valid_conductivity_value_pred_list)
  valid_conductivity_value_val_list_total.append(valid_conductivity_value_val_list)

import pandas as pd

llama3_dataset_test = pd.DataFrame({'valid_0':valid_conductivity_value_pred_list_total[0],'valid_value':valid_conductivity_value_val_list_total[0]})

llama3_dataset_test.to_csv(f'./result/{model_type}_test.csv',index=False)