from unsloth import FastLanguageModel
import torch

max_seq_length = 4096  # Choose any! We auto support RoPE Scaling internally!
load_in_4bit = True # Use 4bit quantization to reduce memory usage. Can be False.

feature = input("Enter feature name: ")
model, tokenizer = FastLanguageModel.from_pretrained(
    model_name="mistralai/Mistral-7B-Instruct-v0.3",
    max_seq_length=max_seq_length,
    dtype=None,
    load_in_4bit=load_in_4bit,
    #bnb_4bit_compute_dtype=torch.bfloat16,
    #token = " ", # use one if using gated models like meta-llama/Llama-2-7b-hf
)

model = FastLanguageModel.get_peft_model(
    model,
    lora_alpha=64,
    lora_dropout=0.05,
    target_modules=['q_proj','k_proj','v_proj','o_proj','gate_proj','up_proj','down_proj'],
    bias='none',
    use_gradient_checkpointing='unsloth',
    random_state=123,
    use_rslora=False,
    loftq_config=None)
    
from datasets import load_dataset

EOS_TOKEN = tokenizer.eos_token

alpaca_prompt = """
<s>[INST]
{}

{}
[/INST]

{}</s>
"""

def formatting_prompts_func(examples):
    instructions = examples['instruction']
    inputs = examples[feature]
    outputs = examples['CONDUCTIVITY']
    texts =[]
    for instruction, input_, output in zip(instructions, inputs, outputs):
        text = alpaca_prompt.format(instruction, input_, output) + EOS_TOKEN
        texts.append(text)
    return {"text": texts,}

dataset = load_dataset("json", data_files={"train": "./data/train.jsonl"})
dataset = dataset.map(
    formatting_prompts_func,
    batched=True,
)

from trl import SFTTrainer
from transformers import TrainingArguments

tokenizer.padding_side = 'right'

trainer=SFTTrainer(
    model=model,
    tokenizer=tokenizer,
    train_dataset=dataset['train'],
    dataset_text_field='text',
    max_seq_length=max_seq_length,
    dataset_num_proc=2,
    packing=False,
    args=TrainingArguments(
        per_device_train_batch_size=4,
        gradient_accumulation_steps=2,
        warmup_steps=5,
        num_train_epochs=30,
        logging_strategy ='epoch',
        learning_rate=2e-4,
        fp16= not  torch.cuda.is_bf16_supported(),
        bf16 = torch.cuda.is_bf16_supported(),
        optim='adamw_8bit',
        weight_decay=0.01,
        lr_scheduler_type='cosine',
        seed=123,
        output_dir='outputs',)
)

trainer_stats = trainer.train() 
print(f'train metrics : {trainer_stats}')

# upload to local
model.save_pretrained('./trained_model/best_model')
