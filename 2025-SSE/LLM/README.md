## Example files for LLM Fine-tuning
1. ./cif/train_list.txt : add .cif list for train
2. ./cif/test_list.txt : add .cif list for test
3. Get features : python3 code/data.py
4. Pretrain : python3 code/inference_pretrain.py 
(modify model_type & input feature)
5. ./result/ : results