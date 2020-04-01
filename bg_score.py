import pandas as pd
import os
import sys
from props import *

import time

args = sys.argv

model_path = args[1] #path to bandgap model

def polymerize(s):
    '''
    Turn monomer SMILES into polymer SMILES
    '''
    return '[*]' + s + '[*]'

def save_smiles_df(l):
    '''
    Save df from list of SMILES
    '''
    n_smiles = len(l)
    ID = list(range(n_smiles))
    pd.DataFrame({"ID": ID, "SMILES": l}).to_csv('./smiles_df.csv')

def create_fp_input():
#     print('file_dataset = ./smiles_df.csv')
#     print('col_smiles = SMILES')
#     print('col_X = aT bT m e')
#     print('col_id = ID')
#     print('file_fingerprint = fp_df.csv')
#     print('polymer_fp_version = 2')
#     print('ismolecule = 0')
#     print('drop_failed_rows = 0')
#     print('ncore = 18')
    f = ('file_dataset = ./smiles_df.csv\n'
    'col_smiles = SMILES\n'
    'col_X = aT bT m e\n'
    'col_id = ID\n'
    'file_fingerprint = fp_df.csv\n'
    'polymer_fp_version = 2\n'
    'ismolecule = 0\n'
    'drop_failed_rows = 0\n'
    'ncore = 18\n')
    text_file = open('fp_input', "w")
    text_file.write(f)
    text_file.close()

def run_fp():
    os.system('fp fp_input')

def fix_fp():
    df = pd.read_csv('fp_df.csv')
    df = df.iloc[df.dropna().index].reset_index().drop('index', axis=1)
    df.to_csv('fp_df_fixed.csv')
    
def create_pred_input(model_path):
    f = ('file_model = ' + model_path + '\n'
    'file_fingerprint  = fp_df_fixed.csv\n'
    'file_output = output.csv\n'
    )
    text_file = open('pred_input', "w")
    text_file.write(f)
    text_file.close()    

def run_pred():
    try:
        os.system('predict pred_input > log')
    except:
        pass

def get_pred():
    
    df = pd.read_csv('output.csv')
    normed = df['y'].tolist()
    return [(val*1.67668) + 4.34115 for val in normed]

def get_all_preds(l):
    poly_l = [polymerize(i) for i in l]
    save_smiles_df(poly_l)

    create_fp_input()

    run_fp()

    fix_fp()
    
    create_pred_input(model_path)

    run_pred()

    return get_pred()    

#perform tasks
ys = []
xs = []
sim2Ds = []
for line in sys.stdin:
    x,y = line.split()
    ys.append(y)
    xs.append(x)
    if y == "None": y = None
    sim2Ds.append(similarity(x, y))

outs = get_all_preds(ys)

for ind, out in enumerate(outs):
    x = xs[ind]
    sim2D = sim2Ds[ind]
    y = ys[ind]
    try:
        print x, y, sim2D, out
    except:
        print x, y, sim2D, 0.0