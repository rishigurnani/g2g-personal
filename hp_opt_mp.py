import os
import sys
from multiprocessing import Pool

args = sys.argv
n_core = int(args[1])

params = {'lr': [.001, .005, .0002],
         'batch_size': [8, 32, 128],
          'depthT': [3, 6, 10],
          'depthG': [3, 5, 8]
         }

epochs = 5

n_decode = 10

bg_path = '~/CS6250_project/models/trial2/all_features/model.pkl'

best_results = []

cwd = os.getcwd() + '/'

round_tup = []

def validate(data):
    lr = data[0]
    batch_size = data[1]
    depthT = data[2]
    depthG = data[3]
    path = cwd + 'lr_%s_bs_%s_depthT_%s_depthG_%s/' %(lr, batch_size, depthT, depthG)
    models_dir = path+'newmodels'
    data_dir = cwd + 'data/'
    processed_dir = cwd + 'processed/'
    results_dir = path+'results/'
    
    #make appropriate directories if they don't already exist
    if not os.path.isdir(path):
        os.makedirs(path)
    if not os.path.isdir(models_dir):
        os.makedirs(models_dir)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    ##########################################################
    
    os.chdir(results_dir)    
    
    #find best model from all epochs having hyperparameters listed above
    os.system('python /home/rgur/CS6250_project/scripts/val_script.py %s %s %s %s %s > %sbest_of_round.txt' %(models_dir, str(epochs), n_decode, bg_path, data_dir, path) )
    ####################################################################
    
    with open('%sbest_of_round.txt' %path, 'r') as f:
                    for line in f:
                        l = line.strip()
                        if 'Epoch' in l:
                            n_epoch = l.split(': ')[1]
                        elif 'Accuracy' in l:    
                            best_acc  = l.split(': ')[1]
    return (path, n_epoch, best_acc)
                
for lr in params['lr']:
    for batch_size in params['batch_size']:
        for depthT in params['depthT']:
            for depthG in params['depthG']:
                
                round_tup.append((lr, batch_size, depthT, depthG))
                
pool = Pool(n_core)
best_results = pool.map(validate, round_tup)
pool.close()
pool.join()                

srt = sorted(best_results, key=lambda x: x[2], reverse=True)
print srt[0]