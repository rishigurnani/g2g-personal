import sys
import os

args = sys.argv

DIR=args[1] #model directory
NUM=int(args[2]) #number of models to test
N_DECODE=args[3] #number of monomers per test monomers #number of test monomers
BG_PATH=args[4] #path to bandgap predictor
DATA_DIR=args[5]

with open('%stest.txt' %DATA_DIR, 'r') as f:
    lines = f.readlines()
    total_n = len([l for l in lines if l.strip(' \n') != ''])

max_acc = 0.0
best_epoch = 0

for i in range(NUM):
    #i += 7 #comment out later
    f="%s/model.iter-%s" %(DIR, str(i))
    print(f)
    if os.path.isfile(f):
        os.system('python ~/iclr19-graph2graph/diff_vae/decode.py --num_decode %s --test %stest.txt --vocab %svocab.txt --model %s --use_molatt | python ~/CS6250_project/scripts/bg_score.py %s > results.%s' %(N_DECODE, DATA_DIR, DATA_DIR, f, BG_PATH, str(i)))
        os.system('python ~/CS6250_project/scripts/bg_analyze.py --num_decode %s --sim_delta .2 --prop_delta 6 --total_n %s < results.%s > analyze.%s' %(N_DECODE, total_n, str(i), str(i)) )
        with open('analyze.%s' %(str(i)) ) as f:
            for line in f:
                pass
                last_line = line   
        
        print last_line
        acc = float(last_line.strip().split('success rate ')[1])
        if acc > max_acc:
            max_acc = acc
            best_epoch = i

print "Epoch with best model: ", best_epoch
print "Accuracy for best model: ", max_acc
    