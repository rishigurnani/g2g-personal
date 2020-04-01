import sys
import argparse
import io

parser = argparse.ArgumentParser()
parser.add_argument('--num_decode', type=int, default=20)
parser.add_argument('--sim_delta', type=float, default=0.4)
parser.add_argument('--prop_delta', type=float, default=0.9)
parser.add_argument('--total_n', type=int, default=0)
args = parser.parse_args()

stdin = sys.stdin
num_decode = args.num_decode
sim_delta = args.sim_delta
prop_delta = args.prop_delta

data = []
start_append = False

len_d = args.total_n
for line in stdin:
    if 'Done' in line:
        start_append = True
    elif start_append:
            data.append(line.split())
#data = [line.split() for line in sys.stdin]
data = [(a,b,float(c),float(d)) for a,b,c,d in data]

n_mols = len_d

n_succ = 0.0
for i in xrange(0, len(data), num_decode):
    set_x = set([x[0] for x in data[i:i+num_decode]])
    assert len(set_x) == 1

    good = [(sim,prop,new) for _,new,sim,prop in data[i:i+num_decode] if 1 > sim >= sim_delta and prop >= prop_delta]
    if len(good) > 0:
        n_succ += 1
    for i in good:
        print i[2]
print 'Evaluated on %d samples' % (n_mols,)
print 'success rate', n_succ / n_mols