import shutil
import sys
import os
from glob import glob

seq_dir = sys.argv[1]
all_seqs = glob("{}/*.fa".format(seq_dir))
groups = sys.argv[2]

per_group = (len(all_seqs) / int(groups)) + 1
print("{} sequences".format(len(all_seqs)))
count = 1

for seq in all_seqs:
    seq_name = seq.split('/')[-1]
    dir_name = seq.replace(seq_name, '')
    if count < per_group:
        folder = "1"
    else:
        for i in range(1,int(groups)+1):
            if count < per_group * i:
                folder = str(i)
                break
    new_dir = "{}run{}/{}".format(dir_name, folder, seq_name)
    print(new_dir)
    if not os.path.exists("{}run{}".format(dir_name, folder)):
        os.makedirs("{}run{}".format(dir_name, folder))
    shutil.move(seq, new_dir)
    count += 1

