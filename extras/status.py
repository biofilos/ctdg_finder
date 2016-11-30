from glob import glob
import sys
from termcolor import colored
total = 0

in_dir = sys.argv[1]
out_dir = sys.argv[2]
finished = [x.split('/')[-1] for x in glob("{}/{}/*".format(in_dir, out_dir))]
for run in glob("{}/run*/".format(in_dir)):
    fastas = [x.split('/')[-1].replace(".fa", '') for x in glob("{}*.fa".format(run))]
    run_finished = len(set(finished).intersection(set(fastas)))
    total_run = len(fastas)
    run_name = run.split("/")[-2]
    perc_finished = round(run_finished / total_run,4) * 100
    if perc_finished > 99:
        color = 'green'
    elif perc_finished > 50:
        color = 'yellow'
    else:
        color = 'red'
    if len(fastas) == run_finished:
        msg1 = '{}: Done'.format(run_name)
    else:
        msg1 = "{}: {} / {} ({:}%)".format(run_name, run_finished, total_run, str(perc_finished)[:5])
    print(colored(msg1, color))
    total += total_run

perc_done = round(len(finished) / total, 4) * 100
print("\nTotal done: {:>2} ({}%)".format(len(finished), str(perc_done)[:5]))
print("Remaining: {:>6}".format(total - len(finished)))
print("Total: {:>10}".format(total))
