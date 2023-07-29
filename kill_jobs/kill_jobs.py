import os

job_list = []
with open('tmp.txt') as f:
    for line in f:
        sl = line.split()
        job_list.append(sl[0])

for job in job_list:
    os.system('alien.py kill ' + job)
