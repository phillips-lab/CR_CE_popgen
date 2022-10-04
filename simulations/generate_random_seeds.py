#' Usage:
#' python generate_random_seeds.py
#' it will output 12000 random seeds

import rstr, itertools, operator, random

#get uniq seeds
def sort_uniq(sequence):
        return map(
             operator.itemgetter(0),
             itertools.groupby(sorted(sequence)))

rlist=[]
n=0

#generate 1500 random numbers with 12 digits
while n<16000:
    rlist.append(rstr.digits(12))
    n=n+1


rlist=[ x for x in rlist if not x.startswith('0')] #remove seeds which begin with '0'
rlist=random.sample(list(sort_uniq(rlist)),12001) #get 1000 unique seeds
rlist.sort()

with open("seeds.txt","w") as txt:
    txt.write('\n'.join(rlist))
