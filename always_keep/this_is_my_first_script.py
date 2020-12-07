# -*- coding:utf-8 -*_
'''
'''
# convert csv file to dict
# build a dict of list like {key:[...element of lst_inner_value...]}
import csv
def row_csv2dict(csv_file):
    dict_club={}
    with open(csv_file) as f:
        reader=csv.reader(f,delimiter=',')
        for row in reader:
            dict_club[row[0]]=row[1]
    return dict_club

curaw_dict = row_csv2dict('curaw.csv')
#print curaw_dict

input_file = open('curaw1.txt',"r").read()

a =input_file.split('\n')

def replace(lst, dictionary):
    for k,v in enumerate(lst):
        if v in dictionary:
         lst[k] = dictionary[v]
    return lst

#for i in replace(a,curaw_dict):
#    print i

with open('cucumber_raw.txt', 'r') as f1, open('cucumber_results.txt', 'w') as f2:
    for line in f1.readlines():
        line = line.strip()
#print line
        if 'mi' in line:
            f2.write(line + '\n')








