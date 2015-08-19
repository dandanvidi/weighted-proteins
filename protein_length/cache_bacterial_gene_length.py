import json

upac2length = {}
for row in open('cache/uniprot_sprot_bacteria.dat', 'r'):
    if row.startswith('ID'):
        length = filter(None, row.strip().split(' '))[-2]
    elif row.startswith('AC'):
        upac = filter(None, row.strip().replace(';', '').split(' '))[-1]
    else:
        continue
    upac2length[upac] = length

with open('cahce/upac2length_bacteria.json', 'wb') as fp:
    json.dump(upac2length, fp)