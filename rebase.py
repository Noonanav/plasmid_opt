import re
import regex
import json

f = open('rebase.txt')
file = (f.read())

enz_list = []

name = regex.finditer(r'<1>.*', file)

for i in name:
    org = regex.search(r'<3>.*', file[i.end():])
    site = regex.search(r'<5>[A-Z\?\^]*', file[i.end():])
    typ = regex.search(r'<8>.*', file[i.end():])
    enzyme = [org.group(0)[3:], i.group(0)[3:], site.group(0)[3:], typ.group(0)[3:]]
    enz_list.append(enzyme)

open('./enz_list.json', 'w').write(json.dumps(enz_list))
