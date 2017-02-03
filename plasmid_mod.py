import json

with open('enz_list.json') as data_file:
	enz_list = json.load(data_file)

# TARGET CODE CHANGE

init_site = []
target_site = []

for i in enz_list:
	if 'Arthrospira platensis' in i[0]:
		if '?' not in i[2]:
        		init_site.append(i[2])

for i in init_site:
	i = i.replace('R', '[GA]')
	i = i.replace('Y', '[CT]')
	i = i.replace('M', '[AC]')
	i = i.replace('K', '[GT]')
	i = i.replace('S', '[GC]')
	i = i.replace('W', '[AT]')
	i = i.replace('B', '[GCT]')
	i = i.replace('D', '[AGT]')
	i = i.replace('H', '[ACT]')
	i = i.replace('V', '[ACG]')
	i = i.replace('N', '[ACGT]')
	print i
	# for x in i:
	# 	if 'A' in x:
	# 		x == '2'
	# 		print i
# print init_site
