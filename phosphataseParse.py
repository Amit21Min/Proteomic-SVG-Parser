import re

# Import txt, split the text file
filename = "/Users/amitmin/Documents/research/dataParsing/phosphataseModify.txt"
infile = open(filename, 'r')
lines = infile.readlines()
#kinaseID and branchCoords
lines2 = lines[217:2316]
#nodeX and nodeY values
lines3 = lines[15:215]
#textX and textY values
lines4 = lines[2317:2519]
# Dictionaries
#pathID = {}

kinase = {}

patternPath = re.compile('F_(.*?)"')
for line in lines2:
	if re.match(r'\t', line):
		if re.search(patternPath, line):
			# Make sure to use the join method to convert the returned List from findall into a string
			tempPath = ''.join(re.findall(patternPath, line))
			#pathID[tempPath] = tempPath
			kinase[tempPath] = {'pathID': tempPath}

kinase['PPP3CC_1_'] = {}
kinase['PHLPP1_1_'] = {}
kinase['PTPN11'] = {}
kinase['SCP2'] = {}
kinase['SYNJ2'] = {}
kinase['DUPS21'] = {}

# Parse branchCoords into the dictionary - Sloppy b/c coords span several lines
patternPath2 = re.compile(r' d="(.*?)$')
patternPath3 = re.compile(r'(?<=\t\t)(.*?)$')
coordBuilder = ""
# old and new is necessary to make my parse work for the first kinase ID
oldID = ""
newID = ""
for line in lines2:
	if re.match(r'\t<path', line):
		oldID = newID
		# Make sure to do this join thing to convert the returned List from findall into a string
		newID = ''.join(re.findall(patternPath, line))
		if coordBuilder != "":
			# cut the end characters off and add it to the dictionary
			coordBuilder=coordBuilder[:-3]
			coordBuilder = str(coordBuilder)
			kinase[oldID]['coords'] = coordBuilder
			coordBuilder = ""
			
	# If the line starts with d=, store the rest of the line
	if re.search(patternPath, line):
		coordBuilder += ''.join(re.findall(patternPath2, line))
	# If the line starts with a character? (^< I forgot what this does) then store the rest of the line
	if re.search('[^<]', line):
		coordBuilder += ''.join(re.findall(patternPath3, line))
# This stores the last kinase - no \t<path is matched for the final kinase
coordBuilder=coordBuilder[:-3]
coordBuilder = str(coordBuilder)
kinase[newID]['coords'] = coordBuilder


# Parse node.x and node.y values into dictionary
patternPathKinaseID = re.compile(r'F_(.*?)"')
patternPathNodeX = re.compile(r'cx="(.*?)"')
patternPathNodeY = re.compile(r'cy="(.*?)"')
for line in lines3:
	if re.match(r'\t', line):
		if re.search(patternPathNodeX, line):
			#Find KinaseID, find nodeX and store in dictionary
			kinaseID = ''.join(re.findall(patternPathKinaseID, line))
			tempNodeX = ''.join(re.findall(patternPathNodeX, line))
			#print(kinaseID)
			#print(tempNodeX)
			kinase[kinaseID]['nodeX'] = tempNodeX

			#kinase[kinaseID]['nodeX'] = tempNodeX
		if re.search(patternPathNodeY, line):
			kinaseID = ''.join(re.findall(patternPathKinaseID, line))
			tempNodeY = ''.join(re.findall(patternPathNodeY, line))
			kinase[kinaseID]['nodeY'] = tempNodeY

# Parse text.x and text.y values into dictionary
patternPathKinaseID = re.compile(r'F_(.*?)"')
patternPathTextX = re.compile(r'1\s0\s0\s1\s(.*?)\s')
#patternTextY = re.compile(r'1 0 0 1 (.*?) "')
for line in lines4:
	if re.match(r'\t', line):
		if re.search(patternPathTextX, line):
			kinaseID = ''.join(re.findall(patternPathKinaseID, line))
			tempTextX = ''.join(re.findall(patternPathTextX, line))
			kinase[kinaseID]['textX'] = tempTextX
			#The line below took me over an hour to make :/ I had to put in a \s and trim off the extra, I can't
			#figure out how to tell it to stop at the parantehsis
			patternPathTextY = tempTextX + r'\s' + r'(.*?)\s'
			if re.search(patternPathTextY, line):
				tempTextY = ''.join(re.findall(patternPathTextY, line))
				tempTextY = tempTextY[:-2]
				kinase[kinaseID]['textY'] = tempTextY

#print(kinase)


tsv = ""
oneList = ["id.coral","id.uniprot","id.ensembl","id.entrez","id.HGNC","id.longname","kinase.group","kinase.family","kinase.subfamily","branch.coords","branch.val","branch.group","branch.col","node.x","node.y","node.group","node.val.col","node.val.radius","node.radius","node.col","node.strokewidth","node.strokecol","text.x","text.y","text.font","text.size","text.label","ode.opacity","node.selected","nodeorder","branchorder"]
tsv += '\t'.join(oneList)
tsv += '\n'


for key in kinase:
	var1 = kinase.get(key)
	varPath = var1.get('pathID')
	varCoord = var1.get('coords')
	varNodeX = var1.get('nodeX')
	varNodeY = var1.get('nodeY')
	varTextX = var1.get('textX')
	varTextY = var1.get('textY')
	if varPath is None:
		print(key + " doesn't have a proper id.Coral value")
	if varCoord is None:
		print(key + " doesn't have a proper coordinate value")
	if varNodeX is None:
		print(key + " doesn't have a proper nodeX value")
	if varNodeY is None:
		print(key + " doesn't have a proper nodeY value")
	if varTextX is None:
		print(key + " doesn't have a proper textX value")
	if varTextY is None:
		print(key + " doesn't have a proper textY value")

print(kinase)
for key in kinase:
	if key == 'DUSP21' or key == 'PTPN1_1_' or key == 'DUPS21' or key == 'SYNJ2' or key == 'SCP2' or key == 'PTPN11' or key == 'PHLPP1_1_' or key == 'PPP3CC_1_' or key == 'ACP7' or key == 'MINPP1' or key == 'SPC2':
		continue
	print(key)
	oneList = kinase[key]['pathID'],'.','.','.','.','.','.','.','.',kinase[key]['coords'],'.','.','.',kinase[key]['nodeX'],kinase[key]['nodeY'] + \
				 '.','.','.','.','.','.','.','.',kinase[key]['textX'],kinase[key]['textY'],'.','.','.','.','.','.','.'
	tsv += '\t'.join(oneList)
	tsv += '\n'
	

text_file = open("PhosphataseTree.tsv", "w")
text_file.write(tsv)