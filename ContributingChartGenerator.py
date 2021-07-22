font = 'sans'
font_size = 7       #px
font_colour = '#000000'
box_size = 30       #px
element = 'fm105'
box_stroke = '#000000'
stroke_width = 0.5
box_fill = '#FFFFFF'
box_space = 2
offset_x = 30
offset_y = 30

import pandas as pd
import naming
import numpy as np
import matplotlib.pyplot as plt
import math
import os

path = '.'
outdir=path+'/A-Z_plot/'
if (os.access(outdir,os.F_OK))==0:
    os.mkdir(outdir)

pd.set_option('display.max_rows', 650)

def _channel_to_hex(color_val: int):
	raw: str = hex(color_val)[2:]
	return raw.zfill(2)

def round_decimals_down(number:float, decimals:int=2):
	"""
	Returns a value rounded down to a specific number of decimal places.
	"""
	if not isinstance(decimals, int):
		raise TypeError("decimal places must be an integer")
	elif decimals < 0:
		raise ValueError("decimal places has to be 0 or more")
	elif decimals == 0:
		return math.floor(number)

	factor = 10 ** decimals
	return math.floor(number * factor) / factor

def rgb_to_hex(red: int, green: int, blue: int):
	return "#" + _channel_to_hex(red) + _channel_to_hex(green) + _channel_to_hex(blue)

def gen_svg(contributing_database, running_percent):
	
	colour_weightings = contributing_database.groupby("Father")["Max CDF contrib"].sum()
	colour_titles = list(colour_weightings.index)
	colour_values = colour_weightings.values
	weights = pd.DataFrame(colour_values, columns = ['Weights'])
	weights.index = colour_titles
	weights['Normalised'] = -1*np.log(weights['Weights'])
	weights = weights.replace([np.inf, -np.inf], np.nan).dropna(subset=["Normalised"], how="all")
	weights['Normalised'] = weights['Normalised']/max(weights['Normalised'])
	bins = list(np.linspace(min(colour_values), max(colour_values), 11))
	bins.reverse()
	colours = [(int(round(0+bins[i]*200)),int(round(255-bins[i]*250)),0) for i in range(len(bins)-1)]
	bins.reverse()
	lines = ['<?xml version="1.0" encoding="utf-8"?>\n','<!DOCTYPE svg\n',"  PUBLIC '-//W3C//DTD SVG 1.1//EN'\n","  'http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd'>\n",'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="3000" height="1920">\n','  <g id="layer0" fill="none">\n','  </g>\n','</svg>\n']
	# print(lines)
	fiss_data = pd.read_csv('./databases/fyu235thermal.txt', sep="	", header=0)
	fiss_data = fiss_data.sort_values(by='Z',ignore_index=True)
	N = fiss_data['A']-fiss_data['Z']
	min_N = min(N)  #X axis
	max_Z = max(fiss_data['Z']) #Y axis
	min_Z = min(fiss_data['Z'])
	max_N = max(N)
	# print(fiss_data)
	fatherCounter = 0
	boxes = []
	title = []
	for i,row in fiss_data.head(n=fiss_data.shape[0]).iterrows():
		try:
			shortName = naming.short_readable(int(row['Z']*10000000+row['A']*10000))
			name = naming.readable(int(row['Z']*10000000+row['A']*10000))
			if name in colour_titles:
				try:
					value = weights.loc[name, 'Normalised']
					print(bins)
					colour_index = [bins.index(bins[i]) for i in range(len(bins)-1) if value>bins[i] and value<=bins[i+1]][0]
					rgb_colour = colours[colour_index]
					box_fill = rgb_to_hex(rgb_colour[0],rgb_colour[1],rgb_colour[2])
					font_colour = '#000000'
					fatherCounter += 1
				except KeyError:
					box_fill = '#FFFFFF'
					font_colour = '#000000'
			else:
				box_fill = '#FFFFFF'
				font_colour = '#000000'
			boxes.append('    <rect id="'+str(shortName)+'" width="'+str(box_size)+'" height="'+str(box_size)+'" stroke="'+box_stroke+'" stroke-width="'+str(stroke_width)+'" fill="'+box_fill+'" x="'+str(int(offset_x+(N[i]-min_N)*(box_size+box_space)))+'" y="'+str(int(((max_Z-min_Z)*(box_size+box_space)-(row['Z']-min_Z)*(box_space+box_size))+offset_y))+'"/>\n')
			title.append('    <text text-anchor="middle" font-family="'+str(font)+'" style="font-size:'+str(font_size)+'px; fill:'+str(font_colour)+'" x="'+str(int(offset_x+(N[i]-min_N)*(box_size+box_space))+box_size/2)+'" y="'+str(int(((max_Z-min_Z)*(box_size+box_space)-(row['Z']-min_Z)*(box_space+box_size))+offset_y+box_size/2))+'">'+str(naming.short_readable(int(row['Z']*10000000+row['A']*10000)))+'</text>\n')
		except ValueError:
			pass

	# print(new_lines)
	for j, line in enumerate(boxes):
		lines.insert(6+j, line)
	lines.insert(-1,'  <g id="layer1" fill="none">\n')
	insert_text_at = len(lines)-1
	for j, line in enumerate(title):
		lines.insert(insert_text_at+j, line)
	lines.insert(-1,'  </g>\n')
	lines.insert(-1,'  <g id="layer2" fill="none">\n')
	lines.insert(-1,'    <text text-anchor="middle" font-family="'+str(font)+'" style="font-size:'+str(font_size+30)+'px; fill:#000000" x="'+str(round((max_N-min_N)*(box_size+box_space)/2))+'" y="100" >Contributing Isomers with '+str(round(running_percent,2))+'% of total</text>\n')
	lines.insert(-1,'  </g>\n')    
	with open('.'+outdir+'/contributing_Isomers'+str(round(running_percent,2))+'.svg','w+') as f:
		f.writelines(lines)
