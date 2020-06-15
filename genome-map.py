#!/usr/bin/env python3
# Draw genome coverage map from mapped reads

# Usage
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Draw covered regions of genome',
	usage='\n  %(prog)s --svg map.svg <DEPTH-FILE>')
parser.add_argument('depth', metavar='DEPTH-FILE', help='samtools depth output file (required)')
parser.add_argument('--mindepth', metavar='INT', default=0, type=int, help='set threshold for min depth cutoff (default=0)')
parser.add_argument('--out', metavar='FILE', default='map.svg', help='save SVG output as specified file (default=map.svg)')
parser.add_argument('--size', metavar='WIDExHIGH', default='800x600', help='specify width and height of SVG in pixels (default="800x600")')
# parser.add_argument('--order', metavar='FILE', help='specify file containing list of genomes/chromosomes (1 per line) in desired order')
parser.add_argument('--minlen', metavar='LEN', default=0, type=int, help='specify minimum length of segment in bp to display (default=0)')
parser.add_argument('--colour', metavar='COLOUR', default='black', help='specify colour of recombination regions in HEX format (default=black)')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()

import os
import sys
import csv
import pandas as pd
from itertools import groupby
from operator import itemgetter
import svgwrite

# Functions
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

genomeLIST = []
nocovCOUNT = {}
covCOUNT = {}
genomeLEN = {}

# Import depth counts into pandas dataframe
df = pd.read_csv(args.depth, sep = '\t', names = ['chr', 'pos', 'cov'])

genomeDICT = df['chr'].unique()
for genome in genomeDICT:
	df_CHR = df[df['chr']==genome]
	df_COV = df_CHR[df_CHR['cov'] > args.mindepth]
	# Calculate number and percentage of covered sites
	total_sites = len(df_CHR)
	cov_sites = len(df_COV)
	nocov_sites = total_sites - cov_sites
	perc_cov = "{:.1%}".format(float(cov_sites/total_sites))
	covCOUNT[genome] = tuple([cov_sites, nocov_sites, total_sites, perc_cov])
	# Obtain list of consecutive covered sites to find start and end of covered segments
	loci_list = df_COV['pos'].tolist()
	for k, g in groupby(enumerate(loci_list), lambda ix: ix[0] - ix[1]):
		loci = list(map(itemgetter(1), g))
		genomeLIST.append(tuple([genome, min(loci), max(loci)]))
	genomeLEN[genome] = len(df_CHR)
# Find maximum length of genome
max_loci = max(genomeLEN.values())

# Draw SVG
svgsize = args.size.split('x',1)		# calibrate desired size of image
width = int(svgsize[0])			# width of image in pixels
height = int(svgsize[1])		# height of image in pixels
numseqs = len(genomeDICT)		# set number of sequences
h = height/float(numseqs)		# calibrate height based on number of sequences
s = width/float(int(max_loci))	# calibrate width to maximum sequence length
fsize = int(h*0.9)				# set font size for label
font = 'font-size:{}px; font-family:Arial'.format(fsize)	# set font
colour = args.colour			# set colour
main_colour = 'black'			# set colour
interval = 500000				# Tick interval in bp

dwg = svgwrite.Drawing(args.out)

def rect(x,p,w,colour):		# Draw regions
	dwg.add(dwg.rect(insert=((x*s)+5, (p*h)+5), size=(w, h*0.5), fill=colour))

def border(p,w):		# Box with 3px margin
	dwg.add(dwg.rect(insert=(0+5, (p*h)+5), size=(w, h*0.5), stroke='black', fill='none'))

def label(n,p):			# Sequence ID labels
	dwg.add(dwg.text(n, insert=((max_loci*s)+15, (((p+1)*h)-(0.2*h))+5), fill=main_colour, style=font))

def ticks(q):			# Draw tick marks
	count = interval
	while count < max_loci:
		dwg.add(dwg.line((count*s,((q)*h)), (count*s,((q)*h)-3), stroke=main_colour))
		dwg.add(dwg.text((count/1000000.0), insert=((count*s)-6, ((q)*h)+13), fill=main_colour, style=font))
		count = count + interval

def box(width,h):		# Box with 5px margin
	dwg.add(dwg.line((0,0), ((width+10),0), stroke=main_colour))
	dwg.add(dwg.line(((width+10),0), ((width+10),((p+1)*h)), stroke=main_colour))
	dwg.add(dwg.line(((width+10),((p+1)*h)), (0,((p+1)*h)), stroke=main_colour))
	dwg.add(dwg.line((0,((p+1)*h)), (0,0), stroke=main_colour))

# Draw covered regions and labels
msg('Drawing SVG to {} ...'.format(args.out))

# Draw a rectangle for each segment with coverage
for seg in genomeLIST:
	chr = seg[0]
	p = [*genomeDICT].index(chr)
	x = int(seg[1])		# locus start
	y = int(seg[2])		# locus stop
	w = (y-x)*s			# calculate width of bars
	# Set minimum segment length in bp
	if (y-x) > args.minlen:
		rect(x,p,w,colour)

# Draw a border for each chromosome
for chr in genomeDICT:
	p = [*genomeDICT].index(chr)
	border(p,genomeLEN[chr]*s)
	label(chr + ': ' + covCOUNT[chr][3],p)			# labels

# Draw bounding box
box(width,h)
# Add tick marks
ticks(p+1)
dwg.save()

# Exit
msg('Done.')
sys.exit(0)
