#!/usr/bin/env python
from __future__ import print_function, division
from optparse import OptionParser
import os
import sys
import matplotlib.pyplot as plt

__author__ = "Louis Dijkstra"

usage = """%prog [options] <histogram-file> <x-label>	Label for the x-axis. Spaces are allowed. 

	<histogram-file> 	Contains the histogram data. First
				column are the labels. Second
				column contains the counts.
	<x-label>		Label for the x-axis. Spaces are allowed. 

Plots the histogram data in the given file. The file should be 
organized in two column (tab-seperated): 

	x_1	c_1
	x_2	c_2
	...	...
	x_n	c_n

where x_1 is the minimal value found and x_n is the maximum value 
found. (Note: x_{i+1} = x_i + 1). c_i is the count for x_i. 
"""

class HistogramData:
	
	def __init__(self, histogram_filename, normalize=False):
		inputfile 	= open(histogram_filename, 'r')
		self.x_values 	= []
		self.count	= []
		self.max_count 	= float('-Inf')
		for line in inputfile:
			values = map(int, line.split())
			self.x_values.append(values[0])
			self.count.append(values[1])
			if self.max_count < values[1]:
				self.max_count = values[1]
		if normalize: 
			total = float(sum(self.count))
			if total != 0: 
				for i in range(len(self.count)):
					self.count[i] /= total  
				self.max_count /= total 

	def returnMinimumX(self):
		if self.x_values == []:
			return None
		return self.x_values[0]

	def returnMaximumX(self):
		if self.x_values == []:
			return None
		return self.x_values[-1]

	def extendValues(self, min_x, max_x):
		if min_x != self.x_values[0]:
			diff = self.x_values[0] - min_x 
			self.x_values = range(min_x, self.x_values[0]) + self.x_values
			self.count = [0] * diff + self.count
		if max_x != self.x_values[-1]:
			diff = max_x - self.x_values[-1]
			self.count += [0] * diff
			self.x_values += range(self.x_values[-1] + 1, max_x + 1) 

	def setToZero (self, min_x, max_x):
		self.x_values = range(min_x, max_x + 1)
		self.count = [0] * (max_x - min_x + 1)


def main():
	parser = OptionParser(usage=usage)
	parser.add_option("--normalize", action="store_true", dest="normalization", default=False,  
				help="Normalizes the data.")
	#parser.add_option("-k", action="store", dest="min_x", default=None, type=int,  
	#			help="Minimum x-value (Default: no minimum)")
	#parser.add_option("-l", action="store", dest="max_x", default=None, type=int,   
	#			help="Maximum x-value (Default: no maximum)")
	parser.add_option("-w", action="store", dest="width", default=1.0, type=float, 
					help="Bar width. (Default=1.0)")
	(options, args) = parser.parse_args()

	if (len(args)<2):
		parser.print_help()
		return 1

	histogram_filename = args[0]
	x_label		= ''
	for i in range(1, len(args) - 1):
		x_label += args[i] + ' '
	x_label += args[-1] 

	y_label = "# occurences"
	if options.normalization: 
		y_label = "probability"

	histogram_data = HistogramData(histogram_filename, normalize=options.normalization)

	plt.bar(histogram_data.x_values, histogram_data.count, align='center', width=options.width)
	plt.show()
	

if __name__ == '__main__':
	sys.exit(main())
