#!/usr/bin/env python
from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np

__author__ = "Louis Dijkstra"

"""
	Contains the functionality needed to plot the results of a ternary classification task,
	e.g., somatic/germline/not present or absent/heterozygous/homozygous.
"""

def normalizeTable(table):
	"""Normalizes the rows in a 3x3 table, i.e., every row sums up to 100%."""
	for row in [0,1,2]:
		row_total = 0.0
		for column in [0,1,2]: # compute row total
			row_total += table[row][column]
		if row_total != 0:
			for column in [0,1,2]: # normalize
				table[row][column] = table[row][column] / row_total * 100.0
	return table

def returnFancyIntervalString(interval):
	"""Returns the interval in the form of a string. Suitable for plotting."""
	if interval[0] is None:
		if interval[1] is None: # unbounded
			return r".-."
		else:	# <= interval[1]
			return r"$\leq$" + str(interval[1])
	elif interval[1] is None: # >= interval[0]
			return r"$\geq$" + str(interval[0])
	return str(interval[0]) + '-'  + str(interval[1])


def plotTernaryClassification (class1, class2, class3, length_ranges, class_names = ['not present', 'heterozygous', 'homozygous'], width = .9):
	layer1, layer2, layer3, labels = [], [], [], []

	for i in range(len(length_ranges)):
		layer1.append(class1[i][0])
		layer2.append(class1[i][1])
		layer3.append(class1[i][2])
		labels.append(returnFancyIntervalString(length_ranges[i]))
	for i in [0]:
		layer1.append(0)
		layer2.append(0)
		layer3.append(0)
		labels.append('')
	for i in range(len(length_ranges)):
		layer1.append(class2[i][0])
		layer2.append(class2[i][1])
		layer3.append(class2[i][2])
		labels.append(returnFancyIntervalString(length_ranges[i]))
	for i in [0]:
		layer1.append(0)
		layer2.append(0)
		layer3.append(0)
		labels.append('')
	for i in range(len(length_ranges)):
		layer1.append(class3[i][0])
		layer2.append(class3[i][1])
		layer3.append(class3[i][2])
		labels.append(returnFancyIntervalString(length_ranges[i]))

	ind 	= np.arange(len(labels))
	bottom 	= layer1
	plot_layer1 = plt.bar(ind, layer1, width, color = 'b')
	plot_layer2 = plt.bar(ind, layer2, width, color = 'r', bottom = bottom)
	for i in range(len(bottom)):
		bottom[i] += layer2[i]
	plot_layer3 = plt.bar(ind, layer3, width, color = 'g', bottom = bottom)

	x_min,x_max,y_min,y_max = plt.axis()
	plt.axis((x_min, x_max, 0, 100))

	plt.xticks(ind + width / 2.0, labels)

	plt.yticks(np.arange(0,101,20), ('0%', '20%', '40%', '60%', '80%', '100%'))
	plt.grid(True, axis='y')

	plt.legend( (plot_layer1[0], plot_layer2[0], plot_layer3[0]), class_names, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
	plt.show()
