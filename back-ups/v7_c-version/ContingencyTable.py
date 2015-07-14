#!/usr/bin/env python
from __future__ import print_function, division

__author__ = "Louis Dijkstra"

"""
	Classes for working with nxn contingency tables.
"""

class Table: 

	def __init__(self, n_rows, n_columns):
		self.n_rows = n_rows
		self.n_columns = n_columns
		self.table = [[0 for i in range(n_columns)] for j in range(n_rows)] 
		self.total = 0 
		
	def add(self, row,column):
		self.table[row][column] += 1
		self.total += 1

	def print(self):
		print(self.table) 
	
class SquareTable(Table):
	
	def __init__(self, n_classes):
		Table.__init__(self, n_classes, n_classes)

	def returnPercentageCorrect(self):
		total_diagonal = 0
		if self.total == 0:
			return '-'
		for i in range(self.n_rows):
			total_diagonal += self.table[i][i]
		
		return float(total_diagonal) / float(self.total)

	def print(self):
		print(self.table) 

class Table3x3(SquareTable):
	
	def __init__(self):
		SquareTable.__init__(self, 3)

	def print(self, title = "", labels = ["0.0", "0.5", "1.0"]):
		"""Prints a 3x3 table to standard output"""
		row_total, column_total, total = [0,0,0], [0,0,0], 0
		# compute the marignals and the total
		for r in [0,1,2]:
			for c in [0,1,2]:		
				row_total[r] 	+= self.table[r][c]
				column_total[c] += self.table[r][c]
				total 		+= self.table[r][c]

		n_correct = 0 # number of correctly classified item
		for i in [0,1,2]: 
			n_correct += self.table[i][i] # on the diagonal

		print('------------------------' + title + '------------------------')
		print('\t\t\t\t\t\ttruth')
		print('\t\t\t|\t', labels[0], '\t|\t', labels[1], '\t|\t', labels[2], '\t|\ttotal')
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\t', labels[0], '\t|\t',self.table[0][0], '\t|\t', self.table[0][1],'\t|\t', self.table[0][2],'\t|\t', row_total[0])
		print('assign.\t\t', labels[1], '\t|\t',self.table[1][0], '\t|\t', self.table[1][1],'\t|\t', self.table[1][2],'\t|\t', row_total[1])
		print('\t\t', labels[2], '\t|\t',self.table[2][0], '\t|\t', self.table[2][1],'\t|\t', self.table[2][2],'\t|\t', row_total[2])
		print('\t--------------------------------------------------------------------------------------------')
		print('\t\ttotal\t|\t', column_total[0], '\t|\t', column_total[1], '\t|\t', column_total[2], '\t|\t', total)
		if total != 0:
			print('\ncorrectly classified:\t', n_correct, ' /', n_correct / float(total) * 100.0, '%\n')
		else: 
			print('\n')

class Table2x2(SquareTable):

	def __init__(self):
		SquareTable.__init__(self, 2)

	def returnRecallPrecision(self): 
		"""Returns the recall and precision. Returns 'None' when not defined"""
		recall, precision = None, None
		column_total	= self.table[0][0] + self.table[1][0]
		row_total 	= self.table[0][0] + self.table[0][1]
		if column_total != 0:
			recall = float(self.table[0][0]) / float(column_total)
		if row_total != 0:
			precision = float(self.table[0][0]) / float(row_total)
		return recall, precision

	def print(self, title = "", labels = ["somatic", "not somatic"]):
		"""Prints a 2x2 table to standard output"""
		row_total, column_total, total = [0,0], [0,0], 0
		# compute the marignals and the total
		for r in [0,1]:
			for c in [0,1]:		
				row_total[r] 	+= self.table[r][c]
				column_total[c] += self.table[r][c]
				total 		+= self.table[r][c]
		recall, precision = self.returnRecallPrecision()
		print('------------------------' + title + '------------------------')
		print('\t\t\t\t\ttruth')
		print('\t\t|\t', labels[0], '\t\t|\t', labels[1], '\t|\ttotal')
		print('--------------------------------------------------------------------------------------------')
		print(labels[0], '\t\t|\t', self.table[0][0], '\t\t|\t', self.table[0][1], '\t\t|\t', row_total[0])
		print(labels[1], '\t\t|\t', self.table[1][0], '\t\t|\t', self.table[1][1], '\t\t|\t', row_total[1])	
		print('--------------------------------------------------------------------------------------------')
		print('total\t\t|\t', column_total[0], '\t\t|\t', column_total[1], '\t\t|\t', total)
		print('\nrecall:\t\t', recall)
		print('precision:\t', precision)	
		print('-------------------------------------------------------------------\n\n')
		

