#!/usr/bin/python
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import sys
import ast
import warnings

warnings.filterwarnings("ignore")

DEFAULT_COLOR = "b"
DEFAULT_ALPHA = 0.4
X_LABEL = "M"
Y_LABEL = "N"
Z_LABEL = "H"

def int_list(arg):
	return ast.literal_eval(arg)

def coord_list(arg):
	return ast.literal_eval(arg)

def plot_2D(k, n, m, mu, lengths, widths, color, alpha):
	fig = plt.figure()
	ax = fig.add_subplot(111, aspect="equal")
	ax.set_xlim(0, m)
	ax.set_ylim(0, n)
	ax.set_xlabel(X_LABEL)
	ax.set_ylabel(Y_LABEL)
	for i in range(k):
		length = lengths[i]
		width = widths[i]
		a = mu[i][0]
		b = mu[i][1]
		ax.add_patch(patches.Rectangle((a, b), length, width, color=color, alpha=alpha))
	plt.show()

def plot_3D(k, n, m, h, mu, lengths, widths, heights, color, alpha):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_xlim(0, m)
	ax.set_ylim(0, n)
	ax.set_zlim(0, h)
	ax.set_xlabel(X_LABEL)
	ax.set_ylabel(Y_LABEL)
	ax.set_zlabel(Z_LABEL)
	for i in range(k):
		length = lengths[i]
		width = widths[i]
		height = heights[i]
		a = mu[i][0]
		b = mu[i][1]
		c = mu[i][2]
		rx = [a, a + length]
		ry = [b, b + width]
		rz = [c, c + height]
		X1, Y1 = np.meshgrid(rx, ry)
		X2, Z2 = np.meshgrid(rx, rz)
		Y3, Z3 = np.meshgrid(ry, rz)
		ax.plot_surface(X1, Y1, c, color=color, alpha=alpha)
		ax.plot_surface(X1, Y1, c + height, color=color, alpha=alpha)
		ax.plot_surface(X2, b, Z2, color=color, alpha=alpha)
		ax.plot_surface(X2, b + width, Z2, color=color, alpha=alpha)
		ax.plot_surface(a, Y3, Z3, color=color, alpha=alpha)
		ax.plot_surface(a + length, Y3, Z3, color=color, alpha=alpha)
	plt.show()

def main():
	parser = argparse.ArgumentParser(description="Plot 2D and 3D rectangles")
	parser.add_argument("k", type=int)
	parser.add_argument("n", type=int)
	parser.add_argument("m", type=int)
	parser.add_argument("h", type=int)
	parser.add_argument("mu", type=coord_list)
	parser.add_argument("lengths", type=int_list)
	parser.add_argument("widths", type=int_list)
	parser.add_argument("heights", type=int_list, nargs="?")
	parser.add_argument("-c", "--color", type=str, nargs="?")
	parser.add_argument("-a", "--alpha", type=float, nargs="?")


	args = parser.parse_args(sys.argv[1:])
	if len(args.mu) != args.k : 
		raise argparse.ArgumentTypeError("Size of rectangles positions and k not matching")
	if len(args.lengths) != args.k :
		raise argparse.ArgumentTypeError("Size of lengths array and k not matching")
	if len(args.widths) != args.k :
		raise argparse.ArgumentTypeError("Size of widths array and k not matching")
	if args.h > 0 and len(args.heights) != args.k :
		raise argparse.ArgumentTypeError("Size of heights array and k not matching")
	color = args.color if args.color else DEFAULT_COLOR
	alpha = args.alpha if args.alpha else DEFAULT_ALPHA
	if args.h > 0 :
		plot_3D(args.k, args.n, args.m, args.h, args.mu, args.lengths, args.widths, args.heights, color, alpha)
	else :
		plot_2D(args.k, args.n, args.m, args.mu, args.lengths, args.widths, color, alpha)



if __name__ == "__main__":
	main()