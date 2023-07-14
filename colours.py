import matplotlib.colors
import random

def makeRandomVisibleColourDict(size):
	cD = dict()
	for i in range(size):
	    randi = random.random()
	    hh = matplotlib.colors.hsv_to_rgb((randi, 1, 0.7))
	    hexi = matplotlib.colors.rgb2hex(hh)
	    cD[i] = hexi
	return (cD)
