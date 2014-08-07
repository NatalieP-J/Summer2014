import multiprocessing
import subprocess
from numpy import *
import os

def makeWMfiles(rangeval):
	for i in range(rangeval):
		f = open('WMrateget{0}.py'.format(i+1),'wb')
		f.write('i = {0}\n'.format(i))
		f.close()
		os.system('cat WMrateget.py >> WMrateget{0}.py'.format(i+1))

rangeval = 7

makeWMfiles(7)

def worker(i):
	subprocess.Popen(['python','WMrateget{0}.py'.format(i+1)])

if __name__ == '__main__':
	for i in range(rangeval):
		p = multiprocessing.Process(target=worker(i))
		p.start()
