import multiprocessing
import subprocess
from numpy import *

def worker(i):
	subprocess.Popen(['python','WMrateget{0}.py'.format(i+1)])

if __name__ == '__main__':
	for i in range(8):
		p = multiprocessing.Process(target=worker(i))
		p.start()
