import multiprocessing
import subprocess

f = open('testoutput.dat','w')
f.write('Testfiles\n')

def worker(fname):
	subprocess.Popen([fname])

if __name__ == '__main__':
	files = ['./file1.py','./file2.py','./file3.py']
	for i in files:
		p = multiprocessing.Process(target=worker(i))
		p.start()