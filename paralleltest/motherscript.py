import multiprocessing
import subprocess

f = open('testoutput.dat','w')
f.write('Testfiles\n')

def worker(fname):
	subprocess.Popen(['python',fname])

if __name__ == '__main__':
	files = ['file1.py','file2.py','file3.py','file4.py','file5.py','file6.py','file7.py','file8.py']
	for i in files:
		p = multiprocessing.Process(target=worker(i))
		p.start()
