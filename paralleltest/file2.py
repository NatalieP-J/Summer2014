import time
for i in range(100):
	f = open('testoutput.dat','a')
	f.write('This is file 2 at {0}\n'.format(time.time()))