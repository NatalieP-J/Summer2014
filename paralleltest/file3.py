import time
#for i in range(100):
while True:
	f = open('testoutput.dat','a')
	f.write('This is file 3 at {0}\n'.format(time.time()))
