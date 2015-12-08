datFile = '/media/win/Users/jaymz/Desktop/new_data/reserviorSparklingWater-7-24-15_highTemps_dataDump.csv'

dat = []

file = open(datFile,'r')

for line in file:
    l = str.split(line,',')
    ln = []
    for thing in l:
        ln.append(float(thing))
    dat.append(ln)
    
dat = np.matrix(dat)

data = {'x' : dat[:,0],
        'x(wn)' : float(10**7)/dat[:,0],
        'y' : dat[:,1], 
        'ringdownTime' : -1.0/dat[:,1] }