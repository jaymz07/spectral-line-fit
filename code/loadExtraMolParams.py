moleculeFile = open('molparam.txt')

moles = {}
mol = None
molData = {}

for line in moleculeFile:
    split = line.split()
    if( len(split) == 2):
        if(len(molData)>0):
            moles[mol] = molData
        mol = split[0]
        molData = {}
    if(len(split) == 5):
        molData[ split[0] ] =   { 'abundance' :  float(split[1]),
                                'Q(296K)' :    float(split[2]),
                                'gj' :         float(split[3]),
                                'mass' :       float(split[4])     }

moles[mol] = molData