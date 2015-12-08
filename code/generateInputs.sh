echo "Shell script for example usage of the extract.py script"
python extract.py -o test/H2O -m H2O -x 1597,1603 -X nm ../../../hitran-master/hitran-master/par/01_hit12.par
python extract.py -o test/CO2-C12 -i 1 -m CO2 -x 1597,1603 -X nm ../../../hitran-master/hitran-master/par/02_hit12.par
python extract.py -o test/CO2-C13 -i 2 -m CO2 -x 1597,1603 -X nm ../../../hitran-master/hitran-master/par/02_hit12.par
