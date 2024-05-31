#entering the specific directory
cd Documents/PDC/Project

# compiling the code
mpicc -o Project Project.c -fopenmp

echo -e "Press 1 and enter when asked"
#mpiexec -n <num of processes> <filename>
mpiexec -n 2 ./Project doctorwho.csv

sleep 5
echo -e "Press 2 and enter when asked"
#sleep function so that the whole code can work
#mpiexec -n <num of processes> <filename>
mpirun -np 2 ./Project Email-EuAll.txt 265214

sleep 5
echo -e "Press 3 and enter when asked"
#sleep function so that the whole code can work
#mpiexec -n <num of processes> <filename>
mpirun -np 2 ./Project Email-Enron.txt 36692



