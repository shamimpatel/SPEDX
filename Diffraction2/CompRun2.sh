g++ -Wall main2.cpp -o main2 || exit 1
echo
echo Compile complete
time ./main2
#gnuplot -persist "Plotscript.txt"