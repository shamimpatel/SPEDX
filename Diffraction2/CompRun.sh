g++ -Wall main.cpp -o main || exit 1
echo
echo Compile complete
time ./main
#gnuplot -persist "Plotscript.txt"