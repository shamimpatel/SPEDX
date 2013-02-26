cd DiffractionProbabilities
"bin/Release/DiffractionProbabilities.exe"
python ProcessScatterResults.py
xcopy ProcessedResults.txt "../Diffraction2" /Y
cd "../Diffraction2"
bin\Release\Diffraction2.exe
gnuplot.exe plotscript2.txt
start plot2.pdf
xcopy AdvDiffractResults.txt "../PinholeDiffraction" /Y
xcopy AdvFluoResults.txt "../PinholeDiffraction" /Y
cd "../PinholeDiffraction"
bin\Release\PinholeDiffraction.exe
gnuplot.exe plotscript2.txt
start plot2.pdf
python ProduceSpectrumPostPinhole.py
PlotSpectrum.bat
pause