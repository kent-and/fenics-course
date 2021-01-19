# This script generates a presentation containing all the slides
# found in the slides subdirectory for easy inspection
#
# Anders Logg, 2009-06-15 -- 2009-06-16
# Modified by: Andre Massing 2012-06-24

# Clear and create temporary directory
rm -rf tmp
mkdir tmp

# Copy files
cp slides/*.tex tmp
cp *.sty tmp
cp *.cls tmp
cp -r eps tmp/
cp -r png tmp/
cp -r pdf tmp/
cp -r figures tmp/

# Enter tmp dir
cd tmp

# Generate slide collection
FILE="all.tmp"
rm -f $FILE
#echo "\\input{../preamble.tex}" >> $FILE
echo "\\documentclass{fenicscourse}" >> $FILE
echo "" >> $FILE
echo "\\begin{document}" >> $FILE
echo "" >> $FILE
for f in *.tex; do
    echo "{$f}"
    #Add file url to slide
    sed -i "s/end{frame}/color{red}\\\\url{slides\/$f}\n\\\\end{frame}/" $f
    echo "\\input{$f}" >> $FILE
done
echo "" >> $FILE
echo "\\end{document}" >> $FILE

# LaTeX file
mv all.tmp all.tex
pdflatex all

# View file
cp all.pdf ../all.pdf
cd ..
evince all.pdf
