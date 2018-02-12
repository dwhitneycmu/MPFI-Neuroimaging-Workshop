Steps to getting MIJ to work with MATLAB:
1.) Download and install Fiji here: https://imagej.net/Fiji/Downloads
2.) In the "MATLAB-Required Files" folder, copy mij.jar and ij.jar to the "java" subfolder of your MATLAB root path. 
Examples: Copy the mij.jar to '/Applications/MATLAB_R2017a.app/java/mij.jar'
Examples: Copy the ij.jar to '/Applications/MATLAB_R2017a.app/java/ij.jar'
If you are unsure what your MATLAB root path is, you can find it out with the "matlabroot" function in MATLAB.
3.) Next increase the Java Heap space allowed for Fiji in MATLAB. First, go to your Preferences-->MATLAB-->General-->Java Heap Memory. Increase the Java Heap Size by sliding the slider to the right. You may note that the maximal allowed RAM for the Java Heap Memory is 25%. If you want to increase the maximal RAM further, then copy the java.opts from the "MATLAB-Required Files" into the "bin\win64" subfolder of your MATLAB root path. Next open the java.opts file, and edit the maximal allowed memory. The file should have a single line: "-Xmx4084m". The 4096 corresponds to 4GB=4096KB. Change this number to the desired amount of RAM you want to allocate to your Java Heap Memory.
4.) Finally to get Fiji to run in MATLAB, we need to add the mij.jar and ij.jar to MATLAB's java classpath:
Type: javaaddpath '/Applications/MATLAB_R2017a.app/java/mij.jar'
Type: javaaddpath '/Applications/MATLAB_R2017a.app/java/ij.jar'
If you don't want to type these commands everytime MATLAB launches, add them to your startup.m file.
5.) Congratulations you're done! Start Fiji in MATLAB with "MIJ.start".

Steps to getting Theo Walker's Cell Magic Wand installed:
1.) Copy the Cell_Magic_Wand_Tool.jar from the "ImageJ Plugins" folder into the "plugins\Tools" subfolder of your ImageJ root path. You may need to create the "Tools" directory.
Example: '/Applications/Fiji.app/plugins/Tools/Cell_Magic_Wand_Tool.jar'
2.) Cell Magic wand should now either appear in your ImageJ toolbar and/or be listed in the Plugins menu. 
3.) Optionally you can directly download the Cell Magic Wand Tool at https://www.maxplanckflorida.org/fitzpatricklab/software/cellMagicWand/

P.S. Daniel Sage offers a fantastic writeup for getting MIJ installed: http://bigwww.epfl.ch/sage/soft/mij/