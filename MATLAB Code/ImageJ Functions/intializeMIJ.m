% Set plot defaults
set(0,'defaultAxesTickDir','out')
set(0,'defaultAxesTickDirMode','manual')
set(0,'defaultAxesBox','off')
set(0,'DefaultFigureWindowStyle','normal') %docked or normal
set(0,'DefaultFigureColor','w')
warning('off','images:imshow:magnificationMustBeFitForDockedFigure');

% Add ImageJ to working directory
javaaddpath([matlabroot '\java\mij.jar']);
javaaddpath([matlabroot '\java\ij.jar']);

% Add ImageJ plugins to the current path
fijiPath = 'C:\Fiji.app\';
javaaddpath([fijiPath 'plugins\plugins'])
javaaddpath([fijiPath 'plugins\plugins\bij.jar'])
javaaddpath([fijiPath 'plugins\Cell_Magic_Wand_Tool.jar'])
javaaddpath([fijiPath 'plugins\Image_Stabilizer\'])
javaaddpath([fijiPath 'plugins\bUnwarpJ_-2.6.2.jar']);
addpath([fijiPath 'scripts\']);

% Startup ImageJ 
Miji