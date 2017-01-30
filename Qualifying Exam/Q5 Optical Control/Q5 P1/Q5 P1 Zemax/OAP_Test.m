% This script opens an existing simple OAP Zemax model and changes
% parameters.
clear; close all; clc;

% Initialize Zemax comm
zDDEInit;
% Open desired file
zLoadFile('C:\Users\TyCoon\Google Drive\AFIT\PhD Research\Qualifying Exam\Q5 Optical Control\Q5 P1\Q5 P1 Zemax\SingleOAP.zmx')

% Copy lens in the ZEMAX DDE server into the Lens Data Editor (LDE)
zPushLens(1)

% Get the value of the mirror's semi-diameter
sD1 = zGetSurfaceData(3,5);

% Change the semi-diameter of the mirror
zSetSurfaceData(3,5,375)

% Get the value of the mirror's semi-diameter
sD2 = zGetSurfaceData(3,5);

