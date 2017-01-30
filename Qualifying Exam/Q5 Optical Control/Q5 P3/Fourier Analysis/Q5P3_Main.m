% Tim Coon
% Qualifying Exam Q5 - Optical Control
% Characterize the influence of a surface roughness term on the focal point
% distribution.

% Plot frauhofer diffraction pattern in the focal plane for the single
% parabolic mirror. The input is a collimated beam and we may reference the
% scenario of an input placed against a lens as detailed in Goodman section
% 5.2.1 validating the math to show the diffraction pattern in the focal
% plane is approximated by the Fraunhofer pattern applied to the input
% field.

clear; close all; clc;

% Plot the Fraunhofer diffraction pattern for an ideal parabolic mirror
% focusing a collimated input beam along the mirror principal axis
run = 1;
isCentered = true;
isIdeal = true;
% Q5P3_diffPat
SM_diff(run,isCentered,isIdeal)

clear
% Plot the Fraunhofer diffraction pattern for an ideal parabolic mirror
% focusing a collimated input beam offset from the mirror principal axis.
% Essentially, this is an off-axis parabolic mirror.
run = 2;
isCentered = false;
isIdeal = true;
% Q5P3_diffPat
SM_diff(run,isCentered,isIdeal)

clear
% Plot the Fraunhofer diffraction pattern for an ideal parabolic mirror
% focusing a collimated input beam offset from the mirror principal axis.
% Essentially, this is an off-axis parabolic mirror.
run = 3;
isCentered = false;
isIdeal = false;
% Q5P3_diffPat
SM_diff(run,isCentered,isIdeal)