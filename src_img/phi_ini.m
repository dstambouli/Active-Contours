clear all
close all
clc

I = imread('3.pgm');

[X,Y] = meshgrid(1:size(I,2), 1:size(I,1));
phi = sqrt((X-50).^2 + (Y-50).^2) - 20;
figure;mesh(phi)

figure; imshow(I,[]); hold on;
contour(1:size(I,2), 1:size(I,1), phi, [0,0], 'r','Linewidth',2);
fileID = fopen('image_geometrique_level_set_function.bin','w');
fwrite(fileID,phi,'double');
fclose(fileID);
