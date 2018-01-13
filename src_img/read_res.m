clear all
close all
clc

%%%%%%%%%%%%%Dimensions 100x100%%%%%%%%%%%%%%%%%

I=imread('world_map.png');

fid=fopen('final_v2.dat','r','l');
mydata = fread(fid,'double');

%Opération de dévectorisation pour revenir à un objet
%matrice
nlig = 345;
ncol = 460;
A = ones(nlig,ncol)*254;


for(i=2:(nlig -1))
    for(j=2:(ncol-1))
        A(i,j)=mydata((i-2)*(ncol-2)+(j-1));
    end
end
figure;imshow(I,[]); hold on;
contour(1:size(A,2) ,1:size(A,1),A,[0,0],'r','Linewidth',2)
