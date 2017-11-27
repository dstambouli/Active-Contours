clear all
close all
clc

%%%%%%%%%%%%%Dimensions 100x100%%%%%%%%%%%%%%%%%

I=imread('3.pgm');
fileID = fopen('image_geometrique.bin','w');
fwrite(fileID,I','double');
fclose(fileID);


fid=fopen('image_geometrique.bin','r','l');

%Opération de dévectorisation pour revenir à un objet
%matrice 

mydata=fread(fid,'double');
for(i=1:100)
    for(j=1:100)
        A(i,j)=mydata((i-1)*100+j);
    end
end
figure;imshow(A,[]);


clear all
close all
clc

%%%%%%%%%%%%%Dimensions 180x150%%%%%%%%%%%%%%%%%

I=imread('image1.pgm');
fileID = fopen('image_cerveau1.bin','w');
fwrite(fileID,I','double');
fclose(fileID);


fid=fopen('image_cerveau1.bin','r','l');
mydata=fread(fid,'double');
for(i=1:180)
    for(j=1:150)
        A(i,j)=mydata((i-1)*150+j);
    end
end
figure;imshow(A,[]);


clear all
close all
clc

%%%%%%%%%%%%%Dimensions 180x150%%%%%%%%%%%%%%%%%

I=imread('left_sag_069.tif');
fileID = fopen('image_cerveau2.bin','w');
fwrite(fileID,I','double');
fclose(fileID);


fid=fopen('image_cerveau2.bin','r','l');
mydata=fread(fid,'double');
for(i=1:128)
    for(j=1:192)
        A(i,j)=mydata((i-1)*192+j);
    end
end
figure;imshow(A,[]);
