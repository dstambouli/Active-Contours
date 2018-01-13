# Active-Contours
A project in parallel computing using MPI. 
Implementantion of "Active contours without edges" method of Chan &amp; Vese to detect objects in an image. You can find the original paper here: http://www.math.ucla.edu/~lvese/PAPERS/IEEEIP2001.pdf.

The method used here is slightly different as we modify the scheme to get an explicit one which is more efficient as it does not require computing the inverse of a matrix. 


# Elements to run the code
In main.c you'll need to modify Nlig and Mcol to match the size of your image.
Please note that in order to run the code you'll need:
  1. to transform your initial image into a .bin file using "ecriture_lecture_images.m" script
  2. modify in main.c the name of the file in "code = MPI_File_open(...)" with the name of the .bin you generate in 1
  3. your solution is stored in final_v2.dat and you can readin using the script "read_res.m"
 
 Compilation using make command and to run it on a single thread with ./CV, if you want to run on multiple threads which is the main purpose of this project use "mpirun -np nb_threads ./CV" where nb_threads is the number of threads you want to use.
 
 # Sorry the comments in the code are in french I'll try to translate them in english soon
