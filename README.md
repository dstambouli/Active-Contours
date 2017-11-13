# Active-Contours
A project in parallel computing using MPI. 
Implementantion of "Active contours without edges" method of Chan &amp; Vese to detect objects in an image. You can find the original paper here: http://www.math.ucla.edu/~lvese/PAPERS/IEEEIP2001.pdf.

The method used here is slightly different as we modify the scheme to get an axplicit one which is more efficient as it does not require computing the invers of a matrix
