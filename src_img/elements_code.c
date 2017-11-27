#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>




#define faux 0
#define ndims 2
#define NB_VOISINS 4
#define N 0
#define E 1
#define S 2
#define W 3
#define NE 0
#define ES 1
#define SO 2
#define ON 3
#define min(a,b) (a<=b?a:b)

/* rang dans le communicateur initial */
int rang;
/* nombre de processus */
int nb_procs;

/* nombre de lignes*/
/*METTRE LE NOMBRE DE LIGNES DE L'IMAGE*/
   int Nlig=100; /*direction des x*/
/* nombre de colonnes*/
/*METTRE LE NOMBRE DE LIGNES DE L'IMAGE*/
   int Mcol=100; /*direction des y*/

  /*********time discretization*******/

   double dt=0.1;

  /*********space discretization********/
   double  h=1.0;

   /********parameter epsilon in the regularized versions of H and delta********/
   double epsilon=1.0;

   /*******coefficients of the fidelity terms****************/
   double a1=1.0;
   double a2=1.0;

   /**********nombre d'itérations maximal*********************/
   int Itermax=400;

   /**********coefficient of the length term******************/
   double eta=0.1*255.0*255.0;


int main(int argc, char *argv[]) {

  /* coordonnées dans la grille */
  int coords[ndims];
  /* tableau des dimensions dans la grille */
  int dims[ndims];
  /* communicateur topologie cartésienne */
  MPI_Comm comm2d;
  int periods[ndims];
  const int reorganisation=faux;
  /* tableau contenant les voisins du sous-domaine courant (haut,bas,gauche,droite)*/
  int voisin[NB_VOISINS];

  /*nombre total de points intérieurs dans la direction x et la direction y*/
  int ntx, nty;
  /* ntx --> direction des lignes -->nombre de lignes total-2 (première ligne et dernière ligne puisqu'on ne fait
  pas de calcul sur ces lignes du fait des conditions au bord de type Neumann homogènes)*/

  /* nty --> direction des colonnes -->nombre de colonnes total-2 (première colonne et dernière colonne puisqu'on ne fait
  pas de calcul sur ces colonnes du fait des conditions au bord de type Neumann homogènes)*/

  int it;
  double t1, t2;




  void initialisation_mpi(int, char**);
  void finalisation_mpi();
  void domaine(MPI_Comm ,int ,int *,int,int ,int * ,int *,int * ,int ,int );
  void voisinage(MPI_Comm,int *,int *, int *);
  double **allocarray(int,int);
  void printarr(double **, int,int, char *);
  void ecrire_mpi(double *,int,int,int *,MPI_Comm);

  /*
        \
 ------- y                      coords[1]/dims[1]
 |      /                     /|\
 |                             |
 |                             |
\ /                            |
 x                             |       \
                               --------- coords[0]/dims[0]
                                       /
*/



  /* Initialisation de MPI */

  initialisation_mpi(argc,argv);

  /* Creation de la topologie cartesienne 2D */
  /*Le nombre de points dans la direction x correspond au nombre de lignes-2 (uniquement les points intérieurs) */
  /*Le nombre de points dans la direction y correspond au nombre de colonnes-2 (uniquement les points intérieurs) */
  ntx=Nlig-2;
  nty=Mcol-2;

  /* Connaître le nombre de processus selon x et le nombre de processus
     selon y en fonction du nombre total de processus */
  if(argc >= 2){
	dims[0] = (int)  *argv[1];
  }else{
  dims[0] =  0;
  }
  dims[1] = 0;
  MPI_Dims_create(nb_procs,ndims,dims);

  /* Creation de la grille de processus 2D sans periodicite */

  periods[0] = periods[1] = faux;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorganisation, &comm2d);
  MPI_Comm_rank(comm2d,&rang);

  /*if(rang == 0) {
    printf("Execution code crack_detection avec %d processus MPI\n"
	   "Taille du domaine : ntx=%d nty=%d\n"
	   "Dimension de la topologie : %d suivant y (colonnes), %d suivant x (lignes)\n"
	   "-----------------------------------------\n",
	   nb_procs, ntx, nty, dims[0], dims[1]);
  }*/

  int tab_bounds[4]; /*sx,ex,sy,ey*/
  /* Determinination des indices de chaque sous domaine */

  domaine(comm2d,rang,coords,ntx,nty,dims,tab_bounds,periods,reorganisation,nb_procs);
  /*printf("Je suis le rang %d\n"
	   "Ma coord 0 dans la direction des colonnes est: coord0=%d\n"
	   "Ma coord 1 dans la direction opposee a la direction des lignes est: coord1=%d\n"
	   "-----------------------------------------\n",
	   rang, coords[0], coords[1]);
  if(rang==3){
  printf("Je suis le rang %d \n"
  "coordonnee en ligne du coin superieur gauche sx=%d\n"
  "coordonnee en colonne du coin superieur gauche sy=%d\n"
  "coordonnee en ligne du coin inferieur droit ex=%d\n"
  "coordonnee en colonne du coin inferieur droit ey=%d\n"
  "-----------------------------------------\n",
	   rang,tab_bounds[0],tab_bounds[2],tab_bounds[1],tab_bounds[3]);
  }*/





  /* Recherche de ses 4 voisins pour chaque processus */
  voisinage(comm2d,voisin,coords,dims);


  int code;
  MPI_File descripteur;
  MPI_Offset deplacement_initial;
  MPI_Status statut;
  /* Ouverture du fichier "image_initiale_f.bin" en lecture */
  code = MPI_File_open(comm2d, "image_geometrique.bin", MPI_MODE_RDONLY,MPI_INFO_NULL, &descripteur);
  /* Test pour savoir si ouverture du fichier est correcte */
  if (code != MPI_SUCCESS) {
    fprintf(stderr, "ATTENTION erreur lors ouverture du fichier");
    MPI_Abort(comm2d, 2);
  }

  MPI_Datatype mysubarray;/****type sous-matrice****/


  double * nu_local=malloc((tab_bounds[1]-tab_bounds[0]+3)*(tab_bounds[3]-tab_bounds[2]+3)*sizeof(double));
  double ** nu_local_mat=allocarray((tab_bounds[1]-tab_bounds[0]+3),(tab_bounds[3]-tab_bounds[2]+3));

  if(nu_local==NULL)
    printf("Erreur dans l'allocation mémoire de nu_local-- \n");

  int starts[2] = {tab_bounds[0]-1,tab_bounds[2]-1};
  int subsizes[2]  = {tab_bounds[1]-tab_bounds[0]+3,tab_bounds[3]-tab_bounds[2]+3};
  int bigsizes[2]  = {Nlig, Mcol};


  MPI_Type_create_subarray(2,bigsizes, subsizes, starts,
                                 MPI_ORDER_C, MPI_DOUBLE, &mysubarray);
  MPI_Type_commit(&mysubarray);


  deplacement_initial=0;

  MPI_File_set_view(descripteur,deplacement_initial,MPI_DOUBLE,mysubarray,"native",MPI_INFO_NULL);
  MPI_File_read(descripteur,nu_local,(tab_bounds[1]-tab_bounds[0]+3)*(tab_bounds[3]-tab_bounds[2]+3),MPI_DOUBLE,&statut);

  MPI_File_close(&descripteur);


  for(int i=0;i<(tab_bounds[1]-tab_bounds[0]+3);i++){

     nu_local_mat[i]=&(nu_local[i*(tab_bounds[3]-tab_bounds[2]+3)]);

  }

  /* Définition de type dérivé type_ligne*/
  int nb_colonne_res = tab_bounds[3]-tab_bounds[2]+1;
  MPI_Type_contiguous(nb_colonne_res, MPI_DOUBLE, &type_ligne);
  MPI_Type_commit(&type_ligne);

  /* Définition de type dérivé type_colonne*/
  int nb_ligne_res = tab_bounds[1]-tab_bounds[0]+1;
  MPI_Type_vector(nb_ligne_res, 1, nb_colonne_res+2, MPI_DOUBLE, &type_colonne);
  MPI_Type_commit(&type_colonne);

  /* Initialisation d'un grand phi sur tout le domaine */
  if (rang == 0){
    phi_global = init_phi(Nlig, Mcol);
  }
  MPI_Barrier(comm2d);
  /* Allocation du phi_local */
  phi_local = init_phi(nb_ligne_res + 2, nb_colonne_res + 2);
  for(int i=0 ; i < nb_ligne_res + 2 ; i++){
    for(int j=0 ; i < nb_colonne_res + 2 ; j++){
     phi_local[i][j] = phi_global[tab_bounds[0] + i][tab_bounds[2] + j];
    }
  }
  


   /* Mesure du temps en seconde dans la boucle en temps */
   t1 = MPI_Wtime();

     /*BOUCLE PRINCIPALE*/


  /* Mesure du temps a la sortie de la boucle */
   t2 = MPI_Wtime();



  /* Ecriture des resultats pour chaque processus */
  /* FONCTION ecriture_mpi */

  /* Affichage du temps de convergence par le processus 3 */
  //if (rang == 3) {


  printf("Convergence en %f secs\n", t2-t1);
  //}
  /****Libération mémoire****/

  finalisation_mpi();

  return 0;
}



  /**************************************************************************************************/
  /*********************************INITIALISATION***************************************************/
  /**************************************************************************************************/
  void initialisation_mpi(int argc, char* argv[]) {
  /* Initialisation de MPI */
  MPI_Init(&argc, &argv);

  /* Savoir quel processus je suis */
  MPI_Comm_rank(MPI_COMM_WORLD, &rang);

  /* Connaitre le nombre total de processus */
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
 }

  /**************************************************************************************************/
  /*********************************FINALISATION***************************************************/
  /**************************************************************************************************/

  void finalisation_mpi() {
  /* Desactivation de MPI */
  MPI_Finalize();
  }


  /**************************************************************************************************/
  /*********************************CREATION DOMAINE*************************************************/
  /**************************************************************************************************/

  void domaine(MPI_Comm comm2d,int rang,int * coords,int ntx,int nty,int * dims,int * tab_bounds,int * periods,int reorganisation,int nb_procs) {



  int nx,ny,positionx,positiony,restex,restey;
  int sx,ex,sy,ey;

  /* Connaître mes coordonnees dans la topologie */
  MPI_Cart_coords(comm2d,rang,ndims,coords);


  /* Calcul pour chaque processus de ses indices de debut et de fin suivant x */

  /*Nombre de points dans la direction x*/

  nx=ntx/dims[1];
  restex=ntx % dims[1];
  positionx=coords[1];
  sx=1+positionx*nx+min(restex,positionx);/*Indice de depart dans la direction des x*/
  if(positionx<restex){
  nx=nx+1;
  }


  ex=sx+nx-1; /*Indice de fin dans la direction des x*/
  ex=min(ex,ntx+1);/*correction si besoin pour le dernier bloc*/


  /**********************************************************/
  /*Pour renumeroter selon la direction opposee a dims[1]   */
  /**********************************************************/



  sx=ntx-ex+1;
  ex=sx+nx-1;

  /*Nombre de points dans la direction y*/

  ny=nty/dims[0];
  restey=nty % dims[0];
  positiony=coords[0];
  sy=1+positiony*ny+min(restey,positiony);/*Indice de depart dans la direction des y*/
  if(positiony<restey){
  ny=ny+1;
  }

  ey=sy+ny-1;/*Indice de fin*/
  ey=min(ey,nty+1);
 /*correction si besoin pour le dernier bloc*/

  tab_bounds[0]=sx;
  tab_bounds[1]=ex;
  tab_bounds[2]=sy;
  tab_bounds[3]=ey;


  }

  /**************************************************************************************************/
  /*********************************VOISINAGE********************************************************/
  /**************************************************************************************************/

  void voisinage(MPI_Comm comm2d,int * voisin,int * coords, int * dims) {


  /* Recherche des voisins Nord et Sud */
  MPI_Cart_shift(comm2d, 0, 1, &(voisin[W]), &(voisin[E]));

  /* Recherche des voisins Ouest et Est */
  MPI_Cart_shift(comm2d, 1, 1, &(voisin[S]), &(voisin[N]));



  }



   double **allocarray(int Nlig,int Mcol) {
    double **array2 = malloc( Nlig* sizeof( double * ) );
    int i;

    if( array2 != NULL ){
        array2[0] = malloc(Nlig * Mcol * sizeof( double ) );
        if( array2[ 0 ] != NULL ) {
            for( i = 1; i < Nlig; i++ )
                array2[i] = array2[0] + i * Mcol;
        }

        else {
            free(array2);
            array2 = NULL;
            printf("Erreur dans l'allocation mémoire -- phase 2\n");
        }
    }
    return array2;
    }

   void printarr(double **data, int nlig,int mcol, char *str) {
    printf("-- %s --\n", str);
    for (int i=0; i<nlig; i++) {
        for (int j=0; j<mcol; j++) {
            printf("%3f ", data[i][j]);
        }
        printf("\n");
    }
   }



 void ecrire_mpi(double *v2_local_vect,int ntx,int nty,int * tab_bounds,MPI_Comm comm2d){

  int code;
  MPI_File descripteur;
  int profil_tab[ndims], profil_sous_tab[ndims], coord_debut[ndims];
  MPI_Datatype type_sous_tab, type_sous_tab_vue;
  int profil_tab_vue[ndims], profil_sous_tab_vue[ndims], coord_debut_vue[ndims];
  MPI_Offset deplacement_initial;
  MPI_Status statut;

  /* Ouverture du fichier "final_v2.dat" en écriture */
  code = MPI_File_open(comm2d, "final_v2.dat", MPI_MODE_WRONLY+MPI_MODE_CREATE,
		MPI_INFO_NULL, &descripteur);

 /* Test pour savoir si ouverture du fichier est correcte */
  if (code != MPI_SUCCESS) {
    fprintf(stderr, "ATTENTION erreur lors ouverture du fichier");
    MPI_Abort(comm2d, 2);
  }

  /* Creation du type derive type_sous_tab qui definit la matrice
   * sans les cellules fantomes */
  profil_tab[0] = tab_bounds[1]-tab_bounds[0]+3;
  profil_tab[1] = tab_bounds[3]-tab_bounds[2]+3;

  /* Profil du sous tableau */
  profil_sous_tab[0] =tab_bounds[1]-tab_bounds[0]+1;
  profil_sous_tab[1] =tab_bounds[3]-tab_bounds[2]+1;

  /* Coordonnees de depart du sous tableau */
  coord_debut[0] = 1;
  coord_debut[1] = 1;

  /* Creation du type_derive type_sous_tab */
  MPI_Type_create_subarray(ndims, profil_tab, profil_sous_tab, coord_debut,
			    MPI_ORDER_C, MPI_DOUBLE, &type_sous_tab);

  /* Validation du type_derive type_sous_tab */
  MPI_Type_commit(&type_sous_tab);

  /* Creation du type type_sous_tab_vue  pour la vue sur le fichier */
  /* Profil du tableau global */
  /*On ne tient pas compte des bords de l'image pour plus de simplicité*/
  profil_tab_vue[0] = ntx;
  profil_tab_vue[1] = nty;

  /* Profil du sous tableau */
  profil_sous_tab_vue[0] = tab_bounds[1]-tab_bounds[0]+1;
  profil_sous_tab_vue[1] = tab_bounds[3]-tab_bounds[2]+1;

  /* Coordonnees de depart du sous tableau */
  coord_debut_vue[0] = tab_bounds[0]-1;
  coord_debut_vue[1] = tab_bounds[2]-1;

  /* Creation du type_derive type_sous_tab_vue */
  MPI_Type_create_subarray(ndims, profil_tab_vue, profil_sous_tab_vue, coord_debut_vue,
			   MPI_ORDER_C, MPI_DOUBLE, &type_sous_tab_vue);

  /* Validation du type_derive type_sous_tab_vue */
  MPI_Type_commit(&type_sous_tab_vue);

  /* Définition de la vue sur le fichier a partir du debut */
  deplacement_initial = 0;
  MPI_File_set_view(descripteur, deplacement_initial, MPI_DOUBLE,
		    type_sous_tab_vue, "native", MPI_INFO_NULL);

  /* Ecriture du tableau u par tous les processus avec la vue */
  MPI_File_write_all(descripteur, v2_local_vect, 1, type_sous_tab, &statut);

  /* Fermeture du fichier */
  MPI_File_close(&descripteur);
}

double **init_phi(int Nlig,int Mcol) {
 double **array2 = malloc( Nlig* sizeof( double * ) );
 int i;

  if( array2 != NULL ){
    for( i = 0; i < Nlig; i++ ){
      array2[i] = malloc(Mcol*sizeof( double ));
      for( int j = 0; j < Mcol; j++){
        /* On intialise un carre */
        if (i == Nlig -1 || i == 0 || j == Mcol -1 || j == 0){
          array2[i][j] = 1.0;
        }
        else{
            array2[i][j] = 0.0;
        }
      }
    }
  }
  else {
     free(array2);
     array2 = NULL;
     printf("Erreur dans l'allocation mémoire -- phase 2\n");
  }

 return array2;
 }
