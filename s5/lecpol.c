#include <stdio.h>
//#include <termio.h>
#include <fcntl.h>

FILE *fichier;

ouvre_(nomfich)
char nomfich[];
{ fichier=fopen(nomfich, "r"); }

ferme_()
{ fclose(fichier); }

lecpol_(test,ntype,natt,i_att,nsom,poly)
char  *ntype;
int   *natt,*nsom,*test;
double i_att[];
float poly[];

{
int ib,ic;
char  fin;
lecture:

	*test=fscanf(fichier, "%c",ntype); if(*test==-1) return;

	fscanf(fichier, "%d", natt);
	for (ib=0; ib<*natt; ib++) fscanf(fichier, "%lf", &i_att[ib]);
	fscanf(fichier, "%d", nsom);
	ic=0;for (ib=0;ib<*nsom;ib++) {fscanf(fichier, "%f%f%f", &poly[ic],&poly[ic+1],&poly[ic+2]);ic+=3;}

	do { fscanf(fichier, "%c", &fin); } while(fin!='p' && !feof(fichier));

	if(!feof(fichier)) fseek(fichier,-sizeof(char),SEEK_CUR);
	
	/* 
	printf("fichier=%d \n", fichier);
	printf("ntype=%c natt=%d\n", *ntype, *natt);
	for (ib = 0; ib < *natt; ib++) printf("iatt(%d) = %ld \n",ib, i_att[ib]);
	printf("nsom=%d\n", *nsom);
	ic=0;for (ib=0;ib<*nsom;ib++) {printf("sommet %d %d  %f %f %f\n",ib,ic,poly[ic],&poly[ic+1],&poly[ic+2]);ic+=3;}
	*/
}








