c	S5 CALCULE LA REPARTITION DANS UNE GRILLE 3D DES CARACTERISTIQUES D'UNE MAQUETTE INFORMATIQUE
c	DENSITE VOLUMIQUE DES SURFACES DES TRIANGLES
c	DISTRIBUTION D'ORIENTATION DES NORMALES AUX TRIANGLES
c	DENSITE VOLUMIQUE DES VARIABLES SURFACIQUES ASSOCIEES (ATTRIBUTS).
	
c	AUTEUR:
c	Bruno Andrieu 
c	INRA Bioclimatologie 78850 Thiverval-Grignon
c	tel #33 1 30815527  Fax #33 1 30815527
c	Email andrieu@bcgn.inra.fr

c	nb1: la maquette peut comprendre plusieurs especes vegetales.
c	Dans ce cas l'analyse est faite pour chacune des especes presentes.
c	La variable espece est codee via le label associ� a chaque triangle
c	(espece= label/10**11)

c	nb2: La structure peut etre de type periodique dans le plan horizontal.
c	(typiquement une periode = un interrang).
c	Une maquette periodique est constituee de la repetition d'un motif elementaire.
c	Dans ce cas l'ensemble des triangles est analyse pour calculer les variations
c	moyennes a l'interieur du motif. Dans le plan horizontal, le motif est divise en
c	njx * njy cellules de dimension dx et dy.
c	(njx*dx et njy*dy representent la periode en x et y) 

c            ****              ****              ****              ****
c           ******            ******            ******            ******
c          ********          ********          ********          ********
c             **                **                **                **              unite de longueur
c             **                **                **                **              >-----<   
c             **                **                **                **
c                                ! ! ! ! ! ! ! ! ! !                         
c                               >------------------<                          ! ! dy = 0.4   et  njy = 9
c                                 longueur =  njy*dy   
c     
c              >-----------------------------------------------------<            
c                  longueur =  yl
c
c
c	nb3:  Le label est decode pour determiner si le triangle lu appartient a une
c	feuille (element surfacique) ou a une tige (element volumique). Dans ce dernier cas, la
c	surface du triangle est ponderee par 0.5 pour les calculs de densite de surface foliaire.
c	Ainsi la surface foliaire associee a un volume est la moitie de sa surface developpee.
c	Ceci est coherent avec le fait que la surface prise en compte pour une feuille est calculee
c	pour une seule face de la feuille.
c
c	nb4: Pour ce qui est de la densite volumique des attributs (variables associees aux triangles):
c	On suppose que les attributs representent des donnees surfaciques. Si les triangles 
c	representent la surface englobante d'un volume, le traitement d'attributs relatifs a ce volume
c	ne peut se faire que si ces attributs sont presentes comme des variables surfaciques 
c	reparties sur la surface englobante.
c
c	nb5: les dimensions maximales des tableaux sont definies dans le code
c	par des commandes parameter. Modifier ces lignes si necessaire pour augmenter 
c	le nombre de cellules ou de classes d'angles.


c	Donnees en entree:
c	******************

c	Fichier fort.51 contenant les triangles  au format can
c	(verifier :lecture non formatee ou format libre)
c	1 enregistrement par triangle  comprenant 
c	-type de primitive, nombre d'attributs, liste des attributs
c	-nombre de sommets et coordonnees des trois sommets du triangle (reels)

c	Fichier de parametres (lecture format libre)
c	sur lequel doit etre redirige l'entree standard 
c	ligne 1 : nje nji nja njz njs
c		nje :	nombre d'especes
c		nji :	nombre de classes de zenith
c		nja :	nombre de classes d'azimuth
c		njz :	nombre de tranches d'altitude
c		njs :	nombre d'attributs dont on veut calculer la repartition spatiale
c			(en plus du calcul de la repartition de la densite volumique)
c	ligne 2 :  dz(njz) 
c		dz(njz):epaisseur des tranches d'altitude, numerotees du haut
c			vers la bas (classe 1= la plus haute)
c	ligne 3 : lisel(njs)
c		lisel(njs): liste des numeros d'attributs dont on veut calculer 
c			    la repartition spatiale (le label est le numero 1)
c	ligne 4 :  xl,njx,dx,yl,njy,dy
c		xl et yl: dimensions totales de la maquette en x et y  
c		njx njy	: nombre de divisions dans un motif en x et en y
c		dx, dy	: dimensions d'une cellule selon x et y 


c	Fichier sortie:
c	***************

c	Les resultats sont edites dans un fichier fort.60 comprenant

c		Dimensions de la maquette 
c		Nombre de fois ou le motif est repete dans la maquette.

c	STATISTIQUES GLOBALES DE CHAQUE ESPECE:
c	Statistiques calculees sur toute la maquette:
c		Surface foliaire totale et indice foliaire 
c		Distribution des orientations des normales aux feuilles
c		(frequence en zenith et en azimuth )

c	STATISTIQUES PAR CELLULES:
c	3 lignes par cellule et par espece:
c		l1 -> je,jx,jy,jz,xvol(,,,,i)
c		l2 -> f(ji) 
c		l3 -> f(ja)
c		je	 : numero d'espece
c		jx,jy,jz : numero de cellule en x, y, et z
c		xvol(,,,,1)  : densite volumique de surface foliaire.
c		xvol(,,,,i+1): densite volumique du ieme attribut selectionne par lisel(i)
c		f(ji)	 : distribution d'angle zenital dans la cellule
c		f(ja)	 : distribution d'azimuth dans la cellule



c	Parametres 
c	**********

c	levelmax: nombre max de niveaux de subdivision d'un triangle
c	njemax	: nombre max d'especes
c	njxmax	: nombre max de subdivisions spatiales selon x
c	njymax	: nombre max de subdivisions spatiales selon y
c	njzmax	: nombre max de subdivisions spatiales selon z
c	njimax	: nombre max de classes de zenith
c	njamax 	: nombre max de classes d'azimuth
c	nattmax : nombre max d'attributs (le label compte pour 1)


c	Compilation et Execution:
c	*************************

c	f77 -extend_source -O s5.f lecpol.o -o s5
c	s5 <s5.par 	(s5.par etant le fichier de parametres).



c	Definition des noms de variables dans le code
c	*********************************************

c	nje : nombres d'especes
c	njs : nombre d'attributs selectionnes
c	nji,nja :nombre de classes de zenith et d'azimuth
c	njx,njy,njz: nombre de tranches en x, y, z
c	dx,dy,dz(jz): �paisseur des tranches en x, y, z
c	p1, p2, p3 : sommets du triangle 	




	parameter (njemax=20,njxmax=20,njymax=20,njzmax=27)
        parameter (njimax=9,njamax=12,nattmax=3)
	parameter (levelmax=6, nsmax=4**levelmax)
	parameter (nx=njemax*njxmax*njymax*njzmax)
        parameter (nxladi=nx*njimax,nxladia=nxladi*njamax)
	parameter (nxlada= nx*njamax,nxvol=nx*nattmax)
        parameter (ndisti=njemax*njimax,ndista=njemax*njamax)
	dimension p1(3),p2(3),p3(3),pg(3),jp1(3),jp2(3),jp3(3),jpg(3)
	dimension p1s(3,nsmax),p2s(3,nsmax),p3s(3,nsmax)
        dimension p13(3),p12(3),p23(3)
	dimension xvol(njemax,njxmax,njymax,njzmax,nattmax)
	dimension xladi(njemax,njxmax,njymax,njzmax,njimax)
	dimension xlada(njemax,njxmax,njymax,njzmax,njamax)
	dimension xladia(njemax,njxmax,njymax,njzmax,njimax,njamax)
	dimension dz(njzmax),bz(njzmax),volume(njzmax)
	dimension dista(njemax,njamax),disti(njemax,njimax)
	dimension xlai(njemax),levels(nsmax)
	dimension lisel(nattmax)
	double precision x_att(nattmax)
	double precision surft(njemax)
	integer*8 litp
	character*1 ntype
	character*7 file_triangles

	common/dxdydz/dx,dy,dz,njz,bz
	common/tab/xladia,xvol
	common/sommet/p1,p2,p3
	common/convert/proscal_to_surf
	common/attribut/x_att
	common/depassZ/ifirst

	external lecpol   ! $pragma C (lecpol)
	external ouvre    ! $pragma C (lecpol)
	external ferme    ! $pragma C (lecpol)
   
c	initialisation
c	**************

	file_triangles = 'fort.51'

	data volume,xlai,surft,levels/njzmax*0.0,
     &       njemax*0.0,njemax*0.d0,nsmax*0.0/
	data xvol,xladi,xlada,xladia/nxvol*0.0,
     &       nxladi*0.0,nxlada*0.0,nxladia*0.0/
	data disti,dista/ndisti*0.0,ndista*0.0/
	data ifirst/0/
	open(unit=70,file='leafarea')
	open(unit=50,file='s5.par')
	read(50,*) nje,nji,nja,njz,njs
	read(50,*) (dz(j),j=1,njz)
	if(njs.ge.1) read(50,*) (lisel(j),j=1,njs)
	read(50,*) xl,njx,dx,yl,njy,dy
c	write(6,*) nje,nji,nja,njz,njs
c	write(6,*) (dz(j),j=1,njz)
c	write(6,*) (lisel(j), j=1,njs)
c	write(6,*) xl,njx,dx,yl,njy,dy

	if((nji.gt.njimax).or.(nja.gt.njamax).
     &  or.(njx.gt.njxmax).or.(njy.gt.njymax).
     &	or.(njz.gt.njzmax).or.(nje.gt.njemax))
     &    stop 'ERREUR DE DIMENSIONNEMENT DES TABLEAUX'
 
	di = 90.0/nji
	da = 360.0/nja
	dxdy = dx*dy
	xymaille = (xl*yl)/(njx*njy*dx*dy)
	ns = 0

c	calcul des bornes superieures des classes d'altitudes
c       *****************************************************
	bz(njz) = dz(njz)
	do 20 j=  njz-1, 1, -1
	bz(j) = dz(j) + bz(j+1)
  20	continue


c	lecture d'un triangle, calcul de sa classe d'orientation et de la cellule de chaque sommet
c	******************************************************************************************

	call ouvre(file_triangles)
   1	level = 0
	call lecpol(ltest,ntype,natt,x_att,nsom,p1,p2,p3)
	if(ltest.eq.-1) goto 999
        je=(x_att(1)/10**6)
	je=int(je/10**5)

        if((je.le.0).or.(je.gt.nje).or.(natt.le.0)
     &     .or.(natt.gt.nattmax).or.(nsom.ne.3))
     &     stop 'FORMAT INCORRECT OU DONNEES ERRONEES SUR FORT.51'

	proscal_to_surf=0.5
	litp = (x_att(1)/10**6)
	jfeuille = int(x_att(1)-litp*10**6)/10**3
        if(jfeuille.eq.0) proscal_to_surf=0.25

	call normal(p1,p2,p3,xi,a) 
	ji = 1+int(xi/di)
	ja = 1+int(a/da)
	ji = min(nji,max(1,ji))
	ja = min(nja,max(1,ja))
	call calcjp(p1,p2,p3,jp1,jp2,jp3,itest)

   2    continue


c	Si les trois sommets du triangle appartiennent a une meme cellule, le triangle est affecte a xladia.
c	On prend en compte le triangle suivant, soit sur le disque, soit dans la liste des sous-triangles. 
c	**************************************************************************************************
	if (itest.eq.0) then
	  call affect(je,jp1,ji,ja,njx,njy,njs,lisel)
	  if(ns.eq.0) go to 1
	  level = levels(ns)
	  do 11 i= 1,3
	    p1(i) = p1s(i,ns)
	    p2(i) = p2s(i,ns)
	    p3(i) = p3s(i,ns)
  11	  continue
	  call calcjp(p1,p2,p3,jp1,jp2,jp3,itest)
	  ns = ns - 1
	  go to 2
	endif


c	le triangle est a cheval sur plusieurs cellules, mais si le niveau de subdivision est superieur a
c	levelmax, on l'affecte a la cellule de son centre de gravite.
c	On prend en compte le triangle suivant, soit sur le disque, soit dans la liste des sous-triangles. 
c	**************************************************************************************************

	if(level.ge.levelmax) then
	  do 30 i=1,3
	    pg(i) = (p1(i)+p2(i)+p3(i))/3.
  30	  continue
	  jpg(1) = nint(pg(1)/dx)
	  jpg(2) = nint(pg(2)/dy)
	  call class(pg(3),jpg(3))
	  call affect(je,jpg,ji,ja,njx,njy,njs,lisel)
	  if(ns.eq.0) go to 1
	  level = levels(ns)
	  do 12 i= 1,3
	    p1(i) = p1s(i,ns)
	    p2(i) = p2s(i,ns)
	    p3(i) = p3s(i,ns)
  12	  continue
	  call calcjp(p1,p2,p3,jp1,jp2,jp3,itest)
	  ns = ns - 1
	  go to 2
	endif


c	le triangle courant n'a pas ete affecte. On le subdivise en quatre sous-triangles. 
c	Trois sous-triangles sont stockes dans les tableaux pis. Le dernier est analyse
c	*********************************************************************************

	do 10 i= 1,3
	  p12(i) = 0.5*(p1(i)+p2(i))
	  p13(i) = 0.5*(p1(i)+p3(i))
	  p23(i) = 0.5*(p2(i)+p3(i))
  10	continue

	ns1 = ns+1
	ns2 = ns+2
	ns3 = ns+3
	ns = ns3

	level = level + 1
	levels(ns1) = level
	levels(ns2) = level
	levels(ns3) = level

	do 40 i = 1,3
	  p1s(i,ns1) = p2(i)
	  p2s(i,ns1) = p12(i)
	  p3s(i,ns1) = p23(i)
	  p1s(i,ns2) = p3(i)
	  p2s(i,ns2) = p23(i)
	  p3s(i,ns2) = p13(i)
	  p1s(i,ns3) = p12(i)
	  p2s(i,ns3) = p13(i)
	  p3s(i,ns3) = p23(i)
  40	continue

	  do 13 i= 1,3
	    p2(i) = p12(i)
	    p3(i) = p13(i)
  13	  continue
	call calcjp(p1,p2,p3,jp1,jp2,jp3,itest)
	go to 2

999	continue
	
	call ferme()

c	Calcul des distribution spatiales de surface foliaire et des frequences d'angle
c	******************************************************************************

	do 52 jz = 1,njz 
   	volume(jz) = xymaille*dxdy*dz(jz)
  52	continue

	do 53 je = 1,nje
	do 53 jx = 1,njx 
	do 53 jy = 1,njy
	do 53 jz = 1,njz

	  do 54 ji = 1, nji
	  do 54 ja = 1,nja
	      tp = xladia(je,jx,jy,jz,ji,ja)
	      xladi(je,jx,jy,jz,ji) = xladi(je,jx,jy,jz,ji) + tp
	      xlada(je,jx,jy,jz,ja) = xlada(je,jx,jy,jz,ja) + tp
	      disti(je,ji) = disti(je,ji) + tp
  	      dista(je,ja) = dista(je,ja) + tp
  54	  continue

c	  passage aux frequences d'angles....
	  tp2 = xvol(je,jx,jy,jz,1)
	  if (tp2.gt.0.)  then
	  do 56 ji = 1,nji
  56	  xladi(je,jx,jy,jz,ji) = xladi(je,jx,jy,jz,ji)/tp2
	  do 57 ja = 1,nja
  57	  xlada(je,jx,jy,jz,ja) = xlada(je,jx,jy,jz,ja)/tp2
	  endif

	  surft(je) = surft(je) + tp2
  53	continue

c	Calcul des distributions globales d'orientation des feuilles et de l'indice foliaire
	do 61 je = 1,nje
	tp3=real(surft(je))
  	xlai(je) = tp3/(xl*yl)
	do 62 ja = 1,nja
  62	dista(je,ja) = dista(je,ja)/tp3
	do 63 ji = 1,nji
  63	disti(je,ji) = disti(je,ji)/tp3
  61	continue



c	Editions
c	********
	write (60,507) xl,yl
	write (60,508) xymaille
	write (60,505)
	do 82 je = 1,nje
	write(60,502) je, surft(je), xlai(je)
	write(60,503) (disti(je,ji),ji = 1,nji)
	write(60,504) (dista(je,ja),ja = 1,nja)
  82	continue

	write (60,506)
	do 81 je = 1,nje
	do 81 jx = 1,njx 
	do 81 jy = 1,njy
	do 81 jz = 1,njz
c	write(6,*) xvol(je,jx,jy,jz,1), xladi(je,jx,jy,jz,1)	
	if(xvol(je,jx,jy,jz,1).gt.0.) then 
	write(60,500) jx,jy,jz,je,
     &       ((xvol(je,jx,jy,jz,js)/volume(jz)),js=1,njs+1)
	write(60,501) (xladi(je,jx,jy,jz,ji),ji =1,nji)
	write(60,501) (xlada(je,jx,jy,jz,ja),ja =1,nja)
	write(70,509) jx,jy,jz,je,
     &       (xvol(je,jx,jy,jz,1)/volume(jz)),
     &       (xladi(je,jx,jy,jz,ji),ji =1,nji),
     &       (xlada(je,jx,jy,jz,ja),ja =1,nja)
	endif
  81	continue

	close(60)
	close(70)

	stop
 500	format (4(1x,i2),10(2x,f9.3)) 
 501	format (99(1x,f6.4)) 
 502	format (1x,'espece: ',i2,4x, 'surface foliaire: ',d8.3,
     &          4x, 'lai : ',f7.4)
 503	format (1x,'distribution en zenith:  ',99(f6.4,1x))
 504	format (1x,'distribution en azimuth: ',99(f6.4,1x))
 505	format (1x,/,1x,'STATISTIQUES GLOBALES DE CHAQUE ESPECE',/)
 506	format (1x,//,1x,'STATISTIQUES PAR CELLULE',/)
 507    format (1x,'dimensions de la maquette (x,y):',2x, 2(f7.3,2x))
 508	format (1x,'nombre de repetitions du motif:',2x,f5.1)
 509    format (4(1x,i2),(1x,f10.3),30(1x,f6.3))
 510    format (4(1x,i2),10(1x,f9.3))
 511    format (4(1x,i2),10(1x,f9.3))

	end
c***********************************************************************************

c***********************************************************************************
	subroutine calcjp(p1,p2,p3,jp1,jp2,jp3,itest)

c	Calcul des positions en unites dx dy dz 
c	nb: le mode d'affectation est tel que la position (0,0,1)
c       est centree sur (0, 0, 0.5*bz(njz))

	parameter (njzmax=27) 
	dimension p1(3),p2(3),p3(3),jp1(3),jp2(3),jp3(3),
     &            dz(njzmax),bz(njzmax)
	common/dxdydz/dx,dy,dz,njz,bz

	jp1(1) = nint(p1(1)/dx)
	jp2(1) = nint(p2(1)/dx)
	jp3(1) = nint(p3(1)/dx)   

	jp1(2) = nint(p1(2)/dy)
	jp2(2) = nint(p2(2)/dy)
	jp3(2) = nint(p3(2)/dy)

	call class(p1(3),jp1(3))
	call class(p2(3),jp2(3))
	call class(p3(3),jp3(3))

	itest12 = (jp1(1)-jp2(1))**2 + (jp1(2)-jp2(2))**2
     &            + (jp1(3)-jp2(3))**2 
	itest13 = (jp1(1)-jp3(1))**2 + (jp1(2)-jp3(2))**2
     &            + (jp1(3)-jp3(3))**2 
	itest = itest12+itest13

	return
	end
c***********************************************************************************

c***********************************************************************************
	subroutine class(z,jz)

c	Classement de la valeur z. Les classes sont definies par leur borne
c	sup�rieures, le numero de classe est dans le sens des z decroissants.
c	nb.: les z negatifs sont classes dans la classe la plus basse...

	parameter (njzmax=27)
	dimension bsup(njzmax),dz(njzmax)
	common/dxdydz/dx,dy,dz,njz,bsup
	common/depassZ/ifirst

	if ((z.ge.bsup(1)).and.(ifirst.ne.1)) then
	  print *,"subroutine class : depassement en z",z
	  ifirst=1
	endif

	do 1 j= 2,njz
	if(z.gt.bsup(j)) go to 2
   1	continue
	jz = njz 
	return
   2	jz = j-1 
	return
	end

c***********************************************************************************

c***********************************************************************************
	function proscal(a,b)
	dimension a(3),b(3)
	proscal=0.0
	do 1 i=1,3
1	proscal=proscal+a(i)*b(i)
	return
	end
c***********************************************************************************

c***********************************************************************************
	subroutine provec(a,b,c)
	dimension a(3),b(3),c(3)
	c(1)=a(2)*b(3)-a(3)*b(2) 
	c(2)=-a(1)*b(3)+a(3)*b(1) 
	c(3)=a(1)*b(2)-a(2)*b(1)
	return
	end
c***********************************************************************************

c***********************************************************************************
	subroutine norme(x,xn)
	dimension x(3),xn(3)
	xnor=sqrt(proscal(x,x))
	if (xnor.lt.1e-10) goto 99
	do 1 i=1,3
1	xn(i)=x(i)/xnor
	return
99	stop 9
	end
c***********************************************************************************

c***********************************************************************************
	subroutine affect(je,jp,ji,ja,njx,njy,njs,lisel)

c	mise a jour des tableaux xladia et xvol 
c	avec un triangle dont les 3 sommets appartiennent a une meme cellule

c	parameter (njemax=2,njxmax=11,njymax=11,njzmax=27)
c        parameter (njimax=9,njamax=12,nattmax=3)
	parameter (njemax=20,njxmax=20,njymax=20,njzmax=27)
        parameter (njimax=9,njamax=12,nattmax=3)
	dimension a(3),b(3),c(3),p1(3),p2(3),p3(3),jp(3)
	dimension xladia(njemax,njxmax,njymax,njzmax,njimax,njamax)
	dimension xvol(njemax,njxmax,njymax,njzmax,nattmax)
	double precision x_att(nattmax)
	dimension lisel(nattmax)
	common/tab/xladia,xvol
	common/sommet/p1,p2,p3
	common/convert/proscal_to_surf
	common/attribut/x_att
	
	jx = mod(jp(1),njx)
	jy = mod(jp(2),njy)
	jz = jp(3)

c	decalage de 1 et prise en compte de possibles coordonnees negatives
	jx = 1+mod((jx+njx),njx)
	jy = 1+mod((jy+njy),njy)

	do 1 i=1,3
          a(i)=p2(i)-p1(i)
          b(i)=p3(i)-p1(i)
  1	continue

	call provec(a,b,c)
	tp=sqrt(proscal(c,c))
	
	surftri_e = proscal_to_surf*tp
	xladia(je,jx,jy,jz,ji,ja) = xladia(je,jx,jy,jz,ji,ja) + surftri_e
	xvol(je,jx,jy,jz,1)       = xvol(je,jx,jy,jz,1)       + surftri_e

	surftri_r = 0.5*tp
	if(njs.le.0) return
	do 2 js=1,njs
	  js1=js+1
	  xvol(je,jx,jy,jz,js1)=xvol(je,jx,jy,jz,js1)
     &        + surftri_r*x_att(lisel(js))
	 
  2	continue

	return
	end
c***********************************************************************************

c***********************************************************************************

	subroutine normal (p1,p2,p3,zen,az)

	dimension p1(3),p2(3),p3(3),p3p2(3)

	do 1 i=1,3
	p3p2(i) = p2(i)-p3(i) 
   1	continue

	u =      p1(2)*p3p2(3) - p1(3)*p3p2(2) + p2(2)*p3(3)-p3(2)*p2(3)
	v = -1.*(p1(1)*p3p2(3) - p1(3)*p3p2(1) + p2(1)*p3(3)-p3(1)*p2(3))
	w =      p1(1)*p3p2(2) - p1(2)*p3p2(1) + p2(1)*p3(2)-p3(1)*p2(2)
	r2 = u**2+v**2
	if(r2.eq.0.) go to 2
	r1 = sqrt(r2)
	zen = 57.29577951*acos(abs(w/sqrt(r2+w**2)))
	az  = 57.29577951*atan2(u/r1,v/r1)
	if(az.lt.0.) az = 360.0 + az
	return
   2    zen = 0.
	az  = 0.
	return
	end

c***********************************************************************************

c***********************************************************************************
