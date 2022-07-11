# Gestion de la lumière entre des modèles de plantes

## Notes sur le fonctionnement
Le package PyRATP est appelé en local et fait parti du package de lightvegemanager. Son développement est en parallèle de l'outil de gestion.

## Utilisation sur le meso@LR avec un conteneur singularity comme environnement d'exécution
* avoir un pyratp.so compilé par le conteneur pour être utilisable
* les scripts python à exécuter doivent avoir la permission d'exécution : `chmod +x script.py`
* tous les fichiers compatibles pour unix 

### compilation de la partie fortran de PyRATP
Il y a un script pour linux 64 bits et Windows 64 bits pour compiler la partie fortran avec f2py.
Ces scripts sont à lancer dans un environnement avec numpy d'installé (normalement installé dans le conteneur), avant de run les simulations.

1) d'abord rentrer dans le conteneur
    ```bash
    singularity shell mycontainer.sif
    ```

2) puis lancer la compilation fortran -> python 
    ```bash
    cd PyRATP/f90
    bash compile_pyratp_lin64.sh
    ```

3) vérifier que la nouvelle librairie a bien été installée dans PyRATP/pyratp sous le nom pyratp.so
    ```bash
    ls -lt PyRATP/pyratp/*.so
    ```

### mise à jour du dépôt 
1) revenir au dernier commit
    ```bash
    git reset --hard HEAD
    ```
2) lancer dernière ligne d'un fichier
    ```bash
    # lance dos2unix pour convertir le format des fichiers du depot
    tail -1 myinstall_lightvegemanager.sh | bash
    ```