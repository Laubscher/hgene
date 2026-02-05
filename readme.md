# Marche à suivre pour convertir Docker en Singularity

## Étapes de création de l'image

Dans le dossier du projet :

1. **Construire l'image Docker** :

   ```
   sudo docker build -t hgene .
   ```

2. **Sauvegarder l'image Docker dans un fichier `.tar`** :

   ```
   sudo docker save hgene -o hgene.tar
   ```

3. **Convertir l'image `.tar` en `.sif` avec Singularity** :

   ```
   sudo singularity build hgene.sif docker-archive://hgene.tar
   ```

## Usage du conteneur Singularity

Exécutez le conteneur Singularity avec les options suivantes :

```
singularity exec hgene.sif Hgene -v [VIRUS] [fastqPrexif] [CPU]
```

### Options disponibles

- **-h** : Affiche cette aide.
- **-v** : Spécifie le virus (par exemple, HHV1, HHV2, CMV).
- **[fastqPrefix]** : Préfixe des fichiers FASTQ.
- **[nb de CPU]** : Nombre de CPU à utiliser.
