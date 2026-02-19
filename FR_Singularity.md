# Marche à suivre pour convertir Docker en Singularity

## Étapes de création de l'image

Dans le dossier du projet :

1. **Construire l'image Docker** :

   `sudo docker build -t hgene:1.0.X \
  --build-arg VERSION=1.0.X \
  --build-arg VCS_REF="$(git rev-parse --short HEAD 2>/dev/null || echo unknown)" \
  --build-arg BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)" `

lofreq dois être dans bin

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
singularity exec hgene.sif hgene -v <VIRUS> -c [CPU] <fastq>
```

### Options disponibles

- **-v** : Spécifie le virus (par exemple, HHV1, HHV2).
- **-c** : Nombre de CPU à utiliser.
- **<fastq>** : Préfixe du fichier FASTQ.
