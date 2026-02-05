#!/usr/bin/env python3

import sys

vcf=open(str(sys.argv[1]), "r")

for l in vcf:
  if l[0]!="#":
    ligne=l.rstrip().split("\t")
    if len(ligne[3]) > len(ligne[4]):  # si la longueur de la REF est plus grande que l'ALT il s'agit d'une deletion 
      if float(ligne[-1].split(";")[1].split("=")[1]) >= 0.4:
        print(l.rstrip())
      else:
        pass #on peut compter les entry supr pour les log
    else:
      print(l.rstrip())   # si pas deletion print
  else: 
    print(l.rstrip())     # si commentaire /header print 
 
vcf.close()

