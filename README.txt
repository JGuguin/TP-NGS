##################################################
##TRAVAUX PRATIQUES - NEXT GENERATION SEQUENCING##
##################################################

Encadrant: Thibault Latrille
Dates :
	Première session : 17.09.2018 - 19.09.2018
	Deuxième session : 22.10.2018 - 24.10.2018
	
Objectifs : Etduier le polymorphisme, notamment celui au niveau des mutations délétères, au sein des populations humaines avec reproduction de l'analyse de Sohail et al.

Article :
Sohail, Mashaal, et al. Negative selection in humans and fruit flies involves synergistic epistasis. Science 356.6337 (2017): 539-542.
https://doi.org/10.1126/science.aah523874
 -> Les auteurs montrent dans cet article que les mutations délétères au sein des populations humaines et de drosophiles, sont majoritairement soumises à un phénomène d'épistasie synergétique. Ce phénomène est associé à de la répulsion entre les mutations : deux mutations délétères ont un effet combinatoire et s'augmentent l'une l'autre leur effet délétère (épistasie synergétique). C'est pourquoi, on les trouve rarement associées au sein d'un même individu (répulsion)

--------------------
 
Déroulement de l'analyse :
	1°/ Installation des modules nécessaires à l'analyse
	2°/ Alignement du séquençage de l'exome pour tous les individus sur un génome de référence
	3°/ Appel de variants afin de savoir où se situent les polymorphismes dans le génome des individus (fichier .vcf)
	4°/ Combinaison du fichier .vcf contenant les polymorphismes avec un fichier .gtf de référence et un fichier .fasta contenant le séquence des transcrits afin de connaître les mutations associées à chaque SNP
	5°/ Classification des mutations en trois catégories : synonymes, non-synonymes et perte de fonction (stop)
	6°/ Tri des individus par population
	7°/ Pour chaque population, calcul du ratio sigma²/Va pour les mutations synonymes, non synonymes et stop afin de voir l'effet des mutations délétères entre elles
	
Les étapes 1 à 3 ont été réalisées lors de la première session. Les étapes 4 à 7 ont été réalisées lors de la deuxième session.
L'ensemble de l'analyse a été réalisée sur le chromosome 20 uniquement pour des raisons de temps d'analyse et de quantité de données.
Le jeu de données utilisé (génomes et individus) provient du projet 1000 Genome Project qui regroupe 2504 individus provenant d'une trentaine de populations différentes.

Tous les scripts se trouvent sur le github : https://github.com/JGuguin/TP-NGS

--------------------
--------------------

####################
##Première session##
####################

Objectif : Obtenir un fichier .vcf contenant l'emplacement des SNPs et le génotype à ces SNPs (variant calling) pour une famille de trois individus (trio mère-père-fille)

--------------------

1°/ Installation des modules nécessaires à l'analyse

Scripts nécessaires à la première session -> Voir le script installation.sh

$ chmod a+x installation.sh
$ installation.sh

Scripts nécessaires à la deuxième session -> Voir le script installation.sh dans le dossier src

$ chmod a+x src/installation.sh
$ src/installation.sh

---------------------

2°/ Alignement du séquençage de l'exome pour tous les individus sur un génome de référence

Voir le script mapping.sh

$ chmod a+x mapping.sh
$ mapping.sh

---------------------

3°/ Appel de variants afin de savoir où se situent les polymorphismes dans le génome des individus (fichier .vcf)

Voir le script variant-calling.sh

$ chmod a+x variant-calling.sh
$ variant-calling.sh

On peut ensuite analyser les variants chez les individus afin de voir où se trouvent les différences entre les individus et entre les individus et la référence.
Comme les individus sont apparentés, il faut aussi l'indiquer au programme car cela peut changer l'analyse.
Tout cela se trouve dans le script trio-analysis.sh

$ chmod a+x trio-analysis.sh
$ trio-analysis.sh

---------------------

L'analyse de la première session nous a ainsi permis de savoir comment obtenir les fichiers .vcf qui seront utilisés lors de la deuxième session, même si ces fichiers seront plus denses car contiendront beaucoup plus d'individus.
De nombreux variants (SNP différents de la référence) sont présents dans la famille 


---------------------
---------------------


####################
##Deuxième session##
####################

Objectif : Faire une comparaison des polymorphismes dans chaque population afin de voir l'interaction entre les mutations délétères

Trois types d'interaction sont possibles :
	- Epistasie antagonistique (phénomène d'attraction) : les mutations ont un effet bénéfique entre elles et diminuent ensemble le poids d'une mutation individuelle, d'où une attraction -> Avoir n mutations sera plus bénéfique qu'en avoir une seule
	- Indépendance : l'effet des mutations ne s'additionne pas -> avoir une ou n mutations ne change rien
	- Epistasie synergétique (phénomène de répulsion) : les mutations ont un effet délétère accru entre elles et augmentent ensemble le poids d'une mutation individuelle, d'où une répulsion -> Avoir n mutations sera moins bénéfique, pire qu'en avoir une seule
	
Dans le papier de Sohail et al., les auteurs se sont essentiellement intéressés aux SNPs présents une seule fois dans le génome (singleton) ou ayant une faible fréquence. Ils ont trouvé que les mutations délétères étaient ainsi soumises à un phénomène d'épistasie synergétique pour les populations qu'ils ont étudié.
Or, plusieurs questions peuvent être soulevées par leur étude, la première étant que l'analyse de singleton n'est peut-être pas particulièrement pertinente, dans le sens où cela est très sélectif et met un poids assez important sur les mutations ponctuelles de novo (survenant lors de la méiose) qui ne sont pourtant peut-être pas délétères à l'état homozygote. C'est le cas par exemple de KO pour des gènes qui n'affectent pas de fonctions vitales immédiates (comme par exemple, certains individus sont KO pour une enzyme métabolisantl'aspirine. C'est une perte de fonction, mais qui n'est pas délétère)
Egalement, leur étude est restreinte à trois populations humaines seulement. Ainsi, on peut se demander si ce qu'ils observent est réellement généralisable à l'ensemble des populations humaines ou non, et si ce n'est pas le cas, est-ce que des phénomènes de dynamique de populations peuvent être la cause de cette épistasie synergétique.

Dans notre étude, on s'intéresse à une trentaine de populations humaines, ainsi qu'à différent cut-off (nombre de fois qu'un SNP particulier est présent dans la population) afin d'être moins restrictif qu'un simple cut-off de un (correspondant aux singletons).

---------------------

L'ensemble de notre script est présent dans le fichier analysis.sh du dossier src (étapes 4 à 7)

$ chmod a+x src/analysis.sh
$ src/analysis.sh

Le script vcf_analysis.py du dossier src n'est pas utilisé ici car est inclus dans le script vcf_meta_analysis.py

---------------------

Afin de savoir dans quel modèle d'interaction on se trouve dans chaque population, on réalise un calcul de sigma²/Va. Ce ratio de la somme des variances sur la variance de la somme peut ensuite être comparé à 1 afin de savoir le modèle d'interaction.
Afin de s'affranchir des dynamiques de populations influençant l'apparition des mutations, on considère ici que la référence 1 est en fait le ratio des variance obtenu pour les mutations synonymes. En effet, ces mutations ne sont pas soumises à la sélection et permettent donc de s'affranchir des effets populationnels.

L'analyse a révélé que le choix du cut-off est extrêmement important pour l'analyse, car on ne sélectionne pas les mêmes mutations. Ainsi des mutations délétères ayant un effet répulsif entre elles à un cut-off de un peuvent se retrouver noyées dans des mutations indépendantes à un cut-off de 5 et l'effet de répulsion n'est alors plus observables.
De même, des effets d'attraction (épistasie antagonistique) apparaissent à des cut-off élevés, probablement car la probabilité de sélectionner des mutations faisant parties d'un même haplotype augmente, et comme on sélectionne ainsi des haplotypes, l'effet observé est forcément attractif car les mutations ségrègent ensemble.
Mais cela ne s'observe que pour certaines populations, probablement car dans ces populations, ces mutations ont un rôle fondamental. La majorité des populations concernées par cela sont les populations asiatiques. Ainsi, les pertes de fonction dans ces populations ont un effet attractif entre elles, ce qui est peut-être lié à des fonctions métaboliques particulières lié à un régime alimentaire spécifique souvent constitué de riz, de poisson et d'algues (sushis !) 
Egalement, il faut souligner le fait que certaines mutations délétères à l'état homozygote sont peut-être bénéfique à l'état hétérozygote et vont faire de l'épistasie antagonistique entre elles. On peut par exemple cité l'anémie falciforme en Afrique qui, bien que délétère à l'état homozygote, onfère un fort avantage sélectif à l'état hétérozygote.

De plus, prendre un cut-off de 2504 n'était pas vraiment pertinent, dans le sens où cela sélectionne absolument tous les SNP présents dans la population car l'effectif de celle-ci est très inférieur à 2504 (le maximum étant 661 individus dans la population -> a un cut-off de 1322, on sélectionnerait déjà absolument tous les SNPs dans l'hypothèse absolue où certains seraient présents chez tous les individus et de manière homozygote. Cela serait possible si la référence diffère grandement de la population étudiée).
Concernant la référence justement, il est possible de certains SNP aient été identifiés dans des populations et pas d'autres du fait d'une différence majeure avec la référence (qui est européenne). Il faudrait peut-être dans ce cas adapter la référence à chaque population étudiée.

Toutefois, si on s'intéresse aux mutations vraiment rares dans les populations (fréquence inférieure à 0,01%), on remarque que celles-ci sont majoritairement indépendantes et quelques fois répulsives, alors que Sohail et al. ont démontré qu'il y avait majoritairement de la répulsion.
