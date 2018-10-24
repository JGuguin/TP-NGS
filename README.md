## Reproducing analysis (at least trying to) of Sohail _et al_, Science (2017)


Sohail, Mashaal, _et al_. Negative selection in humans and 
fruit flies involves synergistic epistasis.
 Science 356.6337 (2017): 539-542.
 
 https://doi.org/10.1126/science.aah523874
 
 ### Donwload the scripts
```
$ git clone https://github.com/JGuguin/TP-NGS
$ cd TP-NGS
```

 ### Install VcfTools, Bedtools and the custom scripts
 Some steps requires you to be superuser (sudo), but can be skipped 
 if you already have the dependencies installed.
```
$ chmod a+x src/installation.sh
$ src/installation.sh
```
### Run the analysis
```
$ chmod a+x src/analysis.sh
$ src/analysis.sh
```
