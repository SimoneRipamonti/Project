#Codice per la risoluzione del sistema di Darcy in 1D e risoluzione equazione di trasporto e reazione di un tracciante.

#file data.pot
Contiene tutti i dati necessari al problema

#muparser_fun.hpp/cpp
Definendo variabili come la porosità, la sorgente esterna o la permeabilità come funzioni dello spazio (handle function in Matlab) bisogna acquisirle tramite le funzionalità della libreria muparser, utilizzata in questa classe.
 
#parameters.hpp/cpp
Grazie a GetPot si acquisiscono i dati dal file data.pot e li si salvano in opportune variabili.
Per le variabili come la permeabilità, la porosità e la sorgente esterna si utilizza la classe muparser_fun.

#matrix.hpp/matrix.cpp
Si definisce prima di tutto l'Abstract_Matrix class generale per tutte le matrici (operatori) necessari alla risoluzione del sistema di Darcy e del problema di trasporto. In seguito si definisicono le singole matrici: 
A è la matrice di massa dell'equazione di Darcy;
B è la matrice del sistema punto sella (sempre di Darcy)(NOTA: in questo codice viene definita come la Trasposta di quella riportata nelle note, per comodità di implementazione);
C è la matrice di massa del tracciante;
F_piu è la matrice relativa ai flussi uscenti; 
F_meno è la matrice relativa ai flussi entranti;
(F_piu e F_meno saltano fuori dall'integrazione sul bordo del termine di trasporto).

#functions.hpp/cpp
Si definiscono due funzioni:
1) La funzione che assembla il sistema di Darcy, fornendo come output la matrice complessiva M e il suo rhs.
2) La funzione che permette di scrivere su file .csv (che poi verranno plottati) i risultati sia della velocità che della pressione. (In output si hanno i file velocity.csv e pressure.csv).

#velocity.sh,pressure.sh
Due file bash che lanciati dopo aver lanciato il main (cioè dopo che i file .csv sono stati riempiti coi dati) permettono di plottare il risultato tramite gnuplot. 

IMPORTANTE1: PER COMPILARE BISOGNA CAMBIARE LE VARIABILI DEL MAKEFILE INSERENDO I PERCORSI CHE CONDUCONO ALLE LIBRERIE NECESSARIE.

IMPORTANTE2: UNA VOLTA CAMBIATI I PERCORSI NEI MAKEFILE SCRIVE SU TERMINALE "module load boost" PER ATTIVARE LE FUNZIONALITÀ DI GNUPLOT
