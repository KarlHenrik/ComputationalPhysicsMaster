# Prosjekt 2

I dette prosjektet løser vi tre ulike diffligninger ved å diskretisere de og sette de opp som egenverdiproblemer. Den første ligningen finner formen til en bøyende bjelke. Den andre finner den radielle komponenten til bølgefunksjonen til et elektron i et harmonisk oscillator potensial. Den tredje finner avstanden mellom to interagerende elektroner i et harmonisk oscillator potensial.

- `buckling.cpp` løser ligningen for den bøyende bjelken. Resultatene blir lagret i filene som begynner på `0`.
- `bucklingTiming.cpp` tar tiden metoden bruker på å løse ligningen for den bøyende bjelken. Resultatene lagres i `bucklingTiming.txt`
- `jacobiSolve.cpp` har funksjoner som implementerer Jacobis egenverdialgoritme, her er koden for det meste tatt fra https://compphysics.github.io/ComputationalPhysics/doc/pub/eigvalues/html/eigvalues.html , med noen endringer
- `makefile` inneholder instrukser for å kompilere alle c++ programmene
- `oneElectron.cpp` løser ligningen for den ett elektron. Resultatene blir lagret i filene som begynner på `1E`.
- `oneElectronRhoN.cpp` finner feilen i den første egenverdien for løsningen av ett elektron ligningen, for ulike kombinasjoner av N og rhoN. Resultat lagres i `oneElectronRhoN.txt`
- `Project2Notebook.ipynb` inneholder all bearbeiding og plotting av resultater
- `testing.cpp` inneholder enhetstester for `jacobiSolve.cpp` og tilknyttede funksjoner. Tester også løsningen av bestemte problemer.
- `twoElecton.cpp` løser ligningen for den to elektron. Resultatene blir lagret i filene som begynner på `2E`. Her er det mange outputfiler siden mange parametre testes og dette var enklest for et så lite prosjekt.

## Results

I dette prosjektet var det mange beregninger som ble gjort, men bare noen av resultatene blir analysert. Filene som begynner med `0` er resultater for den bøyende bjelken. Filene som begynner med `1E` er resultater for ett elektron, og `2E` for to elektron. Videre skilles filene mellom hvilken metode som regnet ut resultatet, som kan være en analytisk formel, Jacobis egenverdimetode eller funksjoner fra C++ biblioteket Armadillo. Endingene på filene sier om de holder på de sorterte egenverdiere eller egenvektorene tilhørende egenverdiene.

Noen filer skiller seg ut. Dette er `bucklingTiming` som holder tidtakingen for de ulike metodene anvendt på den bøyende bjelken. `2EGS...` filene som holder egenvektoren tilhørende den laveste egenverdien for to elektron problemet regnet ut med ulik $\omega_r$ (en $\omega_r$ per fil). Og til slutt `oneElectronRhoN` som holder absolutt feil for den første egenverdien for ett elektron for ulike kombinasjoner av $N$ og $\rho_N$.

## Figures

Alle plots kommer fra notebooken `Project2Notebook.ipynb`.
