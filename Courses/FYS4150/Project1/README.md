# FYS4150 - Project 1

NB! All koden må bli kjørt i en mappe med en "Results" mappe for at de kan lagre output filer.

Oversikt over hva filene her er og hvordan de brukes:

- Project1Notebook.ipynb er en jupyter notebook som leser inn resultater fra "Results" og lager ulike plots som lagres i "Results"
- lib.cpp er kode tatt fra https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/src/learningcpp som løser ligningssett med LU faktorisering
- tridiag.cpp er den tilpassede Thomas algoritmen. Den tar en input i kommandolinjen som er størrelsesorden til n. f.eks: "tridiag.exe 5"
- tridiagSlow.cpp er den generelle Thomas algoritmen. Den tar en input i kommandolinjen som er størrelsesorden til n. f.eks: "tridiagSlow.exe 5"
- tridiagLU.cpp bruker LU faktorisering til å løse ligningen. Den tar en input i kommandolinjen som er størrelsesorden til n. f.eks: "tridiagLU.exe 3"
- tridiagTiming.cpp kjører de to første metoden opp til n = 10^7 og LU metoden opp til n = 10^5. Skriver resultatet til timing.txt. timingFinal.txt holder nå et fullt resultat som ikke blir overskrivet.
- Results mappa holder resultat fra utregninger og plotting.
- Latex mappa holder latex filer som ble brukt til å lage rapporten til dette prosjektet, samt rapporten.

Output filene får tilnamn: sl{k} (slow), td{k} (tridiagonal), lu{k} (LU), hvor k er størrelsesorden til n, eller input parameteren til programmet.

Resultat fila td7.txt ble brukt til å produsere et plot i Project1Notebook.ipynb, men finnes ikke lengre i "Results" mappa, fordi den er 0.6GB stor. Det tar lang tid å gjenskape den, men det kan greit gjøres med kommandoen "tridiag.exe 7".
