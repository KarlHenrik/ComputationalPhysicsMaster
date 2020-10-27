# Prosjekt 3

Rapporten til prosjektet finnes i Report mappen.

## Koden

Code mappen inneholder mange kodefiler:

- earthSun.cpp er programmet som sammenligner Velocity Verlet med Euler. Her tolket jeg oppgaven som at dette ikke skulle bli gjort objektorientert, så dette er et lite fleksibelt program som produserer resultatene oppgaven spør om. Det kjøres uten inputparametre og produserer filer i en mappe Output/EulerVsVverlet hvis den finnes.
- main.cpp er programmet som produserer alle andre resultatene i prosjektet. Her settes opp og simuleres mange ulike kombinasjoner av planeter, tidssteg og tyngdekrefter. Det kjøres uten inputparametre(siden vi må prøve ut veldig mange spesifikke systemer), og produserer filer i en mappe Output hvis den finnes.
- vec3.cpp er en 3d vektorklasse fra kurssiden som støtter mange matematiske operasjoner.
- celestialBody.cpp er en klasse som lar oss lage celestialBody med en 3d vektor posisjon, hastighet, akselerasjon, navn og masse. Disse holdes orden på i et solarSystem objekt.
- solarSystem.cpp er en klasse som lar oss lage solarSystem objekter som holder på celestialBody objekter, regner ut akselerasjon til disse gitt en tyngdekraft, skriver navnene og massene til en gitt fil, eller posisjonene og hastighetene til en gitt fil. Vi kan også gjøre massesenteret til systemet til origo, eller endre hastigheten til en av planetene slik at bevegelsesmengden til systemet blir 0 (endrer på den første planeten, eller den siste hvis inputparameter er -1). Vi legger enten til en og en planet med parametrene til planeten eller ett celestialBody objekt. Vi kan også lese inn planeter fra en fil på formen til planetData.txt (dette blir gjort i main.cpp for å simulere alle planetene)
- vVerlet.cpp lar oss lage vVerlet objekter som lar oss bevege celestialBody objektene i et solarSystem objekt ett tidssteg frem ved hjelp av Velocity Verlet algoritmen.
- euler.cpp lar oss lage euler objekter som lar oss bevege celestialBody objektene i et solarSystem objekt ett tidssteg frem ved hjelp av Euler metoden.
- UnitTests/testing.cpp inneholder tester for disse klassene. Den kjøres uten inputparametre.

Oppskriftene for kompileringen av disse filene er gitt i makefilen.

## Analysen

Resultatene fra kodekjøringen finnes i Output mappen. Disse resultatene blir gjort om til figurer i de tre Jupyter notebookene i Code mappen. AnalysisTwoBody.ipynb inkluderer også en liten forklaring av lesingen av filene.

## Detaljer om bruk

Det enkleste eksempelet på hvordan disse klassene fungerer vist i dette eksempelet hvor vi simulerer banen til Jorden med en sirkelbane rundt sola i 10 år med Velocity Verlet metoden:
```
int n = 100000; // antall tidssteg
double totalTime = 10; // total tid
double dt = totalTime / n; // størrelse på tidssteg
SolarSystem solarSystem("exampleOutput.txt"); // SolarSystem objekt med output fil
vVerlet solver(dt); // vVerlet objekt med tidssteg dt
solarSystem.addBody(vec3(0, 0, 0), vec3(0, 0, 0), 1, "Sun"); // legger til solen
solarSystem.addBody(vec3(1, 0, 0), vec3(0, 2 * M_PI, 0), 5.97219/ 1988500, "Earth"); // legger til jorden
solarSystem.writeSystem(); // skriver navn og masser til planeter til fil
for (int i = 0; i < n; i++) {
    solarSystem.writePosVel(); // skriver posisjoner og hastigheter til jorda og sola til fil
    solver.step(solarSystem); // finner akselerasjoner og oppdaterer posisjoner og hastigheter
}
```
Hvis du vil at noen av objektene ikke skal flyttes i løpet av simuleringen må du legge disse til til slutt og la det siste parameteret være "true" i addBody() for at de skal bli ignorert i oppdateringen, og for at posisjonen og hastigheten deres ikke skal blir inkludert i output filen. Simuleringen med relativistisk akselerasjon har blitt gjort rask, men også lite fleksibel. Man må legge til planeten først, så solen for at det skal fungere.
