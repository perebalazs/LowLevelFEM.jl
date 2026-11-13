# LowLevelFEM – kód-review (solve, eigen, transient)

*Az új csomag (LLFEM.zip) áttekintése a korábbi verzióhoz képest.*

Megnéztem a diffet a `solve.jl`, `eigen.jl`, `transient.jl` és `transient_modify.jl` fájlokban, a `LowLevelFEMIterativeSolversExt.jl`-t, a `Project.toml`-t, a tesztek szerkezetét és a `test_transient.jl` kulcstesztjeit. A többi forrásfájlt (`contact`, `linear`, `nonlinear`, `spectral`, a `general/` és `multifield/` almappák) nem olvastam végig. Juliát nem futtattam, a numerikus állításaimat Python-replikákkal ellenőriztem.

Ez a verzió lényegesen jobb, mint az előző.

## Visszavont gyanú

A HHT-ban a rögzített DOF-ok tehetetlenségi csatolását (`Mfc*Ad`) most az `(1+α)`/`−α` súlyozású erőtagba teszed. Ezt hibásnak gondoltam, mert a tömegtag a szabványos HHT-ban súlyozatlan. Lemodelleztem egy 2 szabadságfokú rendszeren, előírt `sin(3t)` mozgással. A súlyozott változat másodrendűen konvergál (a hiba 7,7e-4-ről 1,2e-5-re esik, a rend 2,00). A súlyozatlan változat csak elsőrendű (rend ≈ 0,98), α=−0,1-nél és −0,3-nál is. A jelenlegi megoldás tehát a helyes.

## Erősségek

- **A rögzített DOF-ok gyorsulása** most központi differenciákból jön, csak a szabad ismeretlenek használnak Newmark-rekurziót. Ez megszünteti azt a hibát, hogy a gyorsulási hiba nem konvergált (a replikámban a jel 12-szerese volt).
- **A rögzített DOF-ok története** kompakt `nfixed × n` tárolású, MPC esetén is a reprezentatív DOF-okon dolgozik. Ez megoldja a nagy memóriaproblémát.
- **A `SymmetricSystemMatrix` kezelése** (`system_matrix_data`): ez volt a nyitott kérdésem a háromszög-tárolásról, és kiderült, hogy valódi hiba volt. A koordináta-transzformáció a fél mátrixot vetítette át, és ez csendben hibás eredményt adott.
- **Az iteratív solver** `LinearSolveConvergenceError`-t dob az igazi `norm(A*x-b)` reziduummal, és az `on_failure=:warn` explicit kikapcsolás. A nulla főátlós elemre `SingularException` jön.
- **A HHT paraméterei:** az `α`-ból automatikusan a szabványos γ és β következik, a régi viselkedés `δ=0`-val érhető el.
- **A CDM `:central` csillapítása** helyes: `M + Δt/2·C` faktorizálva, a rögzített DOF-oknál konzisztens. A késleltetett csillapítás stabilitási képlete (`2/ω·(√(1+ξ²)−ξ)`) a Jury-kritériumból levezetve stimmel.
- **Az `eigen.jl`:** a csupasz `catch` már `PosDefException`-re szűr, van `dense_limit` védelem, a sajátvektorok fázisa rögzített, a `shift` küszöb dokumentáltan nem merevtest-tolerancia, és újrapróbálkozó ciklus kezeli a túl kevés módust.
- **A tesztelés és a karbantartás:** 3700 sor teszt, külön `review_regressions` fájl, API-dokumentáció-ellenőr (`check_api_docs.py`), precompile-latencia teszt, opcionális függőségek weakdeps-ként. Ez ritkán látott fegyelem.

## Gyengeségek

### 1. A tesztek nem mérnek konvergenciarendet

- A mozgó-támaszos teszt `t²` mozgást használ, amit bármely Newmark-séma pontosan visszaad. Egy hibás deriváltkezelés itt nem bukna meg.
- A "weighted equilibrium" teszt a megvalósítás képletét számolja újra a tesztben, vagyis önmagát igazolja.
- Hiányzik egy Δt-felezéses teszt szinuszos támasszal és α≠0-val. Pontosan ez fogta volna meg a korábbi hibát. A replikám alapján egy ilyen teszt másodrendet vár.
- A tesztek 1–8 elemes hálókon futnak, és nem láttam rajtuk valódi analitikus benchmarkot (rúd hullámterjedése, konzol sajátfrekvenciái). Ezt hedge-eléssel mondom, mert nem néztem meg mindegyiket.

### 2. Nincs kimenet-ritkítás

- A `HHT`/`CDM` mindig teljes `ndof × n` mezőt ad vissza `u`-ra és `v`-re. A `reconstruct_transient_solution` közben további két-három ekkora átmeneti tömböt épít (`X`, `Xfull`, `Q*Xfull`).
- Nagy modellnél, hosszú futásnál ez a fal. Kellene egy `saveEvery=` vagy `outputDofs=` opció, és lehetőség arra, hogy a `v`-t ne adja vissza.
- A rögzített DOF-ok memóriaproblémáját megoldottad, a kimeneté marad.

### 3. Az iteratív solver mindig nulláról indul

Az `IterativeSolvers.cg(A, b)` zéró kezdőértékkel fut, miközben a transient lépésekben az előző gyorsulás kiváló kezdőbecslés lenne (`cg!` a megfelelő hívás). Ez az iterációszámot jelentősen csökkentené.

### 4. Az `_transient_boundary_data` `isapprox`-a abszolút tolerancia nélkül dolgozik

MPC mellett időfüggő BC-nél a `0` és az `1e-16` összehasonlítása hamis "incompatible prescribed values" hibát adhat (pl. `sin(π)`). Kell egy `atol`.

### 5. Az `eigen.jl` két megmaradt gyengesége

- A szimmetria-ellenőrzés még mindig teljes Cholesky-faktorizáció, amit az ARPACK utána újra elvégez.
- A fázisnormálásnál az `argmax(abs.(v))` egyenlő nagyságú komponenseknél (szimmetrikus módusok) kerekítési zajtól függ, tehát ott az előjel instabil maradhat.
- A `shift`-újrapróbálkozás minden körben újra faktorizálja az `A − σB`-t.

### 6. A `transient_modify.jl` még mindig drága és sűrű

- Minden iterációban `_transient_snapshot` épít Field-objektumokat és `p.L*z` szorzást végez, még akkor is, ha a callback lineáris. Lineáris callbacknél legalább két kiértékelés jut egy lépésre.
- A `copyOperators=true` alapértelmezés minden callbacknél teljes mátrixmásolatot jelent. A `false` mellékút jó, de az alapértelmezés a lassú.
- A nemlineáris iterációban nincs lépésvágás vagy adaptív Δt: konvergenciahibánál a szimuláció leáll.
- A stílus javult, de még vegyes (`z[p.free, 1] .=free_values`, üres, szóközös sorok a `_modified_hht` közepén). Egy JuliaFormatter-lépés a CI-ben rendet tenne. CI-konfigurációt a zip-ben nem láttam, lehet, hogy csak nem csomagoltad bele.

### 7. Maradékok és következetlenségek

- Az `FDM_old` alias még ott van.
- Az `FDM(K, C, q, X0, …)` argumentumai nincsenek típusozva, míg a `HHT`/`CDM` igen.
- A CDM alapértelmezése még mindig `:lagged` és `:backward`, vagyis a gyengébb séma. A visszafelé kompatibilitás érthető, de a következő major verzióban megfordítanám.
- A `_project_system_matrix` egységmátrix `T` esetén magát az `A`-t adja vissza, nem másolatot. A dokumentált "mátrixot nem módosítjuk" garancia ezen áll, a `copyOperators=false` mód mellett is.

### 8. Az architektúra méretéből adódó kockázat

- 50 000 sor egy modulban, és több óriásfájl: `linear.jl` 5400, `contact.jl` 5000, `nonlinear.jl` 4200 sor.
- A régi API megmaradt (`∇_old`, `legacy/` tesztek), tehát két felületet kell karbantartani.
- A `general.jl` és a `multifield.jl` szétbontása jó irány. Ugyanezt a `linear`/`nonlinear`/`contact` fájloknál is érdemes folytatni.

## Összkép

A korábbi kritikus pontjaimból a lényegesek javultak: a numerikus hiba és a memóriaprobléma megszűnt, az iteratív solver már nem fogad el konvergálatlan megoldást. A számszerű mag most már megbízhatóbb, mint a tesztek, amik ezt védenék. A következő lépés szerintem az, hogy a tesztek konvergenciarendet és analitikus benchmarkokat mérjenek, utána pedig a kimenet-ritkítás és a transient-iteráció teljesítménye.
