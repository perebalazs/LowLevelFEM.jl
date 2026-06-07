# Kontakt

- [ ] physical group csomópontjainak lekérése (lásd constrainedDoFs függvény)

- [ ] adott csomópontban a mező értékeinek átmásolása

- [ ] `a`-ban tároljuk a mező értéket a kisebb dimenziójú objektumon, `numElem`-ben a csomópontok sorszámát (pl. `VectorField`)

- lehet hogy sparse vektorként kell tárolni hogy lehessen szorozni a transzformációs mátrixszal

- de az is lehet hogy ha lekérdezzük a felület normálvektorait, akkor azokkal is el lehet végezni egy transzformációt.

- lehet hogy egy külön structot kellene létrehozni a kontakt felületnek, és saját műveleteket definiálni a skalár szorzásokhoz, forgatásokhoz, stb.

- kell:
  
  - tenzormezőt vektormezővel szorozni
  
  - vektormezőt vektormezővel (pl. normálvektorok) szorozni
  
  - vektormezőket összeadni
  
  - vektormezőket skalár számmal szorozni

# Mező műveletek

- [ ] spaceToPlane (3D --> 2D vektormező konvertálás)

- [ ] planeToSpace (2D --> 3D vektormező konvertálás)

- [ ] distrubuted To Nodal (megoszló mező redukálása csomópontokba)

- [ ] `type`-ok rendberakása - minimális számú type (:s, :e, :2D, :3D, :scalar)

- [ ] minden szorzásnál és leszűkíteni az elemek számát a nem nullákra

- [ ] minden összeadásnál kibővíteni az elemek számát a legbővebb nem nullára
