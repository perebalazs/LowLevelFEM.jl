# Finite deformation formulation in LowLevelFEM – összefoglaló

Ez az összefoglaló egy **nagy alakváltozású, nemlineáris végeselemes megoldási láncot** ír le a LLFEM keretrendszerben, különös tekintettel:

* a **Total Lagrange** megfogalmazásra,
* az **anyagfüggetlen operátorokra**,
* a **követő (follower) terhelések** konzisztens kezelésére.

A csatolt példa (`Large-deformations_twist.ipynb`) egy egyszerű, de működő demonstrációja ennek az architektúrának.

---

## 1. Alapfilozófia

A LLFEM nagy alakváltozású megoldása az alábbi elvekre épül:

1. **Az anyagtörvény nincs beégetve az operátorokba**
2. **Az operátorok csak kinematikai és integrálási feladatot végeznek**
3. **Minden mező csomóponti**, Gauss-pontba csak integráláskor interpolálunk
4. **A Newton–Raphson linearizáció explicit módon, operátorokra bontva történik**

Ez lehetővé teszi, hogy:

* az anyagtörvény tetszőleges (kézi, AD-alapú, külső csomag),
* a kód alacsony szintű, mégis jól bővíthető maradjon.

---

## 2. Nemlineáris egyenletrendszer

A megoldandó egyenlet:

$$
\boxed{
\mathbf R(\mathbf u)
=
\mathbf f_{\text{int}}(\mathbf P)
-
\mathbf f_{\text{ext}}(\mathbf u)
=
\mathbf 0
}
$$

ahol:

* $ \mathbf u $ : elmozdulás
* $ \mathbf P $ : tetszőleges másodrendű feszültségtensor (értelmezése a felhasználó dolga)

---

## 3. Operátorokra bontott Newton-linearizáció

A Newton-módszerben szükséges tangens:

$$
\boxed{
\mathbf K
=

 \frac{\partial \mathbf R}{\partial \mathbf u}
=
\mathbf K_{\text{mat}}
+
\mathbf K_{\text{geo}}

-

\mathbf K_{\text{ext}}
}
$$

A három tag **külön operátorral** kerül összeállításra.

---

## 4. Belső erő operátor

### `internalForceVector`

* bemenet: `P::TensorField`

* szerep:
  
  $$
  \mathbf f_{\text{int}} = \int_\Omega B_P^T , \mathbf P , d\Omega
  $$

Tulajdonságok:

* stress-measure-agnosztikus (`P`, `S`, `σ` mind elfogadott),
* a feszültség **csomóponti mező**, Gauss-pontban interpolált,
* Total Lagrange típusú integrálás.

---

## 5. Anyagi (konstitutív) tangens

### `materialTangentMatrix`

* bemenet:
  
  * `F::TensorField`
  * `C::6×6` mátrix (Number vagy ScalarField elemekkel)

Szerep:
$$
\mathbf K_{\text{mat}}
=

\int_\Omega B(F)^T , \mathbf C , B(F), d\Omega
$$

Tulajdonságok:

* az anyagtörvény **teljesen külső**,

* `C` lehet:
  
  * konstans,
  * térben változó,

* Green–Lagrange-alapú kinematika,

* Mandel/Voigt jelölés.

Ez az operátor teszi lehetővé, hogy:

* egyszerű anyagtörvények gyorsan kipróbálhatók legyenek,
* később implicit vagy redukált tangensre lehessen váltani.

---

## 6. Geometriai (initial stress) tangens

### `initialStressMatrix`

* bemenet: `stress::TensorField`

Szerep:

$$
\mathbf K_{\text{geo}}
=

\int_\Omega (\nabla N)^T , \mathbf S , (\nabla N), d\Omega
$$

Tulajdonságok:

* kizárólag a meglévő feszültségmezőt használja,
* stress-measure-agnosztikus,
* klasszikus geometriai merevség (buckling, nagy rotáció).

---

## 7. Külső erők és követő terhelések

### Külső erő

A külső erővektor:

* `loadVector(...)`
* lehet **dead load** vagy deformációfüggő (pl. F-fel módosított mező)

### Követő terhelés tangense

### `externalTangentFollower`

Szerep:

$$
\mathbf K_{\text{ext}}
=
\frac{\partial \mathbf f_{\text{ext}}(\mathbf u)}{\partial \mathbf u}
$$

Tulajdonságok:

* kezeli:
  
  * felületi erőket,
  * nyomást,

* figyelembe veszi:
  
  * a terhelés forgását,
  * a Jacobian és $ F^{-T} $ változását,

* Total Lagrange leírásban konzisztens.

Dead load esetén ez a tag **nem jelenik meg**.

---

## 8. Megoldási ciklus (magas szinten)

Egy Newton-iteráció lépései:

1. Elmozdulásból:
   
   * $ \mathbf F $ számítása

2. Anyagtörvény (külső):
   
   * $ \mathbf P $, $ \mathbf C $

3. Operátorok:
   
   * `f_int(P)`
   * `K_mat(F, C)`
   * `K_geo(P)`
   * `K_ext(F)`

4. Newton-lépés:
   
   $$
   \mathbf K \Delta \mathbf u = -\mathbf R
   $$

---

## 9. A csatolt példa szerepe

A `Large-deformations_twist.ipynb` notebook:

* demonstrálja:
  
  * a nagy rotációt,
  * a követő terhelés szükségességét,

* megmutatja:
  
  * mi történik, ha `K_ext` hiányzik,
  * hogyan stabilizálja a megoldást a konzisztens linearizáció,

* szándékosan **kezdetleges**, de:
  
  * matematikailag korrekt,
  * architekturálisan tiszta.

---

## 10. Záró megjegyzés

Ez a felépítés:

* **tankönyvileg korrekt**,

* **kutatásbarát**,

* **nem zárja be** a felhasználót egy anyagtörvénybe,

* és később:
  
  * implicit tangens,
  * Gauss-pontos állapotváltozók,
  * plaszticitás
    irányába is továbbfejleszthető.

Röviden:

> **A LLFEM nem „megoldja” a nagy alakváltozást –
> hanem teret ad arra, hogy helyesen oldd meg.**

Ha szeretnéd, ebből:

* csinálhatunk egy **rövidebb README-verziót**, vagy
* egy **ábrával kísért elméleti összefoglalót** (nagyon ütne a dokumentációban).
