### **Rezultate și interpretare**



#### **Performanța modelului**



Modelul Random Forest antrenat pe datasetul expression\_matrix\_olezhka0809\_topvar1000.csv a obținut performanță perfectă pe setul de test:



* Accuracy: 1.00
* F1-score (macro): 1.00
* Precision și Recall: 1.00 pentru toate clasele



Evaluarea a fost realizată pe un set de test separat (20% din date), folosind un split stratificat, ceea ce indică faptul că rezultatele nu sunt obținute prin supraînvățare trivială pe datele de antrenament.



Genele cele mai importante



Top 10 gene identificate de Random Forest pe baza scorului feature\_importances\_ sunt:



1. CTSS
2. GPX1
3. C1orf162
4. S100A9
5. FCER1G
6. GC
7. AGXT
8. APCS
9. FABP1
10. AHSG



Aceste gene contribuie cel mai mult la separarea claselor și reflectă diferențe biologice clare între grupurile identificate prin pseudo-labeling. Mai multe dintre ele sunt cunoscute pentru roluri în:



* procese inflamatorii și imunitare (ex. S100A9, FCER1G, CTSS),
* metabolism hepatic (FABP1, AGXT, AHSG),
* răspunsuri oxidative (GPX1).



Prezența acestor gene sugerează că Random Forest captează semnale biologice relevante și nu doar variații aleatorii din date.



#### **Analiza erorilor de clasificare**



Matricea de confuzie arată că:



* Nu există confuzii între clase pe setul de test.
* Toate probele au fost clasificate corect.



Această separare perfectă indică faptul că:



* clasele (pseudo-labelurile generate prin KMeans) sunt foarte bine separate în spațiul expresiei genice;
* problema de clasificare este relativ ușoară pentru un model non-liniar precum Random Forest, mai ales după selecția genelor cu varianță mare.



