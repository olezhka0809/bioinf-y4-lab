# Lab 08 – Supervised Learning: Tissue Classification (GTEx)

## Dataset demo petru comoditate

Datasetul complet conține 3631 probe și este prea mare pentru a fi inclus integral.
Am inclus un subset demonstrativ, echilibrat între clase, generat reproductibil (seed=42), exclusiv pentru inspecție și exemplificare.

## Dataset utilizat

Pentru acest laborator am utilizat un **subset din GTEx v10 (bulk RNA-seq)**, conținând expresia genică la nivel de genă (read counts), obținută cu RNASeQC v2.4.2.

Datasetul final (`tissue_gene_expression.csv`) are următoarele caracteristici:
- **3631 probe**
- **100 gene** selectate aleator
- **3 clase distincte de țesut**:
  - Brain (3234 probe)
  - Liver (282 probe)
  - Kidney (115 probe)

Datasetul este dezechilibrat, cu o clasă dominantă (Brain).

## Preprocesare și splitarea datelor

- Nu s-a aplicat downsampling.
- Dezechilibrul dintre clase a fost tratat folosind parametrul `class_weight="balanced"` în Random Forest.
- Etichetele de țesut au fost codificate numeric cu `LabelEncoder`, în ordinea:
  - Brain
  - Kidney
  - Liver

Împărțirea datelor s-a realizat prin **train/test split stratificat**:
- **Train:** 2904 probe
- **Test:** 727 probe

Această abordare asigură păstrarea proporțiilor claselor în ambele seturi.

## Model și optimizare

Pentru clasificare am utilizat un **Random Forest Classifier**.

Optimizarea hiperparametrilor a fost realizată cu **GridSearchCV**, folosind:
- validare încrucișată cu 5 fold-uri,
- metrica de evaluare **F1-macro**, potrivită pentru seturi dezechilibrate.

Parametrii testați:
- `n_estimators`: 100, 200, 400
- `max_depth`: None, 10, 20, 30
- `min_samples_split`: 2, 5, 10

### Cei mai buni parametri obținuți:
- `n_estimators = 200`
- `max_depth = None`
- `min_samples_split = 2`

Scorul mediu obținut în cross-validation:
- **F1-macro = 1.00**

## Rezultate pe setul de test

Modelul optimizat a obținut performanță perfectă pe setul de test:

- **Accuracy:** 1.00
- **F1-macro:** 1.00

Pentru fiecare clasă (Brain, Liver, Kidney), precision, recall și F1-score au fost egale cu 1.00.  
Matricea de confuzie indică o clasificare corectă pentru toate probele din setul de test.

## Analiza importanței genelor

Random Forest permite evaluarea importanței caracteristicilor.  
Cele mai importante **10 gene** identificate de model sunt:

- ENSG00000104760.17  
- ENSG00000240935.7  
- ENSG00000186529.16  
- ENSG00000144645.15  
- ENSG00000161082.13  
- ENSG00000118322.14  
- ENSG00000113296.14  
- ENSG00000131044.18  
- ENSG00000172031.7  
- ENSG00000006125.18  

Aceste gene contribuie cel mai mult la separarea claselor de țesut.

## Reantrenare folosind Top 10 gene

Modelul a fost reantrenat folosind **doar cele 10 gene cele mai importante**.

Rezultatul obținut:
- **Accuracy = 1.00**
- **F1-macro = 1.00**

Performanța perfectă se menține chiar și după reducerea drastică a dimensionalității, indicând faptul că un număr foarte mic de gene conține suficient semnal biologic pentru clasificarea țesuturilor analizate.

## Reflecție

Performanța foarte ridicată obținută nu indică supraînvățare, deoarece:
- evaluarea a fost realizată pe un set de test separat,
- hiperparametrii au fost selectați prin validare încrucișată,
- rezultatele sunt consistente și după reducerea numărului de gene.

Scorul perfect reflectă **separabilitatea biologică puternică** a țesuturilor Brain, Liver și Kidney în datele GTEx. Aceste țesuturi prezintă profiluri de expresie genică fundamental diferite, ceea ce face problema relativ ușoară din punct de vedere al clasificării.

Totuși, pentru țesuturi mai apropiate biologic, performanța ar putea scădea semnificativ, fiind necesare modele mai fine sau analize suplimentare.

În concluzie, acest laborator demonstrează eficiența Random Forest pentru clasificarea tipului de țesut pe baza expresiei genice și evidențiază importanța interpretării critice a rezultatelor obținute.
