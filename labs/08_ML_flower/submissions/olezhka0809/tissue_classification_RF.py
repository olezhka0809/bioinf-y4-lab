import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder

import matplotlib.pyplot as plt
import seaborn as sns


# =========================
# CONFIG
# =========================
CSV_FILE = "tissue_gene_expression.csv"
RANDOM_STATE = 42

# Dacă vrei să echilibrezi datasetul prin downsampling:
USE_DOWNSAMPLING = False   # True / False
# Dacă USE_DOWNSAMPLING=True, fiecare clasă va fi redusă la mărimea celei mai mici (Kidney).
# =========================


def downsample_to_min_class(df: pd.DataFrame, label_col: str, random_state: int = 42) -> pd.DataFrame:
    """Downsample fiecare clasă la mărimea celei mai mici clase."""
    rng = np.random.RandomState(random_state)
    counts = df[label_col].value_counts()
    min_n = counts.min()

    parts = []
    for cls, n in counts.items():
        cls_df = df[df[label_col] == cls]
        parts.append(cls_df.sample(n=min_n, random_state=rng))

    out = pd.concat(parts, axis=0).sample(frac=1, random_state=rng).reset_index(drop=True)
    return out


def main():
    print("1) Încarc datasetul...")
    data = pd.read_csv(CSV_FILE)

    if "Tissue_Type" not in data.columns:
        raise ValueError("CSV-ul trebuie să conțină coloana 'Tissue_Type'.")

    print("   Shape:", data.shape)
    print("   Distribuție clase inițială:")
    print(data["Tissue_Type"].value_counts(), "\n")

    #Creez un demo dataset pentru exemplu
    
    demo_df = (
        data
        .groupby("Tissue_Type", group_keys=False)
        .apply(lambda x: x.sample(n=10, random_state=42))
    )

    demo_df.to_csv("tissue_gene_expression_demo.csv", index=False)

    # Verifică valori lipsă
    na_count = int(data.isna().sum().sum())
    if na_count > 0:
        print(f"[AVERTISMENT] Datasetul are {na_count} valori lipsă. Le elimin (dropna).")
        data = data.dropna().reset_index(drop=True)

    # Downsampling (opțional)
    if USE_DOWNSAMPLING:
        print("2) Fac downsampling pentru echilibrare clase...")
        data = downsample_to_min_class(data, "Tissue_Type", random_state=RANDOM_STATE)
        print("   Shape după downsampling:", data.shape)
        print("   Distribuție clase după downsampling:")
        print(data["Tissue_Type"].value_counts(), "\n")
    else:
        print("2) Nu folosesc downsampling (voi folosi class_weight='balanced').\n")

    # Separă features și label
    X = data.drop(columns=["Tissue_Type"])
    y = data["Tissue_Type"]

    # Encode labels
    label_encoder = LabelEncoder()
    y_encoded = label_encoder.fit_transform(y)
    class_names = label_encoder.classes_
    print("3) Clase (ordine LabelEncoder):", list(class_names), "\n")

    # Split stratificat (important la dezechilibru)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, test_size=0.2, random_state=RANDOM_STATE, stratify=y_encoded
    )

    print("4) Train/Test split:")
    print("   Train:", X_train.shape, " Test:", X_test.shape, "\n")

    # GridSearchCV
    print("5) GridSearchCV pentru RandomForest (scoring=f1_macro)...")
    param_grid = {
        "n_estimators": [100, 200, 400],
        "max_depth": [None, 10, 20, 30],
        "min_samples_split": [2, 5, 10]
    }

    base_rf = RandomForestClassifier(
        random_state=RANDOM_STATE,
        class_weight=None if USE_DOWNSAMPLING else "balanced",
        n_jobs=-1
    )

    grid = GridSearchCV(
        estimator=base_rf,
        param_grid=param_grid,
        cv=5,
        scoring="f1_macro",
        n_jobs=-1,
        verbose=1
    )

    grid.fit(X_train, y_train)
    best_model = grid.best_estimator_

    print("\n   Best params:", grid.best_params_)
    print("   Best CV f1_macro:", grid.best_score_, "\n")

    # Predict & Evaluation
    print("6) Evaluare pe test set (model optimizat)...")
    y_pred = best_model.predict(X_test)

    print("Classification Report:")
    print(classification_report(y_test, y_pred, target_names=class_names))

    # Confusion matrix
    conf_matrix = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        conf_matrix, annot=True, fmt='d', cmap='Blues',
        xticklabels=class_names, yticklabels=class_names
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Confusion Matrix (Random Forest - Best Model)")
    plt.tight_layout()
    plt.show()

    # Feature importance (top 20 pentru vizual)
    importances = best_model.feature_importances_
    idx_sorted = np.argsort(importances)[::-1]
    top_n = min(20, len(importances))
    top_features = X.columns[idx_sorted[:top_n]]
    top_importances = importances[idx_sorted[:top_n]]

    plt.figure(figsize=(10, 6))
    plt.barh(top_features[::-1], top_importances[::-1])
    plt.xlabel("Feature Importance")
    plt.title(f"Top {top_n} Feature Importances (Best Model)")
    plt.tight_layout()
    plt.show()

    # Retrain cu Top 10 features
    print("7) Retraining cu Top 10 features...")
    top10 = X.columns[idx_sorted[:10]]
    print("   Top 10 genes/features:", list(top10), "\n")

    X_train_top = X_train[top10]
    X_test_top = X_test[top10]

    rf_top = RandomForestClassifier(
        random_state=RANDOM_STATE,
        n_estimators=best_model.n_estimators,
        max_depth=best_model.max_depth,
        min_samples_split=best_model.min_samples_split,
        class_weight=None if USE_DOWNSAMPLING else "balanced",
        n_jobs=-1
    )
    rf_top.fit(X_train_top, y_train)
    y_pred_top = rf_top.predict(X_test_top)

    print("Classification Report (Top 10 features retrain):")
    print(classification_report(y_test, y_pred_top, target_names=class_names))

    conf_matrix_top = confusion_matrix(y_test, y_pred_top)
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        conf_matrix_top, annot=True, fmt='d', cmap='Blues',
        xticklabels=class_names, yticklabels=class_names
    )
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title("Confusion Matrix (Random Forest - Top 10 Features)")
    plt.tight_layout()
    plt.show()

    print("\n✔ Gata. Ai: GridSearchCV + raport + confusion matrix + feature importances + retrain top10.\n")


if __name__ == "__main__":
    main()
