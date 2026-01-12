import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.feature_selection import RFECV
from sklearn.metrics import accuracy_score
import os
import seaborn as sns

# Set font to Arial for publication-quality plots
plt.rcParams['font.family'] = 'Arial'


def main():
    # --- 1. Data Loading ---
    input_file = 'python.csv'
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found in the current directory.")
    
    print(f"Loading dataset: {input_file}")
    df = pd.read_csv(input_file)
    
    # Extract features (X) and labels (y)
    X = df.iloc[:, 1:].values
    y = df.iloc[:, 0].values
    feature_names = df.columns[1:]  # Save feature names for later interpretation
    
    print(f"Dataset shape: {df.shape} (Features: {X.shape[1]}, Samples: {X.shape[0]})")

    # --- 2. Data Preprocessing ---
    # Standardize features to have zero mean and unit variance
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # --- 3. Model Configuration & Feature Selection ---
    # Define Stratified Cross-Validation strategy
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Initialize SVM Classifier (Linear kernel recommended for RFE)
    svm = SVC(kernel='linear', probability=True, random_state=42)

    # Configure RFECV (Recursive Feature Elimination with Cross-Validation)
    rfecv = RFECV(
        estimator=svm,
        step=1,  # Remove one feature at a time
        cv=cv,
        scoring='accuracy',
        min_features_to_select=1,
        n_jobs=-1,  # Use all CPU cores for parallel processing
        verbose=1
    )

    # --- 4. Train and Select Features ---
    print("\nStarting SVM-RFE feature selection...")
    rfecv.fit(X_scaled, y)

    # --- 5. Extract Results ---
    n_features_optimal = rfecv.n_features_
    feature_ranking = rfecv.ranking_
    selected_features_mask = rfecv.support_

    # Create DataFrame to display feature importance
    feature_importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Ranking': feature_ranking,
        'Selected': selected_features_mask
    }).sort_values('Ranking')
    
    print(f"\nOptimal number of features: {n_features_optimal}")
    print("\nFeatures sorted by importance (Rank 1 is best):")
    print(feature_importance_df.head(10))  # Print top 10 for brevity

    # Save feature ranking to CSV
    csv_path = "svm_rfe_feature_importance.csv"
    feature_importance_df.to_csv(csv_path, index=False)
    print(f"\nFeature importance saved to: {csv_path}")

    # --- 6. Plotting Performance Metrics ---
    # Extract cross-validation results
    cv_scores = rfecv.cv_results_['mean_test_score']
    error_rates = 1 - cv_scores
    accuracy_rates = cv_scores
    feature_counts = np.arange(1, len(feature_names) + 1)[-len(cv_scores):]

    # Plot 1: Error Rate
    plot_performance(
        x=feature_counts, 
        y=error_rates, 
        title='SVM-RFE: Error Rate vs Number of Features',
        xlabel='Number of Features',
        ylabel='Error Rate (1 - Accuracy)',
        optimal_point=n_features_optimal,
        optimal_label='Optimal Features (Min Error)',
        color='#E24A33',
        filename='svm_rfe_error_rate'
    )

    # Plot 2: Accuracy
    plot_performance(
        x=feature_counts, 
        y=accuracy_rates, 
        title='SVM-RFE: Accuracy vs Number of Features',
        xlabel='Number of Features',
        ylabel='Accuracy',
        optimal_point=n_features_optimal,
        optimal_label='Optimal Features (Max Accuracy)',
        color='#348ABD',
        filename='svm_rfe_accuracy'
    )

    # --- 7. Final Model Training & Evaluation ---
    # Select only the important features
    X_selected = X_scaled[:, selected_features_mask]
    
    # Split the selected features into train/test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X_selected, y, test_size=0.3, random_state=42, stratify=y
    )

    # Train final model
    final_model = SVC(kernel='linear', probability=True, random_state=42)
    final_model.fit(X_train, y_train)

    # Evaluate performance
    train_acc = accuracy_score(y_train, final_model.predict(X_train))
    test_acc = accuracy_score(y_test, final_model.predict(X_test))
    
    print(f"\nFinal Model Performance with {n_features_optimal} features:")
    print(f"Training Accuracy: {train_acc:.4f}")
    print(f"Testing Accuracy: {test_acc:.4f}")

    # --- 8. Feature Importance Visualization ---
    # For linear SVM, use the absolute value of coefficients as importance
    if hasattr(final_model, 'coef_'):
        importances = np.abs(final_model.coef_[0])
        importance_df = pd.DataFrame({
            'Feature': feature_names[selected_features_mask],
            'Importance': importances
        }).sort_values('Importance', ascending=False)

        # Save importance scores
        importance_df.to_csv("svm_feature_coefficients.csv", index=False)

        # Plot Top N features
        plot_top_features(importance_df, top_n=10)
        # Plot all selected features
        plot_top_features(importance_df, top_n=len(importance_df), filename_suffix="all_selected")


def plot_performance(x, y, title, xlabel, ylabel, optimal_point, optimal_label, color, filename):
    """Helper function to plot error rate or accuracy."""
    plt.figure(figsize=(10, 8))
    plt.plot(x, y, 'o-', color=color, linewidth=2, markersize=6)
    plt.axvline(x=optimal_point, color='gray', linestyle='--', linewidth=2, label=f'{optimal_label}: {optimal_point}')
    
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.title(title, fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(fontsize=12)
    
    # Annotate the optimal point
    if 'Error' in title:
        opt_value = min(y)
        text_offset = 0.05
    else:
        opt_value = max(y)
        text_offset = -0.05
        
    plt.plot(optimal_point, opt_value, 'ro', markersize=8)
    plt.annotate(f'Optimal: {opt_value:.4f}',
                xy=(optimal_point, opt_value), 
                xytext=(optimal_point + 1, opt_value + text_offset),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
                fontsize=12)

    plt.tight_layout()
    
    # Save as both PNG and PDF
    for ext in ['png', 'pdf']:
        plt.savefig(f"{filename}.{ext}", dpi=300, bbox_inches='tight')
    
    print(f"Plot saved: {filename}.png/.pdf")
    plt.show()


def plot_top_features(importance_df, top_n=10, filename_suffix=""):
    """Plot horizontal bar chart of feature importances."""
    plot_df = importance_df.head(top_n) if top_n < len(importance_df) else importance_df
    
    plt.figure(figsize=(10, max(6, len(plot_df) * 0.4)))
    sns.barplot(data=plot_df, x='Importance', y='Feature', palette='viridis')
    
    suffix_text = f" (Top {top_n})" if top_n < len(importance_df) else " (All Selected)"
    plt.title(f'Feature Importance by SVM Coefficients{suffix_text}', fontsize=16)
    plt.xlabel('Absolute Coefficient Value', fontsize=14)
    plt.ylabel('Features', fontsize=14)
    plt.tight_layout()
    
    # Save as both PNG and PDF
    base_name = f"svm_rfe_feature_importance_top{top_n}" if top_n < len(importance_df) else "svm_rfe_all_selected_features"
    for ext in ['png', 'pdf']:
        plt.savefig(f"{base_name}.{ext}", dpi=300, bbox_inches='tight')
    
    print(f"Feature importance plot saved: {base_name}.png/.pdf")
    plt.show()


if __name__ == "__main__":
    main()
