import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.neural_network import MLPClassifier
import seaborn as sns
import matplotlib.patches as mpatches
import os

def main():
    """
    Main function to execute the machine learning workflow:
    1. Load and preprocess data
    2. Perform feature selection
    3. Train a Neural Network model
    4. Visualize the network architecture
    """
    
    # --- 1. Configuration and Data Loading ---
    # Set working directory and load dataset
    input_file_path = r"C:\Users\xieru\Desktop\ADHD\代谢组\python.csv"
    
    if not os.path.exists(input_file_path):
        raise FileNotFoundError(f"Input file not found: {input_file_path}")
    
    os.chdir(os.path.dirname(input_file_path))
    df = pd.read_csv(os.path.basename(input_file_path))
    
    print(f"Loaded dataset with shape: {df.shape}")

    # --- 2. Data Preprocessing & Feature Selection ---
    # Extract features (X) and labels (y)
    X = df.iloc[:, 1:].values
    y = df.iloc[:, 0].values
    feature_names_all = df.columns[1:].tolist()

    # Select top 20 most important features using ANOVA F-value
    SELECTED_FEATURE_COUNT = 20
    if X.shape[1] < SELECTED_FEATURE_COUNT:
        raise ValueError(f"Dataset has only {X.shape[1]} features, but need {SELECTED_FEATURE_COUNT}.")

    selector = SelectKBest(f_classif, k=SELECTED_FEATURE_COUNT)
    X_selected = selector.fit_transform(X, y)

    # Get the names and scores of the selected features
    selected_indices = selector.get_support(indices=True)
    selected_feature_names = [feature_names_all[i] for i in selected_indices]
    selected_scores = selector.scores_[selected_indices]

    # Save feature importance ranking to CSV
    importance_df = pd.DataFrame({
        'Feature': selected_feature_names,
        'Importance_Score': selected_scores
    }).sort_values('Importance_Score', ascending=False)
    
    importance_df.to_csv('selected_features.csv', index=False)
    print(f"Selected features saved to CSV.")
    
    # --- 3. Model Training ---
    # Standardize the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_selected)

    # Split into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.3, random_state=42, stratify=y
    )

    # Initialize and train the Multi-Layer Perceptron (Neural Network)
    model = MLPClassifier(
        hidden_layer_sizes=(10, 5),
        activation='relu',
        solver='adam',
        alpha=0.0001,
        batch_size='auto',
        learning_rate_init=0.001,
        max_iter=1000,
        random_state=42
    )
    model.fit(X_train, y_train)

    # --- 4. Visualization ---
    # Visualize the Neural Network architecture
    plot_neural_network(model, importance_df)
    print("Neural network structure plot saved successfully.")


def plot_neural_network(model, importance_df):
    """
    Generates and saves a schematic diagram of the trained neural network.
    
    Args:
        model: Trained MLPClassifier model.
        importance_df: DataFrame containing feature names and their importance scores.
    """
    
    # Extract model architecture details
    n_input = len(importance_df)
    hidden_layers = model.hidden_layer_sizes
    n_output = 1  # Binary classification

    # Define layer sizes and positions for plotting
    layer_sizes = [n_input] + list(hidden_layers) + [n_output]
    layer_positions = np.linspace(0.1, 0.9, len(layer_sizes))  # Horizontal spacing

    # Create the plot
    plt.figure(figsize=(16, 12))
    ax = plt.gca()
    
    # Color palette
    COLORS = {
        'input_top': '#1f77b4',   # Dark Blue for Top features
        'input_rest': '#aec7e8',  # Light Blue for other features
        'hidden1': '#ff9896',     # Salmon
        'hidden2': '#ffbb78',     # Orange
        'output': '#2ca02c'       # Green
    }

    # --- Draw Nodes ---
    max_nodes_per_layer = max(layer_sizes)
    
    for layer_idx, (n_nodes, layer_x) in enumerate(zip(layer_sizes, layer_positions)):
        # Vertical positioning of nodes
        v_spacing = 0.8 / (n_nodes - 1) if n_nodes > 1 else 0.5
        node_ys = [0.1 + i * v_spacing for i in range(n_nodes)]

        for node_idx, node_y in enumerate(node_ys):
            # Determine node color based on layer and importance (for input layer)
            if layer_idx == 0:  # Input Layer
                feature_name = importance_df.iloc[node_idx]['Feature']
                # Highlight Top 10 features
                if node_idx < 10: 
                    color = COLORS['input_top']
                    radius = 0.03
                else:
                    color = COLORS['input_rest']
                    radius = 0.025
            elif layer_idx == len(layer_sizes) - 1:  # Output Layer
                color = COLORS['output']
                radius = 0.04
            elif layer_idx == 1:  # First Hidden Layer
                color = COLORS['hidden1']
                radius = 0.035
            else:  # Second Hidden Layer (if exists)
                color = COLORS['hidden2']
                radius = 0.035

            # Draw the node (circle)
            circle = plt.Circle((layer_x, node_y), radius, color=color, ec='black', lw=0.5)
            ax.add_patch(circle)

            # Add labels (only for input and output)
            if layer_idx == 0 and node_idx < 15:  # Label first 15 input nodes
                label_x = layer_x - 0.08
                # Bold font for top features
                weight = 'bold' if node_idx < 10 else 'normal'
                plt.text(label_x, node_y, feature_name[:20] + '...', 
                         va='center', ha='right', fontsize=9, fontweight=weight)
            
            elif layer_idx == len(layer_sizes) - 1:  # Output label
                plt.text(layer_x + 0.05, node_y, 'Disease Risk\n(0 or 1)', 
                         va='center', ha='left', fontsize=10, fontweight='bold')

    # --- Draw Connections (Edges) ---
    # Calculate min/max importance for transparency scaling
    min_imp = importance_df['Importance_Score'].min()
    max_imp = importance_df['Importance_Score'].max()
    imp_range = max_imp - min_imp if max_imp != min_imp else 1

    for layer_idx in range(len(layer_sizes) - 1):
        current_layer_size = layer_sizes[layer_idx]
        next_layer_size = layer_sizes[layer_idx + 1]
        
        current_x = layer_positions[layer_idx]
        next_x = layer_positions[layer_idx + 1]
        
        # Vertical positions for current and next layer nodes
        current_ys = np.linspace(0.1, 0.9, current_layer_size) if current_layer_size > 1 else [0.5]
        next_ys = np.linspace(0.1, 0.9, next_layer_size) if next_layer_size > 1 else [0.5]

        for i, y1 in enumerate(current_ys):
            for j, y2 in enumerate(next_ys):
                # For input layer -> hidden layer: adjust alpha based on feature importance
                if layer_idx == 0:
                    imp_score = importance_df.iloc[i]['Importance_Score']
                    alpha = 0.1 + 0.6 * ((imp_score - min_imp) / imp_range) # Scale alpha by importance
                else:
                    alpha = 0.15  # Default alpha for other connections
                
                plt.plot([current_x, next_x], [y1, y2], 'k-', alpha=alpha, lw=1)

    # --- Styling and Legend ---
    # Add layer description text
    for pos, size, name in zip(layer_positions, layer_sizes, 
                              ['Input\n(Neurons)', 'Hidden\nLayer 1', 'Hidden\nLayer 2', 'Output\n(Node)']):
        plt.text(pos, 0.95, f"{name}\n({size})", ha='center', va='top', 
                 fontsize=11, fontweight='bold', bbox=dict(facecolor='white', alpha=0.5))

    # Create Legend
    legend_elements = [
        mpatches.Patch(color=COLORS['input_top'], label='Top 10 Features (Bold)'),
        mpatches.Patch(color=COLORS['input_rest'], label='Other Selected Features'),
        mpatches.Patch(color=COLORS['hidden1'], label='Hidden Layer 1 (10 neurons)'),
        mpatches.Patch(color=COLORS['hidden2'], label='Hidden Layer 2 (5 neurons)'),
        mpatches.Patch(color=COLORS['output'], label='Output (Prediction)')
    ]
    
    plt.legend(handles=legend_elements, loc='upper center', 
               bbox_to_anchor=(0.5, 0.02), ncol=3, fontsize=10)

    # Final plot setup
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.title('Neural Network Architecture for Disease Prediction', fontsize=20, pad=20)

    # Save the figure
    output_pdf = "neural_network_structure.pdf"
    output_png = "neural_network_structure.png"
    
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    plt.close()

    return os.path.abspath(output_pdf)


# --- Execute the script ---
if __name__ == "__main__":
    main()
