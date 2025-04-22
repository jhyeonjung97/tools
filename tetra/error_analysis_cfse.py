import os
import re
import socket
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Server address setup
hostname = socket.gethostname()
user_name = os.getlogin()
if hostname == 'PC102616':
    root = '/Users/jiuy97/Desktop/7_V_bulk/figures'
elif user_name == 'jiuy97':
    root = '/pscratch/sd/j/jiuy97/7_V_bulk/figures'
elif user_name == 'hailey' or user_name == 'root':
    root = '/Users/hailey/Desktop/7_V_bulk/figures'
else:
    raise ValueError(f"Unknown hostname: {hostname}. Please set the root path manually.")

def extract_info_from_filename(filename):
    """Extract model type and correlation threshold from filename"""
    # Filename pattern: form_pred_modelTypeThreshold_all_result.log
    # Example: form_pred_lgb70_all_result.log -> model_type: lgb, threshold: 0.7
    pattern = r'form_pred_cfse_([a-z]+)(\d+)_all_result\.log'
    match = re.match(pattern, filename)
    
    if match:
        model_type = match.group(1)
        if match.group(2) == '00':
            threshold = 100
        else:
            threshold = int(match.group(2))
        return model_type, threshold
    return None, None

def extract_metrics_from_file(filepath):
    """Extract test MAE and MSE values from log file"""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            # Find MAE and MSE values in Test Metrics section
            mae_pattern = r'Test Metrics:\n.*?\nMAE: ([\d.]+)'
            mse_pattern = r'Test Metrics:\n.*?\nMSE: ([\d.]+)'
            
            mae_match = re.search(mae_pattern, content, re.DOTALL)
            mse_match = re.search(mse_pattern, content, re.DOTALL)
            
            mae = float(mae_match.group(1)) if mae_match else None
            mse = float(mse_match.group(1)) if mse_match else None
            
            return mae, mse
    except Exception as e:
        print(f"Error reading file {filepath}: {str(e)}")
    return None, None

def main():
    # List to store results
    results = []
    
    # Check all files in figures folder
    for filename in os.listdir(root):
        if filename.startswith('form_pred_cfse_') and filename.endswith('_all_result.log'):
            model_type, threshold = extract_info_from_filename(filename)
            if model_type and threshold:
                mae, mse = extract_metrics_from_file(os.path.join(root, filename))
                if mae is not None and mse is not None:
                    results.append({
                        'model': model_type,
                        'threshold': threshold,
                        'mae': mae,
                        'mse': mse
                    })
    
    # Convert results to DataFrame and sort
    df = pd.DataFrame(results)
    
    if df.empty:
        print("No data available for analysis.")
        return
    
    # Sort DataFrame by model and threshold
    df_sorted = df.sort_values(['model', 'threshold'])
    
    # Save as TSV with rounded values (4 decimal places)
    df_tsv = df_sorted.copy()
    df_tsv['mae'] = df_tsv['mae'].round(4)
    df_tsv['mse'] = df_tsv['mse'].round(4)
    tsv_path = os.path.join(root, 'error_analysis_cfse.tsv')
    df_tsv.to_csv(tsv_path, sep='\t', index=False)
    print(f"Results saved as TSV (rounded to 4 decimal places): {tsv_path}")
    
    # Save as CSV with full precision
    csv_path = os.path.join(root, 'error_analysis_cfse.csv')
    df_sorted.to_csv(csv_path, index=False)
    print(f"Results saved as CSV (full precision): {csv_path}")
    
    # Read the sorted data back for plotting
    df_plot = pd.read_csv(tsv_path, sep='\t')
    
    # Create plots (one for MAE, one for MSE)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot each model type with different colors and markers
    for model in df_plot['model'].unique():
        model_data = df_plot[df_plot['model'] == model].sort_values('threshold')
        
        # MAE plot
        ax1.plot(
            model_data['threshold'], 
            model_data['mae'], 
            marker='o', 
            label=model.upper(),
            linewidth=2
        )
        
        # MSE plot
        ax2.plot(
            model_data['threshold'], 
            model_data['mse'], 
            marker='o', 
            label=model.upper(),
            linewidth=2
        )
    
    # Configure MAE plot
    ax1.set_xlabel('Correlation Threshold')
    ax1.set_ylabel('Test MAE (eV)')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
    ax1.set_xlim(60, 105)
    ax1.set_xticks([70, 80, 90, 100])
    ax1.set_xticklabels(['0.7', '0.8', '0.9', '1.0'])
    
    # Configure MSE plot
    ax2.set_xlabel('Correlation Threshold')
    ax2.set_ylabel('Test MSE (eV²)')
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()
    ax2.set_xlim(60, 105)
    ax2.set_xticks([70, 80, 90, 100])
    ax2.set_xticklabels(['0.7', '0.8', '0.9', '1.0'])
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(root, 'error_analysis_cfse.png'), dpi=300, bbox_inches='tight')
    print(f"Analysis results saved to {os.path.join(root, 'error_analysis_cfse.png')}")

if __name__ == '__main__':
    main()
