#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import seaborn as sns
from scipy import stats
import warnings
import os
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
font_path = os.path.join('/insomnia001/home/qc2358/.local/share/fonts', 'arial.ttf')
fm.fontManager.addfont(font_path)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.family'] = 'Arial'

#%% Parallel processing
from multiprocessing import Pool
import multiprocessing as mp
from tqdm import tqdm

num_cores = mp.cpu_count() - 1 # Leave one core for system
print(f"Number of cores: {num_cores}")

#%% Load data files
# Load dose-response genes
print("\n1. Loading dose-response gene data...")
dose_response = pd.read_csv('T_cell_dose_gene_with_score.csv')
print(f"   Shape: {dose_response.shape}")
print(f"   Columns: {list(dose_response.columns)}")

# Load gene copy number variations syn69962258
print("\n2. Loading gene copy number data...")
# copy_number = pd.read_csv('raw/variants.gene_copy_number.csv')
# print(f"   Shape: {copy_number.shape}")
# print(f"   Columns: {list(copy_number.columns)}")
copy_number = pd.read_csv('copy_number_parsed.csv')
print(f"   Shape: {copy_number.shape}")
print(f"   Columns: {list(copy_number.columns)}")

# Load gene expression data syn69961520
print("\n3. Loading gene expression matrix...")
gene_expression = pd.read_csv('raw/gene_tpm_matrix_all_samples.tsv', sep='\t')
print(f"   Genes: {gene_expression.shape[0]}, Samples: {gene_expression.shape[1] - 1}")

# Load clinical data syn69931127
print("\n4. Loading clinical data...")
clinical = pd.read_csv('raw/clinical.surgeries.csv')
print(f"   Shape: {clinical.shape}")
print(f"   Columns: {list(clinical.columns)}")

# Load variant annotations (for IDH status) syn69962263
print("\n5. Loading variant annotations...")
variants_anno_sample = pd.read_csv('raw/variants.anno.csv')
print(f"   Sample shape: {variants_anno_sample.shape}")
print(f"   Columns: {list(variants_anno_sample.columns)}")

#%% Auxiliary functions for sample ID parsing
def parse_sample_id(sample_id):

    parts = sample_id.split('-')
    case_barcode = '-'.join(parts[:3])
    sample_type = parts[3]
    tech = parts[5]
    code = parts[6]

    return case_barcode, sample_type, tech, code

#%% Parse sample IDs for copy number data
with Pool(num_cores) as pool:
    results = list(tqdm(pool.imap(parse_sample_id, copy_number['aliquot_barcode']), 
                       total=len(copy_number['aliquot_barcode']),
                       desc="Parsing sample IDs"))
    
parsed_data = pd.DataFrame(results, columns=['case_barcode', 'sample_type', 'tech', 'code'])
copy_number[['case_barcode', 'sample_type', 'tech', 'code']] = parsed_data

print(f"   Unique cases: {copy_number['case_barcode'].nunique()}")
print(copy_number['sample_type'].value_counts())
print(copy_number['tech'].value_counts())
copy_number.to_csv('copy_number_parsed.csv', index=False)

#%% Identify PDGFRA Copy Number Status
# Extract PDGFRA copy number data
pdgfra_cn = copy_number[copy_number['gene_symbol'] == 'PDGFRA'].copy()
print(f"\n1. PDGFRA copy number records: {len(pdgfra_cn)}")

# Check for conflicting calls
print("\n2. Checking for conflicting copy number calls...")
conflicts = pdgfra_cn.groupby(['case_barcode', 'sample_type'])['hlvl_call'].agg(['nunique', 'count'])
conflicting_cases = conflicts[conflicts['nunique'] > 1]
print(f"   Cases with conflicting calls: {len(conflicting_cases)}")

# Get cases with conflicts
conflicts = pdgfra_cn.groupby(['case_barcode', 'sample_type'])['hlvl_call'].nunique()
conflicting_cases = conflicts[conflicts > 1].index
conflicting_cases


#%% Resolve conflicts
# Resolve conflicts by taking the most extreme alteration
# pdgfra_status = pdgfra_cn.groupby(['case_barcode', 'sample_type']).agg({
#     'hlvl_call': lambda x: x.loc[x.abs().idxmax()],  # Most extreme alteration
#     'wcr': 'mean',  # Average weighted copy ratio
#     'aliquot_barcode': 'first'
# }).reset_index()

# Filter out conflicting cases
pdgfra_cn_filtered = pdgfra_cn[~pdgfra_cn.set_index(['case_barcode', 'sample_type']).index.isin(conflicting_cases)]

# Update pdgfra_status to use filtered data
pdgfra_status = pdgfra_cn_filtered.groupby(['case_barcode', 'sample_type']).agg({
    'hlvl_call': 'first',
    'wcr': 'mean',
    'aliquot_barcode': 'first'
}).reset_index()

print(f"Unique case-sample combinations after filtering: {len(pdgfra_status)}")
print(pdgfra_status['hlvl_call'].value_counts().sort_index())

#%% Create integrated dataset with clinical info
clinical['sample_type'] = clinical['sample_barcode'].str.split('-').str[3]

# Check for duplicate entries
duplicate_check = clinical.groupby(['case_barcode', 'sample_type']).size().reset_index(name='count')
duplicates = duplicate_check[duplicate_check['count'] > 1]
duplicates # No duplicates found

#%% Integrate clinical and molecular data
integrated_data = pdgfra_status.merge(
    clinical[['case_barcode', 'sample_type', 'idh_status', 'codel_status', 'grade',
              'histology', 'mgmt_methylation']],
    on=['case_barcode', 'sample_type'],
    how='left'
)

# integrated_data['sample_category'] = integrated_data['sample_type'].apply(
#     lambda x: 'Primary' if x == 'TP' else ('Recurrent' if str(x).startswith('R') else 'Other')
# )

#%% Select top dose-response genes
# Filter for significant genes
sig_threshold = 0.05
sig_genes = dose_response[dose_response['q_value'] < sig_threshold].copy()

# Select top genes by effect size
n_top_genes = 100
top_upregulated = sig_genes.nlargest(n_top_genes, 'normalized_effect')
top_downregulated = sig_genes.nsmallest(n_top_genes, 'normalized_effect')

# Extract gene names for expression analysis
up_gene_names = top_upregulated['gene_short_name'].tolist()
print(up_gene_names)

down_gene_names = top_downregulated['gene_short_name'].tolist()
print(down_gene_names)

#%% All expression data
gene_expression.set_index('Gene_symbol', inplace=True)
expression_matrix = gene_expression.T
expression_matrix.index.name = 'sample_id'
expression_matrix.reset_index(inplace=True)
expression_matrix['sample_id'] = expression_matrix['sample_id'].str.replace('.', '-')

parsed_data = pd.DataFrame([parse_sample_id(sid) for sid in expression_matrix['sample_id']], 
                         columns=['case_barcode', 'sample_type', 'tech', 'code'])
expression_matrix[['case_barcode', 'sample_type', 'tech', 'code']] = parsed_data

expression_matrix.head(5)
expression_data = expression_matrix.merge(
    integrated_data[['case_barcode', 'sample_type',
                     'hlvl_call', 'wcr', 'idh_status',
                     'grade', 'histology', 'mgmt_methylation']],
    on=['case_barcode', 'sample_type'],
    how='inner'
)
expression_data

#%%
up_gene_names = gene_expression.index[gene_expression.index.isin(up_gene_names)]
down_gene_names = gene_expression.index[gene_expression.index.isin(down_gene_names)]

print(len(up_gene_names))
print(len(down_gene_names))

#%% Auxiliary functions for signature score calculation
def calculate_signature_scores(df, gene_list):
    expr = np.log2(df[gene_list] + 1)
    return expr.mean(axis=1)

def gsva_score(df, gene_sets):
    import gseapy
    
    # Get all genes from gene sets
    gene_list = []
    for gene_set in gene_sets.values():
        gene_list.extend(gene_set)
    
    expr = np.log2(df[gene_list] + 1)
    es = gseapy.gsva(expr.T, gene_sets)
    return es.res2d

#%% Calculate signature scores
up_score = calculate_signature_scores(expression_data, up_gene_names)
down_score = calculate_signature_scores(expression_data, down_gene_names)
expression_data['up_score'] = up_score.astype(float)
expression_data['down_score'] = down_score.astype(float)

#%% Calculate GSVA scores
gsva_results = gsva_score(expression_data, {'Up-regulated': up_gene_names, 'Down-regulated': down_gene_names})
gsva_results['Name'] = gsva_results['Name'].astype(int)

up_regulated_gsva = gsva_results[gsva_results['Term'] == 'Up-regulated']
up_regulated_gsva.set_index('Name', inplace=True)
up_regulated_gsva.index.name = None
up_regulated_gsva.sort_index(ascending=True, inplace=True)

down_regulated_gsva = gsva_results[gsva_results['Term'] == 'Down-regulated']
down_regulated_gsva.set_index('Name', inplace=True)
down_regulated_gsva.index.name = None
down_regulated_gsva.sort_index(ascending=True, inplace=True)

expression_data['up_gsva_score'] = up_regulated_gsva['ES'].astype(float)
expression_data['down_gsva_score'] = down_regulated_gsva['ES'].astype(float)

expression_data['up_gsva_score'] = expression_data['up_gsva_score']
expression_data['down_gsva_score'] = expression_data['down_gsva_score']

#%%
## PLOTTING SCORES
#%%
def plot_scores(df, score_type='signature'):

    from statannotations.Annotator import Annotator

    # Validate score type
    if score_type not in ['signature', 'gsva']:
        raise ValueError("score_type must be either 'signature' or 'gsva'")

    title = {'signature': 'Signature', 'gsva': 'GSVA'}

    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(5, 7))

    # Define color palette and labels
    palette = {'Loss / Deletion': '#ef233c', 'Copy neutral': '#edf2f4', 'Gain / Amplification': '#8d99ae'}

    # Set score columns based on type
    up_col = 'up_score' if score_type == 'signature' else 'up_gsva_score'
    down_col = 'down_score' if score_type == 'signature' else 'down_gsva_score'
    # Get unique status values and create all pairwise combinations
    status_values = df['PDGFRA_status'].unique()
    pairs = [(status_values[i], status_values[j]) 
            for i in range(len(status_values))
            for j in range(i+1, len(status_values))]

    # Count samples per status for labels
    status_counts = df['PDGFRA_status'].value_counts()
    status_labels = [f'{status}\n(n = {status_counts[status]})' for status in status_values]

    # Plot for up-regulated genes
    sns.boxplot(data=df, x='PDGFRA_status', y=up_col, ax=ax1,
                palette=palette, showfliers=False)
    sns.stripplot(data=df, x='PDGFRA_status', y=up_col,
                 ax=ax1, color='black', size=3)


    annotator = Annotator(ax1, pairs, data=df, x='PDGFRA_status', y=up_col, order=df['PDGFRA_status'].unique())
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', line_width=1, line_height=0.01, fontsize=10)
    annotator.apply_and_annotate()

    ax1.grid(False)
    ax1.set_xlabel('PDGFRA Status', fontsize=12, fontfamily='Arial')
    ax1.set_ylabel(f'Up-regulated {title[score_type]} Score', fontsize=12, fontfamily='Arial')
    ax1.set_xticklabels(status_labels,
                        fontsize=8, fontfamily='Arial')
    ax1.tick_params(axis='both', which='major', labelsize=10, length=5, left=True, bottom=True)
    for label in ax1.get_yticklabels():
        label.set_fontfamily('Arial')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(True)
    ax1.spines['bottom'].set_color('black')
    ax1.spines['left'].set_color('black')

    # Plot for down-regulated genes
    sns.boxplot(data=df, x='PDGFRA_status', y=down_col, ax=ax2,
                palette=palette, showfliers=False)
    sns.stripplot(data=df, x='PDGFRA_status', y=down_col,
                 ax=ax2, color='black', size=3)

    # Perform pairwise Mann-Whitney U tests for down-regulated scores
    annotator = Annotator(ax2, pairs, data=df, x='PDGFRA_status', y=down_col, order=df['PDGFRA_status'].unique())
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', line_width=1, line_height=0.01, fontsize=10)
    annotator.apply_and_annotate()

    ax2.grid(False)
    ax2.set_xlabel('PDGFRA Status', fontsize=12, fontfamily='Arial')
    ax2.set_ylabel(f'Down-regulated {title[score_type]} Score', fontsize=12, fontfamily='Arial')
    ax2.set_xticklabels(status_labels,
                        fontsize=8, fontfamily='Arial')
    ax2.tick_params(axis='both', which='major', labelsize=10, length=5, left=True, bottom=True)
    for label in ax2.get_yticklabels():
        label.set_fontfamily('Arial')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['bottom'].set_visible(True)
    ax2.spines['left'].set_visible(True)
    ax2.spines['bottom'].set_color('black')
    ax2.spines['left'].set_color('black')

    plt.tight_layout()
    return fig, (ax1, ax2)


#%%
# Combine hlvl_call categories for plotting
los_amp = {-2: 'Loss / Deletion', -1: 'Loss / Deletion', 0: 'Copy neutral', 1: 'Gain / Amplification', 2: 'Gain / Amplification'}
expression_data['PDGFRA_status'] = expression_data['hlvl_call'].map(los_amp)

#%%
plot_scores(expression_data, score_type='gsva')
plt.savefig('figure/All_samples_gsva_scores.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/All_samples_gsva_scores.svg', dpi=300, bbox_inches='tight')

#%%
plot_scores(expression_data, score_type='signature')
plt.savefig('figure/All_samples_signature_scores.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/All_samples_signature_scores.svg', dpi=300, bbox_inches='tight')

#%%
subset = expression_data.query('histology == "Glioblastoma"')
subset['hlvl_call'].value_counts()

#%%
plot_scores(subset, score_type='gsva')
plt.savefig('figure/Glioblastoma_samples_gsva_scores.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/Glioblastoma_samples_gsva_scores.svg', dpi=300, bbox_inches='tight')

#%%
plot_scores(subset, score_type='signature')
plt.savefig('figure/Glioblastoma_samples_signature_scores.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/Glioblastoma_samples_signature_scores.svg', dpi=300, bbox_inches='tight')

#%%
def plot_gene_exp(df, gene_name):
    
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    # Define color palette and labels
    palette = {'Loss / Deletion': '#ef233c', 'Copy neutral': '#edf2f4', 'Gain / Amplification': '#8d99ae'}

    df['PDGFRA_status'] = df['hlvl_call'].map(los_amp)
    sns.boxplot(data=df, x='PDGFRA_status', y=np.log2(df[gene_name] + 1), ax=ax,
                palette=palette, showfliers=False)
    sns.stripplot(data=df, x='PDGFRA_status', y=np.log2(df[gene_name] + 1),
                 ax=ax, color='black', size=3)
    
    # Count samples per status for labels
    status_counts = df['PDGFRA_status'].value_counts()
    status_values = df['PDGFRA_status'].unique()
    status_labels = [f'{status}\n(n = {status_counts[status]})' for status in status_values]
    
    ax.set_xlabel('PDGFRA Status', fontsize=12, fontfamily='Arial')
    ax.set_ylabel(f'{gene_name} expression ($log_2(TPM + 1)$)', fontsize=12, fontfamily='Arial')
    ax.set_xticklabels(status_labels, fontsize=8, fontfamily='Arial')
    ax.tick_params(axis='both', labelsize=10, length=5, left=True, bottom=True)
    for label in ax.get_yticklabels():
        label.set_fontfamily('Arial')

    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    
    return fig, ax


#%%
plot_gene_exp(expression_data, 'SOD2')
plot_gene_exp(expression_data, 'CD274')
plot_gene_exp(expression_data, 'IDO1')
plot_gene_exp(expression_data, 'CD3D')


#%%
def plot_pdgfra_vs_tcell_signature(df, score_type='gsva', 
                                   pdgfra_gene='PDGFRA'):

    from scipy.stats import gaussian_kde
    
    # Color palette (same as plot_scores)
    palette = {'Loss / Deletion': '#ef233c', 'Copy neutral': '#edf2f4', 'Gain / Amplification': '#8d99ae'}
    
    # Get PDGFRA expression and use log2 transformation
    if pdgfra_gene not in df.columns:
        raise ValueError(f"{pdgfra_gene} not found in expression data")
    
    df['log2_pdgfra'] = np.log2(df[pdgfra_gene] + 1)
    up_col = 'up_gsva_score' if score_type == 'gsva' else 'up_score'
    down_col = 'down_gsva_score' if score_type == 'gsva' else 'down_score'

    title = {'gsva': 'GSVA', 'signature': 'Signature'}

    # Use pre-calculated GSVA scores (same as plot_scores)
    if up_col not in df.columns or down_col not in df.columns:
        raise ValueError(f"{up_col} or {down_col} not found in expression_data. Run GSVA calculation first.")
    
    # Filter out any NaN values
    df_plot = df[['log2_pdgfra', up_col, down_col, 'PDGFRA_status']].dropna()
    
    if len(df_plot) == 0:
        raise ValueError("No valid data points after filtering")
    
    # Create figure with two subplots arranged vertically
    fig = plt.figure(figsize=(5, 8))
    
    # Process each subplot (up and down) - arranged vertically
    for idx, (gene_type, y_col, title_suffix) in enumerate([('Up', up_col, 'Up-regulated'), 
                                                              ('Down', down_col, 'Down-regulated')]):
        # Create grid for subplots with marginal plots
        # Vertical arrangement: top subplot (idx=0) at top, bottom subplot (idx=1) at bottom
        # Leave space for legends below each subplot
        if idx == 0:
            # Top subplot - leave space at bottom for legend
            gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], 
                                  hspace=0.05, wspace=0.05,
                                  left=0.1, right=0.95, 
                                  top=0.98, bottom=0.56)
        else:
            # Bottom subplot - leave space at bottom for legend, ensure no overlap with top
            gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], 
                                  hspace=0.05, wspace=0.05,
                                  left=0.1, right=0.95, 
                                  top=0.44, bottom=0.02)
        
        # Main scatter plot
        ax_main = fig.add_subplot(gs[1, 0])
        
        # Plot scatter points for each PDGFRA status
        status_order = ['Loss / Deletion', 'Copy neutral', 'Gain / Amplification']
        for status in status_order:
            subset = df_plot[df_plot['PDGFRA_status'] == status]
            if len(subset) > 0:
                ax_main.scatter(subset['log2_pdgfra'], subset[y_col],
                               c=palette[status], label=f'{status}\n(n = {len(subset)})',
                               alpha=1.0, s=30, edgecolors='black', linewidth=0.5)
        
        # Add regression line with confidence interval (no label)
        x_all = df_plot['log2_pdgfra'].values
        y_all = df_plot[y_col].values
        
        # Fit regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_all, y_all)
        x_line = np.linspace(x_all.min(), x_all.max(), 100)
        y_line = slope * x_line + intercept
        
        # Calculate confidence interval
        n = len(x_all)
        y_pred = slope * x_all + intercept
        residuals = y_all - y_pred
        mse = np.sum(residuals**2) / (n - 2)
        t_crit = stats.t.ppf(0.975, n - 2)
        
        x_mean = x_all.mean()
        ssx = np.sum((x_all - x_mean)**2)
        se_fit = np.sqrt(mse * (1/n + (x_line - x_mean)**2 / ssx))
        ci = t_crit * se_fit
        
        ax_main.plot(x_line, y_line, 'b-', linewidth=2)
        ax_main.fill_between(x_line, y_line - ci, y_line + ci, alpha=0.2, color='gray')

        # Calculate Pearson correlation
        pearson_r, pearson_p = stats.pearsonr(x_all, y_all)
        
        # Add correlation text
        ax_main.text(0.95, 0.95, f'r = {pearson_r:.3f}\np = {pearson_p:.4f}',
                    transform=ax_main.transAxes, fontsize=10, fontfamily='Arial',
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Set labels
        ax_main.set_xlabel('PDGFRA expression ($log_2(TPM + 1)$)', fontsize=12, fontweight='bold', fontfamily='Arial')
        ax_main.set_ylabel(f'{title_suffix} {title[score_type]} Score', 
                          fontsize=12, fontweight='bold', fontfamily='Arial')
        
        # Legend at bottom (only for status, not regression line)
        handles, labels = ax_main.get_legend_handles_labels()
        # Adjust legend position to prevent overlap
        if idx == 0:
            # Top subplot - legend below, closer to axis
            legend = ax_main.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.16), 
                                   ncol=3, fontsize=10, frameon=False)
        else:
            # Bottom subplot - legend below, closer to axis
            legend = ax_main.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.16), 
                                   ncol=3, fontsize=10, frameon=False)
        for text in legend.get_texts():
            text.set_fontfamily('Arial')
        
        ax_main.tick_params(axis='both', labelsize=10, length=5, left=True, bottom=True)
        for label in ax_main.get_xticklabels():
            label.set_fontfamily('Arial')
        for label in ax_main.get_yticklabels():
            label.set_fontfamily('Arial')
        ax_main.grid(False)  # No grid
        ax_main.spines['top'].set_visible(False)
        ax_main.spines['right'].set_visible(False)
        ax_main.spines['bottom'].set_visible(True)
        ax_main.spines['left'].set_visible(True)
        ax_main.spines['bottom'].set_color('black')
        ax_main.spines['left'].set_color('black')

        # Top marginal plot (PDGFRA expression density)
        ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
        x_range = np.linspace(df_plot['log2_pdgfra'].min(), df_plot['log2_pdgfra'].max(), 200)
        for status in status_order:
            subset = df_plot[df_plot['PDGFRA_status'] == status]
            if len(subset) > 1:
                try:
                    kde = gaussian_kde(subset['log2_pdgfra'])
                    density = kde(x_range)
                    ax_top.plot(x_range, density, color=palette[status], linewidth=2, alpha=0.8)
                    ax_top.fill_between(x_range, density, alpha=0.3, color=palette[status])
                except:
                    ax_top.hist(subset['log2_pdgfra'], bins=30, density=True, alpha=0.6,
                               color=palette[status], edgecolor='black', linewidth=0.5)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['left'].set_visible(False)
        ax_top.tick_params(labelleft=False, labelbottom=False)
        ax_top.grid(False)
        
        # Right marginal plot (T cell signature expression density)
        ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
        y_range = np.linspace(df_plot[y_col].min(), df_plot[y_col].max(), 200)
        for status in status_order:
            subset = df_plot[df_plot['PDGFRA_status'] == status]
            if len(subset) > 1:
                try:
                    kde = gaussian_kde(subset[y_col])
                    density = kde(y_range)
                    ax_right.plot(density, y_range, color=palette[status], linewidth=2, alpha=0.8)
                    ax_right.fill_betweenx(y_range, density, alpha=0.3, color=palette[status])
                except:
                    ax_right.hist(subset[y_col], bins=30, density=True, 
                                 orientation='horizontal', alpha=0.6,
                                 color=palette[status], edgecolor='black', linewidth=0.5)
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)
        ax_right.spines['bottom'].set_visible(False)
        ax_right.spines['left'].set_visible(False)
        ax_right.tick_params(labelleft=False, labelbottom=False)
        ax_right.grid(False)
        
        # Remove top-right subplot (empty)
        ax_empty = fig.add_subplot(gs[0, 1])
        ax_empty.axis('off')
    
    plt.tight_layout()
    return fig

#%%
# Plot PDGFRA expression vs aggregated T cell dose-response gene expression
fig = plot_pdgfra_vs_tcell_signature(expression_data, score_type='gsva')
plt.savefig('figure/All_samples_gsva_vs_pdgfra.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/All_samples_gsva_vs_pdgfra.svg', dpi=300, bbox_inches='tight')

#%%
fig = plot_pdgfra_vs_tcell_signature(expression_data, score_type='signature')
plt.savefig('figure/All_samples_signature_vs_pdgfra.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/All_samples_signature_vs_pdgfra.svg', dpi=300, bbox_inches='tight')

#%%
# Plot PDGFRA expression vs aggregated T cell dose-response gene expression
fig = plot_pdgfra_vs_tcell_signature(subset, score_type='gsva')
plt.savefig('figure/Glioblastoma_samples_gsva_vs_pdgfra.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/Glioblastoma_samples_gsva_vs_pdgfra.svg', dpi=300, bbox_inches='tight')

#%%
fig = plot_pdgfra_vs_tcell_signature(subset, score_type='signature')
plt.savefig('figure/Glioblastoma_samples_signature_vs_pdgfra.png', dpi=300, bbox_inches='tight')
plt.savefig('figure/Glioblastoma_samples_signature_vs_pdgfra.svg', dpi=300, bbox_inches='tight')