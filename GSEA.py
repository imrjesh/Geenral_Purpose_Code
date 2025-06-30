import os
import pandas as pd
import gseapy as gp
from gseapy.plot import gseaplot
import matplotlib.pyplot as plt

# Set global font size
plt.rcParams.update({
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
})

# Define paths
rnk_path = '/data/kumarr9/DMS273/DESeq2_result'
gmt_file = '/data/kumarr9/Embryo/adult_hepatocyte_2.gmt'

# Define the .rnk files to analyze
rnk_files = ['ic4day_vs_dms273control.rnk', 'dc4day_vs_dms273control.rnk']
rnk_files = [f"{rnk_path}/{file}" for file in rnk_files]

# Dictionary to store results
all_results = {}

# Process each .rnk file
for rnk_file in rnk_files:
    cluster_name = rnk_file.split('/')[-1].replace('_vs_all_GSEA.rnk', '')
    print(f"Processing {cluster_name}...")
    
    # Read the .rnk file
    try:
        ranking = pd.read_csv(rnk_file, sep='\t', header=None, names=['Gene', 'stat'])
    except Exception as e:
        print(f"Failed to read {rnk_file}: {e}")
        continue
    
    # Sort by stat in descending order
    ranking = ranking.sort_values('stat', ascending=False).reset_index(drop=True)
    
    # Define output directory
    outdir = f'{rnk_path}/gsea_output_{cluster_name}'
    os.makedirs(outdir, exist_ok=True)
    
    # Run GSEA prerank with LOWERED THRESHOLDS
    try:
        pre_res = gp.prerank(
            rnk=ranking, 
            gene_sets=gmt_file,  
            seed=6, 
            permutation_num=100,
            min_size=5,       # Lowered from default 15
            max_size=2000,     # Increased from default 500
            outdir=outdir
        )
    except Exception as e:
        print(f"GSEA failed for {cluster_name}: {e}")
        continue
    
    # Extract results with LOWERED P-VALUE THRESHOLD
    out = []
    for term in list(pre_res.results):
        # Include terms with higher p-value threshold (0.1 instead of 0.05)
        if pre_res.results[term]['fdr'] < 0.1:  # Increased FDR threshold
            out.append([
                term,
                pre_res.results[term]['fdr'],
                pre_res.results[term]['es'],
                pre_res.results[term]['nes']
            ])
    
    # Create results dataframe
    out_df = pd.DataFrame(out, columns=['Term', 'fdr', 'es', 'nes']).sort_values('fdr').reset_index(drop=True)
    
    # Store results
    all_results[cluster_name] = {
        'prerank_results': pre_res,
        'summary_df': out_df
    }
    
    print(f"\nAll enriched terms for {cluster_name} (FDR < 0.1):")
    print(out_df)
    
    # Define all expected gene sets
    expected_gene_sets = [
        'Adult_Hepatocyte',
        'Adult_Biliary', 
        'Adult_Hepatobiliary_Progenitors',
        'Fetal_Hepatobiliary',
        'Fetal_Liver'
    ]
    
    # Plot ALL gene sets meeting lowered threshold
    print(f"\nPlotting all gene sets for {cluster_name}:")
    plots_created = 0
    
    for gene_set in expected_gene_sets:
        if gene_set in pre_res.results and pre_res.results[gene_set]['fdr'] < 0.1:
            print(f"Plotting {gene_set}...")
            try:
                gseaplot(
                    pre_res.ranking, 
                    term=gene_set, 
                    **pre_res.results[gene_set]
                )
                fig = plt.gcf()
                fig.set_size_inches(4, 4)
                ax = plt.gca()
                ax.set_title(ax.get_title(), fontsize=24, fontweight='bold', color='black')
                ax.set_xlabel(ax.get_xlabel(), fontsize=24, fontweight='bold', color='black')
                ax.set_ylabel(ax.get_ylabel(), fontsize=24, fontweight='bold', color='black')
                plt.tight_layout()
                
                safe_gene_set = gene_set.replace(' ', '_').replace('/', '_')
                plot_file = f"{outdir}/{cluster_name}_{safe_gene_set}_plot.png"
                plt.savefig(plot_file, bbox_inches='tight', dpi=300)
                plt.close()
                print(f"  ✓ Saved plot to {plot_file}")
                plots_created += 1
                
                result = pre_res.results[gene_set]
                print(f"    FDR: {result['fdr']:.4f}, NES: {result['nes']:.3f}, ES: {result['es']:.3f}")
                
            except Exception as e:
                print(f"  ✗ Plotting failed for {gene_set}: {e}")
        else:
            status = "not found" if gene_set not in pre_res.results else "FDR too high"
            print(f"  ⚠ {gene_set} not plotted ({status})")
    
    print(f"Created {plots_created} plots for {cluster_name}")
    print(f"\n{'='*50}\n")

# Save results to CSV
for cluster_name in all_results:
    output_file = f"{rnk_path}/{cluster_name}_GSEA_results.csv"
    all_results[cluster_name]['summary_df'].to_csv(output_file, index=False)
    print(f"Saved {cluster_name} results to {output_file}")

print("\nSUMMARY:")
print("Modified thresholds applied:")
print("- Minimum gene set size: 5 (default: 15)")
print("- Maximum gene set size: 2000 (default: 500)")
print("- FDR threshold for reporting: < 0.1 (default: < 0.05)")
