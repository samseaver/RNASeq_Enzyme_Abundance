#!/usr/bin/env python
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def load_and_pivot(filepath, value_col):
    df = pd.read_csv(filepath, sep='\t')
    split_cols = df['condition'].str.split('_', expand=True)
    df['tissue'], df['treatment'], df['day'] = split_cols[0], split_cols[1], split_cols[2]
    
    df['day'] = df['day'].replace({'02d': '2d', '04d': '4d', '07d': '7d'})
    df = df[df['tissue'] == 'Leaf']
    df = df[df['treatment'].isin(['Control', 'FeLim'])]
    
    return df.pivot_table(index=['reaction_id', 'day'], columns='treatment', values=value_col, aggfunc='first').reset_index()

def prepare_data_group(poplar_path, sorghum_path, val_col):
    dfs = []
    if os.path.exists(poplar_path):
        df_p = load_and_pivot(poplar_path, val_col)
        df_p['species'] = 'Poplar'
        dfs.append(df_p)
    if os.path.exists(sorghum_path):
        df_s = load_and_pivot(sorghum_path, val_col)
        df_s['species'] = 'Sorghum'
        dfs.append(df_s)

    if not dfs: return pd.DataFrame()
    combined = pd.concat(dfs, ignore_index=True)

    global_max = max(combined['Control'].max(), combined['FeLim'].max())
    combined['Control_norm'] = combined['Control'] / global_max
    combined['FeLim_norm'] = combined['FeLim'] / global_max
    combined['signed_dist'] = combined['FeLim_norm'] - combined['Control_norm']
    combined['abs_dist'] = combined['signed_dist'].abs()
    combined['rxn_dist_quantile'] = combined['abs_dist'].rank(pct=True)
    
    max_dist = combined['abs_dist'].max()
    combined['alpha'] = 0.5 + 0.5 * (combined['abs_dist'] / max_dist) if max_dist > 0 else 0.5
    combined['direction'] = np.where(combined['signed_dist'] >= 0, 'Up', 'Down')
    
    return combined

def plot_normalized_dashboard(df_abs, df_rel, out_png, out_html):
    target_days = ['2d', '4d', '7d', '14d', '21d']
    numeric_days = np.array([2, 4, 7, 14, 21])
    
    specs = [[{}, {}, {}, {}, {}, {'secondary_y': True}]] * 4
    
    subplot_titles = [f"Day {d}" for d in target_days] + [""]
    subplot_titles += [""] * 6 * 3 
            
    fig = make_subplots(
        rows=4, cols=6, specs=specs,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.015, # Massively reduced whitespace
        vertical_spacing=0.035
    )
    
    row_configs = [
        (df_abs[df_abs['species'] == 'Poplar'], "Poplar Abs"),
        (df_abs[df_abs['species'] == 'Sorghum'], "Sorghum Abs"),
        (df_rel[df_rel['species'] == 'Poplar'], "Poplar Rel"),
        (df_rel[df_rel['species'] == 'Sorghum'], "Sorghum Rel")
    ]
    
    for row_idx, (df, row_title) in enumerate(row_configs):
        r = row_idx + 1 
        
        print(f"\n{'='*45}\n STATS FOR: {row_title}\n{'='*45}")
        print(f"{'Day':<5} | {'Dir':<5} | {'Mean':<10} | {'SD':<10}\n{'-'*45}")
        
        # --- PANEL A: Scatter Plots (Cols 1 to 5) ---
        for col_idx, day in enumerate(target_days):
            c = col_idx + 1
            day_df = df[df['day'] == day].sort_values(by='rxn_dist_quantile')
            
            if not day_df.empty:
                quantiles = day_df['rxn_dist_quantile'].tolist()
                line_colors = ['black' if q >= 0.95 else 'rgba(0,0,0,0)' for q in quantiles]
                line_widths = [2.0 if q >= 0.95 else 0 for q in quantiles]
                marker_sizes = [9.0 if q >= 0.95 else 8.0 for q in quantiles] 

                customdata = day_df[['Control', 'FeLim', 'rxn_dist_quantile']].values

                fig.add_trace(go.Scatter(
                    x=day_df['Control_norm'], y=day_df['FeLim_norm'],
                    mode='markers', customdata=customdata,
                    marker=dict(
                        color=day_df['signed_dist'], coloraxis="coloraxis", 
                        opacity=day_df['alpha'], size=marker_sizes,
                        line=dict(color=line_colors, width=line_widths) 
                    ),
                    text=day_df['reaction_id'], showlegend=False,
                    hovertemplate="<b>%{text}</b><br>Mean: %{customdata[0]:.2e}<br>SD: %{customdata[1]:.2e}<extra></extra>"
                ), row=r, col=c)
                
            fig.add_trace(go.Scatter(
                x=[0, 1], y=[0, 1], mode='lines',
                line=dict(color='rgba(0,0,0,0.5)', width=1.5, dash='dash'),
                hoverinfo='skip', showlegend=False
            ), row=r, col=c)

            # X-Axes: Tick labels ONLY on the bottom row (r == 4)
            fig.update_xaxes(
                range=[-0.05, 1.05], tickmode='array', tickvals=[0.0, 0.25, 0.5, 0.75, 1.0], 
                tickformat=".2f", showticklabels=(r == 4), row=r, col=c
            )
            # Y-Axes: Tick labels ONLY on the far left column (c == 1). 
            # scaleanchor="x" forces the subplot to be a perfect square!
            fig.update_yaxes(
                range=[-0.05, 1.05], tickmode='array', tickvals=[0.0, 0.25, 0.5, 0.75, 1.0], 
                tickformat=".2f", showticklabels=(c == 1), scaleanchor="x", scaleratio=1, row=r, col=c
            )
            
            if r == 4: fig.update_xaxes(title_text="Control", row=r, col=c)
            if c == 1: fig.update_yaxes(title_text=f"<b>{row_title}</b><br>FeLim", row=r, col=c)

        # --- PANEL B: Mean + Secondary Axis Error Bars (Col 6) ---
        up_group = df[df['direction'] == 'Up'].groupby('day')['abs_dist']
        up_mean, up_std = up_group.mean().reindex(target_days).fillna(0), up_group.std().reindex(target_days).fillna(0)

        down_group = df[df['direction'] == 'Down'].groupby('day')['abs_dist']
        down_mean, down_std = down_group.mean().reindex(target_days).fillna(0), down_group.std().reindex(target_days).fillna(0)

        for day in target_days:
            print(f"{day:<5} | Up    | {up_mean[day]:.6f}   | {up_std[day]:.6f}")
            print(f"{'':<5} | Down  | {down_mean[day]:.6f}   | {down_std[day]:.6f}")

        K = 5
        offset = 0.4 
        
        # 1. Primary Traces: Just the Mean Lines
        fig.add_trace(go.Scatter(
            x=numeric_days - offset, y=up_mean.values,
            mode='lines+markers', name='Up-reg Mean',
            line=dict(color='chocolate', width=3), marker=dict(symbol='circle', size=10),
            showlegend=(r == 1) 
        ), row=r, col=6, secondary_y=False)
        
        fig.add_trace(go.Scatter(
            x=numeric_days + offset, y=down_mean.values,
            mode='lines+markers', name='Down-reg Mean',
            line=dict(color='steelblue', width=3), marker=dict(symbol='square', size=10),
            showlegend=(r == 1)
        ), row=r, col=6, secondary_y=False)

        # 2. Ghost Traces: Invisible markers mapped to Secondary Y
        fig.add_trace(go.Scatter(
            x=numeric_days - offset, y=up_mean.values * K,
            mode='markers', name='Up-reg SD', marker=dict(opacity=0),
            error_y=dict(type='data', array=up_std.values, visible=True, thickness=2.5, width=6, color='chocolate'),
            showlegend=False 
        ), row=r, col=6, secondary_y=True)

        fig.add_trace(go.Scatter(
            x=numeric_days + offset, y=down_mean.values * K,
            mode='markers', name='Down-reg SD', marker=dict(opacity=0),
            error_y=dict(type='data', array=down_std.values, visible=True, thickness=2.5, width=6, color='steelblue'),
            showlegend=False 
        ), row=r, col=6, secondary_y=True)

        if r == 4: fig.update_xaxes(title_text="Time Point (Day)", row=r, col=6)
        
        fig.update_xaxes(range=[0, 23], tickmode='array', tickvals=[2, 4, 7, 14, 21], ticktext=['2d', '4d', '7d', '14d', '21d'], row=r, col=6)

        # Primary Axis (Mean): 2 Decimal Places
        fig.update_yaxes(
            range=[-0.02, 0.06], 
            tickmode='array', tickvals=[-0.02, 0.00, 0.02, 0.04, 0.06], tickformat=".2f", 
            row=r, col=6, secondary_y=False
        ) 
        # Secondary Axis (SD): 1 Decimal Place
        fig.update_yaxes(
            range=[-0.10, 0.30], 
            tickmode='array', tickvals=[-0.10, 0.00, 0.10, 0.20, 0.30], tickformat=".1f", showgrid=False,
            row=r, col=6, secondary_y=True
        )

    # --- Global Publication Layout ---
    fig.update_layout(
        height=1400, width=1950, plot_bgcolor='white',
        margin=dict(t=50, b=50, l=60, r=60),
        font=dict(size=16, color='black', family='Arial'), 
        coloraxis=dict(
            colorscale='icefire', cmin=-1, cmax=1,      
            colorbar=dict(
                title=dict(text="Signed<br>Norm Dist", font=dict(size=18)),
                thickness=20, len=0.7, y=0.5,
                tickmode='array', tickvals=[-1.0, -0.5, 0.0, 0.5, 1.0], tickformat=".2f", tickfont=dict(size=16)
            )
        ),
        legend=dict(
            yanchor="top", y=1.0, xanchor="left", x=0.83, 
            bgcolor="rgba(255, 255, 255, 0.9)", bordercolor="black", borderwidth=1.5,
            font=dict(size=18)
        )
    )
    
    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=22, color='black')

    # Enforce solid black line boundaries around every subplot
    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='black', mirror=True,
        showgrid=True, gridwidth=1, gridcolor='rgba(211,211,211,0.6)',
        title_font=dict(size=18), tickfont=dict(size=16)
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='black', mirror=True,
        showgrid=True, gridwidth=1, gridcolor='rgba(211,211,211,0.6)',
        title_font=dict(size=18), tickfont=dict(size=16)
    )

    try:
        fig.write_image(out_png, scale=2)
        print(f"\n[*] Saved static PNG: {out_png}")
    except Exception as e:
        print(f"\n[!] Could not save PNG. Error: {e}")
        
    fig.write_html(out_html)
    print(f"[*] Saved interactive HTML: {out_html}")

if __name__ == '__main__':
    base_dir = "/Users/seaver/Seaver_Lab/Git_Repos/RNASeq_Enzyme_Abundance/projects/qpsi-plastidial/integration_results/"
    
    poplar_abs = os.path.join(base_dir, 'Poplar_reaction_scores.tsv')
    sorghum_abs = os.path.join(base_dir, 'Sorghum_reaction_scores.tsv')
    poplar_rel = os.path.join(base_dir, 'Poplar_reaction_molar_fractions.tsv')
    sorghum_rel = os.path.join(base_dir, 'Sorghum_reaction_molar_fractions.tsv')

    print("\nGrouping and Normalizing Absolute Data (Poplar & Sorghum)...")
    df_abs = prepare_data_group(poplar_abs, sorghum_abs, val_col='reaction_score')
    
    print("Grouping and Normalizing Relative Data (Poplar & Sorghum)...")
    df_rel = prepare_data_group(poplar_rel, sorghum_rel, val_col='relative_reaction_score')

    out_png = os.path.join(base_dir, 'std_figures', "Plotly_Master_Dashboard_Normalized.png")
    out_html = os.path.join(base_dir, 'std_figures', "Plotly_Master_Dashboard_Normalized.html")
    
    if not df_abs.empty and not df_rel.empty:
        plot_normalized_dashboard(df_abs, df_rel, out_png, out_html)