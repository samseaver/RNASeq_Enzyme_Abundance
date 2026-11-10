#!/usr/bin/env python
"""
Figure 4: Control vs FeLim reaction scores, 4 rows x 5 timepoints.

Rows 1-2 are the objective reaction score r_s (Poplar, Sorghum); rows 3-4 the
relative score r_s-tilde. Both axes are log10 -- scores span roughly five
decades -- and the dashed diagonal is the identity line. Colour is I-dist, the
perpendicular distance to that line in log space; the top 5% by |I-dist| are
outlined in black. Column 6 carries the weighted mean and standard error of
I-dist per day.

Data sources
------------
projects/qpsi-plastidial/integration_results/{Poplar,Sorghum}_reaction_scores.tsv
                                             {Poplar,Sorghum}_reaction_molar_fractions.tsv

Outputs
-------
fig_scatter_rslt.png and fig_scatter_rslt.html, beside this script.

Run: micromamba run -n bf-runtime python figures/fig_scatter_rslt.py
"""
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

LOG_FLOOR_QUANTILE = 0.001   # bottom 0.1% of non-zero scores are floored onto the axis


def prepare_data_group(poplar_path, sorghum_path, val_col, dist_space='log'):
    """Pivot, then place both conditions on a log10 axis.

    Reaction scores are log-normal over ~5 decades, so the original
    divide-by-global-max scaling pushed >90% of points into the bottom-left
    corner. Plotting log10(score) spreads them out; the identity line is still
    y = x.

    dist_space controls how I-dist (distance to the identity line) is measured:
      'log'    — perpendicular distance in log space, i.e. a scaled log fold
                 change. Genuinely independent of the size of the scores.
      'linear' — the original (FeLim - Control) / global_max. Ranks, and hence
                 the top-5% outlines, are identical to the published figure.
    """
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

    positive = pd.concat([combined['Control'], combined['FeLim']]).dropna()
    positive = positive[positive > 0]
    lo = float(np.floor(np.log10(positive.quantile(LOG_FLOOR_QUANTILE))))
    hi = float(np.ceil(np.log10(positive.max())))
    floor_val = 10.0 ** lo

    combined['Control_log'] = np.log10(combined['Control'].clip(lower=floor_val))
    combined['FeLim_log'] = np.log10(combined['FeLim'].clip(lower=floor_val))

    if dist_space == 'log':
        combined['signed_dist'] = (combined['FeLim_log'] - combined['Control_log']) / np.sqrt(2)
    else:
        global_max = max(combined['Control'].max(), combined['FeLim'].max())
        combined['signed_dist'] = (combined['FeLim'] - combined['Control']) / global_max

    combined['abs_dist'] = combined['signed_dist'].abs()
    combined['rxn_dist_quantile'] = combined['abs_dist'].rank(pct=True)

    max_dist = combined['abs_dist'].max()
    combined['alpha'] = 0.5 + 0.5 * (combined['abs_dist'] / max_dist) if max_dist > 0 else 0.5
    combined['direction'] = np.where(combined['signed_dist'] >= 0, 'Up', 'Down')

    combined.attrs['log_range'] = (lo, hi)
    return combined

def plot_normalized_dashboard(df_abs, df_rel, out_png, out_html):
    target_days = ['2d', '4d', '7d', '14d', '21d']
    numeric_days = np.array([2, 4, 7, 14, 21])
    
    specs = [[{}, {}, {}, {}, {}, {}]] * 4
    
    subplot_titles = [f"Day {d}" for d in target_days] + [""]
    subplot_titles += [""] * 6 * 3 
            
    fig = make_subplots(
        rows=4, cols=6, specs=specs,
        subplot_titles=subplot_titles,
        horizontal_spacing=0.015, # Massively reduced whitespace
        vertical_spacing=0.035
    )
    
    abs_range = df_abs.attrs['log_range']
    rel_range = df_rel.attrs['log_range']

    row_configs = [
        (df_abs[df_abs['species'] == 'Poplar'], "Poplar Abs", abs_range),
        (df_abs[df_abs['species'] == 'Sorghum'], "Sorghum Abs", abs_range),
        (df_rel[df_rel['species'] == 'Poplar'], "Poplar Rel", rel_range),
        (df_rel[df_rel['species'] == 'Sorghum'], "Sorghum Rel", rel_range)
    ]

    # Colour scale: symmetric, set by the 98th percentile of |I-dist| so the
    # icefire ramp actually spans the bulk of the data instead of a few outliers.
    cmax = float(pd.concat([df_abs['abs_dist'], df_rel['abs_dist']]).quantile(0.98))

    for row_idx, (df, row_title, log_range) in enumerate(row_configs):
        r = row_idx + 1
        lo, hi = log_range
        # Label every other decade — neighbouring panels are close enough that
        # every decade collides. X ticks skip the two edge decades as well,
        # otherwise adjacent panels print their labels on top of each other.
        y_decades = list(range(int(lo), int(hi) + 1, 2))
        x_decades = list(range(int(lo) + 1, int(hi), 2))
        y_text = [f"10<sup>{k}</sup>" for k in y_decades]
        x_text = [f"10<sup>{k}</sup>" for k in x_decades]
        # The Abs rows (1-2) and Rel rows (3-4) span different decades, so the
        # Abs block needs its own x labels rather than borrowing row 4's.
        show_x = (r in (2, 4))
        
        print(f"\n{'='*55}\n STATS FOR: {row_title}\n{'='*55}")
        print(f"{'Day':<5} | {'Dir':<5} | {'PopMean':<10} | {'SE':<10} | {'n':<5}\n{'-'*55}")
        
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
                    x=day_df['Control_log'], y=day_df['FeLim_log'],
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
                x=[lo, hi], y=[lo, hi], mode='lines',
                line=dict(color='rgba(0,0,0,0.5)', width=1.5, dash='dash'),
                hoverinfo='skip', showlegend=False
            ), row=r, col=c)

            # X-Axes: Tick labels ONLY on the bottom row (r == 4)
            fig.update_xaxes(
                range=[lo, hi], tickmode='array', tickvals=x_decades, ticktext=x_text,
                showticklabels=show_x, row=r, col=c
            )
            # Y-Axes: Tick labels ONLY on the far left column (c == 1).
            # scaleanchor="x" forces the subplot to be a perfect square!
            fig.update_yaxes(
                range=[lo, hi], tickmode='array', tickvals=y_decades, ticktext=y_text,
                showticklabels=(c == 1), scaleanchor="x", scaleratio=1, row=r, col=c
            )

            # Row 2 gets tick labels but no axis title — the title would land on
            # top of row 3's panels.
            if r == 4: fig.update_xaxes(title_text="Control", row=r, col=c)
            if c == 1: fig.update_yaxes(title_text=f"<b>{row_title}</b><br>FeLim", row=r, col=c)

        # --- PANEL B: Count-weighted Mean (sum/n_total) + Analytic SE (col 6) ---
        n_total = df.groupby('day').size().reindex(target_days).fillna(0)

        def _weighted(direction):
            grp = df[df['direction'] == direction].groupby('day')['abs_dist']
            n_sub = grp.size().reindex(target_days).fillna(0)
            sd    = grp.std().reindex(target_days).fillna(0)
            denom = n_total.replace(0, np.nan)
            pop_mean = (grp.sum().reindex(target_days).fillna(0) / denom).fillna(0)
            se       = (np.sqrt(n_sub) * sd / denom).fillna(0)
            return pop_mean, se, n_sub

        up_mean, up_se, up_n     = _weighted('Up')
        down_mean, down_se, down_n = _weighted('Down')

        for day in target_days:
            print(f"{day:<5} | Up    | {up_mean[day]:.6f}   | {up_se[day]:.6f}   | {int(up_n[day]):<5}")
            print(f"{'':<5} | Down  | {down_mean[day]:.6f}   | {down_se[day]:.6f}   | {int(down_n[day]):<5}")

        offset = 0.4

        fig.add_trace(go.Scatter(
            x=numeric_days - offset, y=up_mean.values,
            mode='lines+markers', name='Up',
            line=dict(color='chocolate', width=3), marker=dict(symbol='circle', size=10),
            error_y=dict(type='data', array=up_se.values, visible=True, thickness=2.5, width=6, color='chocolate'),
            showlegend=(r == 1)
        ), row=r, col=6)

        fig.add_trace(go.Scatter(
            x=numeric_days + offset, y=down_mean.values,
            mode='lines+markers', name='Down',
            line=dict(color='steelblue', width=3), marker=dict(symbol='square', size=10),
            error_y=dict(type='data', array=down_se.values, visible=True, thickness=2.5, width=6, color='steelblue'),
            showlegend=(r == 1)
        ), row=r, col=6)

        if r == 4: fig.update_xaxes(title_text="Time Point (Day)", row=r, col=6)

        fig.update_xaxes(range=[0, 23], tickmode='array', tickvals=[2, 4, 7, 14, 21], ticktext=['2d', '4d', '7d', '14d', '21d'], showticklabels=(r == 4), row=r, col=6)

        # Scale the trend panel to its own row: Abs and Rel I-dist live on
        # different scales, and a hard-coded range flattens one of them.
        top = float(np.nanmax(np.concatenate([up_mean.values + up_se.values,
                                              down_mean.values + down_se.values])))
        # Round up to a round number of ~4 ticks
        step = 10.0 ** np.floor(np.log10(top / 4.0))
        for mult in (1, 2, 2.5, 5, 10):
            if top / (step * mult) <= 4.5:
                step *= mult
                break
        top = np.ceil(top / step) * step
        ticks = np.arange(0.0, top + step / 2, step)
        fig.update_yaxes(
            range=[-0.05 * top, top],
            tickmode='array',
            tickvals=list(ticks),
            ticktext=[f"{t:g}" for t in ticks],
            side='right',
            row=r, col=6
        )

    # Dashed separator between scatter cols (1-5) and trend col (6),
    # placed at the midpoint of the gap between col 5's right edge and col 6's left edge
    x_separator = (fig.layout.xaxis5.domain[1] + fig.layout.xaxis6.domain[0]) / 2
    # Legend anchored to the upper-left corner of col 6's top-row plot (row 1, col 6 = xaxis6/yaxis6)
    legend_x = fig.layout.xaxis6.domain[0] + 0.005
    legend_y = fig.layout.yaxis6.domain[1] - 0.005

    fig.add_shape(
        type='line', xref='paper', yref='paper',
        x0=x_separator, x1=x_separator,
        y0=0.04, y1=0.96,
        line=dict(color='black', width=1.5, dash='dash')
    )

    # --- Global Publication Layout ---
    fig.update_layout(
        height=1400, width=1950, plot_bgcolor='white',
        margin=dict(t=50, b=50, l=60, r=120),
        font=dict(size=16, color='black', family='Arial'), 
        coloraxis=dict(
            # RdBu reversed: blue = down, red = up, pale in the middle. The old
            # 'icefire' ramp goes black at both ends, which hid the black
            # top-5% outlines once the scale was tightened onto the data.
            colorscale='RdBu_r', cmin=-cmax, cmax=cmax,
            colorbar=dict(
                thickness=20, len=0.7, y=0.5, x=1.045,
                tickmode='array',
                tickvals=[-cmax, -cmax / 2, 0.0, cmax / 2, cmax],
                tickformat=".2f", tickfont=dict(size=16)
            )
        ),
        legend=dict(
            yanchor="top", y=legend_y, xanchor="left", x=legend_x,
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
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(script_dir, "..", "projects", "qpsi-plastidial", "integration_results")

    poplar_abs = os.path.join(base_dir, 'Poplar_reaction_scores.tsv')
    sorghum_abs = os.path.join(base_dir, 'Sorghum_reaction_scores.tsv')
    poplar_rel = os.path.join(base_dir, 'Poplar_reaction_molar_fractions.tsv')
    sorghum_rel = os.path.join(base_dir, 'Sorghum_reaction_molar_fractions.tsv')

    print("\nGrouping and Normalizing Absolute Data (Poplar & Sorghum)...")
    df_abs = prepare_data_group(poplar_abs, sorghum_abs, val_col='reaction_score')

    print("Grouping and Normalizing Relative Data (Poplar & Sorghum)...")
    df_rel = prepare_data_group(poplar_rel, sorghum_rel, val_col='relative_reaction_score')

    out_png = os.path.join(script_dir, "fig_scatter_rslt.png")
    out_html = os.path.join(script_dir, "fig_scatter_rslt.html")

    if not df_abs.empty and not df_rel.empty:
        plot_normalized_dashboard(df_abs, df_rel, out_png, out_html)