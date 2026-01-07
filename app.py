import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
import pandas as pd
from cyvcf2 import VCF
import re
import numpy as np

import requests, sys
import json
import requests
import sys
import json

TARGET_GENE = "SCN1A"
variant = "NM_001165963.4(SCN1A):c.1060G>C"
VCF_FILE = "/path/to/gnomad_ms.BCSQ.PASS.REVEL.am.mpc2.nochr.cadd.vcf.gz"
CLINVAR_VCF = "/path/to/clinvar_plp_ms.BCSQ.cadd.revel.am.mpc2.nochr.vcf.gz"
GENE_COORD = "/path/to/gene_coord.csv.gz"

# first, the gene coordinates
gene_coord = pd.read_csv(GENE_COORD, compression='gzip')
conv_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

# Now filter only those rows that have a non conventional chromosome
gene_coord = gene_coord[gene_coord['Chromosome/scaffold name'].isin(conv_chroms)]

server = "https://rest.ensembl.org"
# Added fields parameter to ensure we get gene_symbol for all transcripts
ext = "/vep/human/hgvs/" + variant + "?mane=1&canonical=1&REVEL=1&AlphaMissense=1&CADD=1&vcf_string=1&fields=transcript_consequences"

r = requests.get(server + ext, headers={"Content-Type": "application/json"})

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()

# The response is a list with one element
data = decoded[0]

# Get ID from colocated variants if available
rsid = data.get('colocated_variants', [{}])[0].get('id', '.')

if data.get('vcf_string'):
    # Split the string by the hyphen and unpack the parts.
    chrom, pos, ref, alt = data['vcf_string'].split('-')
else:
    # Fallback if vcf_string is not present in the API response.
    print("Warning: 'vcf_string' not found in the API result. Unable to determine REF/ALT.")
    ref, alt = None, None

qual = '.'
filter_col = '.'

# --- 2. Construct INFO Field ---
info_parts = []

# a) Build BCSQ tag
bcsq_entries = []
genomic_change = f"{pos}{ref}>{alt}"

for tc in data.get('transcript_consequences', []):
    gene = tc.get('gene_symbol', '.')
    transcript_id = tc.get('transcript_id', '.')
    biotype = tc.get('biotype', '.')

    is_protein_altering = 'amino_acids' in tc

    for term in tc.get('consequence_terms', []):
        base_entry = f"{term}|{gene}|{transcript_id}|{biotype}"
        if term == "missense_variant":
            base_entry = f"missense|{gene}|{transcript_id}|{biotype}"

        if is_protein_altering:
            strand = '+' if tc.get('strand') == 1 else '-'
            aa_pos = tc.get('protein_start')
            ref_aa, alt_aa = tc.get('amino_acids', '/').split('/')
            protein_change = f"{aa_pos}{ref_aa}>{aa_pos}{alt_aa}"
            full_entry = f"{base_entry}|{strand}|{protein_change}|{genomic_change}"
            bcsq_entries.append(full_entry)
        else:
            bcsq_entries.append(base_entry)

if bcsq_entries:
    info_parts.append("BCSQ=" + ",".join(bcsq_entries))

# b) Add other scores using smarter logic

def get_best_score(transcript_consequences, score_key, target_gene):
    """
    Finds the best score from a list of transcript consequences.
    
    Logic:
    1. Collect all available scores for the given score_key.
    2. If there's only one unique score value, return it.
    3. If there are multiple, return the one associated with the target_gene.
    4. If none are found, or none match the target gene in a conflict, return None.
    """
    scores_with_genes = []
    for tc in transcript_consequences:
        if score_key in tc:
            scores_with_genes.append({
                'gene': tc.get('gene_symbol'),
                'score': tc[score_key]
            })

    if not scores_with_genes:
        return None

    # Use a string representation for uniqueness check, handles dicts like AlphaMissense
    unique_scores_str = {json.dumps(s['score']) for s in scores_with_genes}

    if len(unique_scores_str) == 1:
        return scores_with_genes[0]['score']
    else:
        # Conflict: multiple different scores exist. Prioritize the target gene.
        for item in scores_with_genes:
            if item['gene'] == target_gene:
                return item['score']
    
    # Fallback if multiple scores exist but none match the target gene
    return None

all_tcs = data.get('transcript_consequences', [])

cadd_score = get_best_score(all_tcs, 'cadd_phred', TARGET_GENE)
if cadd_score is not None:
    info_parts.append(f"cadd_v1.7={cadd_score}")

revel_score = get_best_score(all_tcs, 'revel', TARGET_GENE)
if revel_score is not None:
    info_parts.append(f"REVEL={revel_score}")

am_score = get_best_score(all_tcs, 'alphamissense', TARGET_GENE)
if am_score is not None:
    info_parts.append(f"am_pathogenicity={am_score.get('am_pathogenicity', '.')}")
    info_parts.append(f"am_class={am_score.get('am_class', '.')}")

info_string = ";".join(info_parts)

# --- 3. Assemble Final VCF Row ---
vcf_row = "\t".join(map(str, [chrom, pos, rsid, ref, alt, qual, filter_col, info_string]))

VCF_ROW_STRING = vcf_row
vcf_string = data['vcf_string'] # define vcf-string, to be displayed in the legend

# get the date of the clinvar vcf file, as written in the header, to be displayed in the legend
import gzip
with gzip.open(CLINVAR_VCF, 'rt') as f:
    for line in f:
        if line.startswith('##fileDate='):
            clinvar_date = line.strip().split('=')[1]
            break

# Initialize the Dash app
app = dash.Dash(__name__)

def get_gene_coordinates(gene_coord_df, gene_symbol):
    """
    Get chromosome and coordinates for a specific gene
    """
    gene_info = gene_coord_df[gene_coord_df['HGNC symbol'] == gene_symbol]
    if gene_info.empty:
        raise ValueError(f"Gene {gene_symbol} not found in gene_coord dataframe")
    gene_info = gene_info.iloc[0]
    chrom = str(gene_info['Chromosome/scaffold name'])
    start = int(gene_info['Gene start (bp)'])
    end = int(gene_info['Gene end (bp)'])
    return chrom, start, end

def to_float(value):
    """Safely convert a value to float, returning np.nan on failure."""
    if value is None:
        return np.nan
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan

def parse_info_field(info_str):
    """Parses a VCF INFO string into a dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def parse_variant_record(variant, gene_symbol, source_type="vcf"):
    """
    Helper function to parse a single variant record (from cyvcf2 or a dictionary).
    Returns a list of variant data dictionaries (one for each relevant transcript).
    """
    parsed_variants = []
    
    if source_type == "vcf":  # From cyvcf2
        info_dict = variant.INFO
        chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, ','.join(variant.ALT)
    else:  # From our custom row parser
        info_dict = variant  # The dict itself is the info
        chrom, pos, ref, alt = variant['chrom'], variant['pos'], variant['ref'], variant['alt']
    
    bcsq_string = info_dict.get('BCSQ')
    if bcsq_string is None:
        return []
    
    for entry in bcsq_string.split(','):
        fields = entry.split('|')
        if len(fields) < 6:
            continue
        
        consequence, entry_gene, entry_transcript, biotype, strand, aa_change = fields[0:6]
        
        if consequence != 'missense' or entry_gene != gene_symbol or not aa_change:
            continue
        
        match = re.match(r'(\d+)', aa_change)
        if not match:
            continue
        
        aa_position = int(match.group(1))
        
        variant_data = {
            'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
            'gene': entry_gene, 'transcript': entry_transcript, 'biotype': biotype,
            'aa_position': aa_position, 'aa_change': aa_change,
            'AC_joint': info_dict.get('AC_joint', 0),
            'AC_genomes': info_dict.get('AC_genomes', 0),
            'nhomalt_joint': info_dict.get('nhomalt_joint', 0),
            'nhomalt_genomes': info_dict.get('nhomalt_genomes', 0),
            'REVEL': info_dict.get('REVEL'),
            'am_pathogenicity': info_dict.get('am_pathogenicity'),
            'cadd_v1.7': info_dict.get('cadd_v1.7'),
            'MPC': info_dict.get('MPC')
        }
        parsed_variants.append(variant_data)
    
    return parsed_variants

def parse_gene_variants_region(vcf_file, gene_coord_df, gene_symbol):
    """
    Parse VCF file for a specific gene using coordinates.
    """
    print(f"Looking up {gene_symbol} coordinates...")
    chrom, start, end = get_gene_coordinates(gene_coord_df, gene_symbol)
    print(f"Found {gene_symbol} at {chrom}:{start}-{end}")
    print(f"Querying VCF region in {vcf_file}...")
    
    variants = []
    vcf = VCF(vcf_file)
    region = f"{chrom}:{start}-{end}"
    
    for variant in vcf(region):
        variants.extend(parse_variant_record(variant, gene_symbol, source_type="vcf"))
    
    print(f"Found {len(variants)} {gene_symbol} missense variants in {vcf_file}")
    return pd.DataFrame(variants)

def parse_vcf_row(vcf_row_string, gene_symbol):
    """
    Parses a single tab-separated VCF row string.
    """
    try:
        print("Parsing custom VCF row string...")
        fields = vcf_row_string.strip().split('\t')
        if len(fields) < 8:
            print("Warning: Custom VCF row is malformed. Skipping.")
            return pd.DataFrame()
        
        chrom, pos, _, ref, alt, _, _, info_str = fields[0:8]
        info_dict = parse_info_field(info_str)
        
        # Add required fields for the parser
        info_dict['chrom'] = chrom
        info_dict['pos'] = int(pos)
        info_dict['ref'] = ref
        info_dict['alt'] = alt
        
        parsed_variants = parse_variant_record(info_dict, gene_symbol, source_type="dict")
        print(f"Found {len(parsed_variants)} missense variants in custom row.")
        return pd.DataFrame(parsed_variants)
    except Exception as e:
        print(f"Error parsing custom VCF row: {e}")
        return pd.DataFrame()

# --- Load all data at startup ---
# Load gnomAD data
try:
    print(f"Loading gnomAD variants from {VCF_FILE}...")
    GNOMAD_DATA = parse_gene_variants_region(VCF_FILE, gene_coord, TARGET_GENE)
    print(f"Loaded {len(GNOMAD_DATA)} gnomAD {TARGET_GENE} variants")
except FileNotFoundError:
    print(f"ERROR: gnomAD VCF file not found at '{VCF_FILE}'. App cannot start.")
    exit()
except Exception as e:
    print(f"An error occurred loading gnomAD VCF: {e}")
    GNOMAD_DATA = pd.DataFrame()

# Load ClinVar data (optional, always displayed when available)
try:
    print(f"Loading ClinVar variants from {CLINVAR_VCF}...")
    CLINVAR_DATA = parse_gene_variants_region(CLINVAR_VCF, gene_coord, TARGET_GENE)
    print(f"Loaded {len(CLINVAR_DATA)} ClinVar {TARGET_GENE} variants")
except FileNotFoundError:
    print(f"INFO: ClinVar VCF file not found at '{CLINVAR_VCF}'. Skipping.")
    CLINVAR_DATA = pd.DataFrame()
except Exception as e:
    print(f"An error occurred loading ClinVar VCF: {e}")
    CLINVAR_DATA = pd.DataFrame()

# Parse custom variant row (optional, always displayed when available)
try:
    CUSTOM_VARIANT_DATA = parse_vcf_row(VCF_ROW_STRING, TARGET_GENE)
    if not CUSTOM_VARIANT_DATA.empty:
        print(f"Loaded custom variant for {TARGET_GENE}")
except Exception as e:
    print(f"An error occurred parsing custom variant row: {e}")
    CUSTOM_VARIANT_DATA = pd.DataFrame()

if GNOMAD_DATA.empty:
    print("WARNING: No gnomAD variants found! The plot will be empty.")
    TRANSCRIPTS = []
else:
    TRANSCRIPTS = sorted(GNOMAD_DATA['transcript'].unique().tolist())
    print(f"Found {len(TRANSCRIPTS)} transcripts: {TRANSCRIPTS}")

# App layout
app.layout = html.Div([
    html.H1(f"{TARGET_GENE} Variant Pathogenicity Analysis", 
            style={'textAlign': 'center', 'marginBottom': 30}),
    
    html.Div([
        html.Div([
            html.Label("Select Transcript:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='transcript-dropdown',
                options=[{'label': t, 'value': t} for t in TRANSCRIPTS],
                value=TRANSCRIPTS[0] if TRANSCRIPTS else None,
                style={'width': '100%'}
            )
        ], style={'width': '30%', 'display': 'inline-block', 'marginRight': '2%'}),
        
        html.Div([
            html.Label("Pathogenicity Score:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='score-dropdown',
                options=[
                    {'label': 'REVEL', 'value': 'REVEL'},
                    {'label': 'AlphaMissense', 'value': 'am_pathogenicity'},
                    {'label': 'CADD v1.7', 'value': 'cadd_v1.7'},
                    {'label': 'MPC2', 'value': 'MPC'}
                ],
                value='REVEL',
                style={'width': '100%'}
            )
        ], style={'width': '30%', 'display': 'inline-block', 'marginRight': '2%'}),
        
        html.Div([
            html.Label("gnomAD Filter Field:", style={'fontWeight': 'bold'}),
            dcc.Dropdown(
                id='threshold-field-dropdown',
                options=[
                    {'label': 'Allele Count (Exomes + Genomes)', 'value': 'AC_joint'},
                    {'label': 'Allele Count (Genomes)', 'value': 'AC_genomes'},
                    {'label': 'Homozygous Count (Exomes + Genomes)', 'value': 'nhomalt_joint'},
                    {'label': 'Homozygous Count (Genomes)', 'value': 'nhomalt_genomes'}
                ],
                value='AC_genomes',
                style={'width': '100%'}
            )
        ], style={'width': '30%', 'display': 'inline-block'}),
    ], style={'marginBottom': 20}),
    
    html.Div([
        html.Label("gnomAD Filter Threshold (show variants with count > this value):", 
                   style={'fontWeight': 'bold'}),
        dcc.Input(
            id='threshold-value',
            type='number',
            value=0,
            min=0,
            step=1,
            style={'width': '200px', 'marginLeft': '10px'}
        ),
        html.Span(" (Set to 0 to show all gnomAD variants)", 
                  style={'marginLeft': '10px', 'fontStyle': 'italic', 'color': '#666'})
    ], style={'marginBottom': 30}),
    
    html.Div(id='info-display', 
             style={'marginBottom': 20, 'padding': '10px', 
                    'backgroundColor': '#f0f0f0', 'borderRadius': '5px'}),
    
    dcc.Graph(id='pathogenicity-plot', style={'height': '600px'}),
    
    html.Div([
        html.H3("Legend:", style={'marginTop': 5}),
        html.Ul([
            html.Li([html.Span("★", style={'color': 'purple', 'fontSize': '20px'}), 
                    " Your Variant: (hg38) "+ vcf_string]),
            html.Li([html.Span("◆", style={'color': 'darkred', 'fontSize': '20px'}), 
                    " ClinVar P/LP Variants (date: " + clinvar_date + ")"]),
            html.Li([html.Span("●", style={'color': '#17BECF', 'fontSize': '20px'}), 
                    " gnomADv4.1 Variants"])
        ])
    ], style={'marginTop': 5, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
              'borderRadius': '5px'})
])

def create_hover_text(df, score_column):
    """Helper function to create hover text"""
    hover_text = []
    for _, row in df.iterrows():
        # Safely get values with defaults
        score_val = row.get(score_column)
        score_str = f"{float(score_val):.3f}" if pd.notna(score_val) else 'N/A'
        
        text = (
            f"<b>Position:</b> {row.get('aa_position', 'N/A')}<br>"
            f"<b>Change:</b> {row.get('aa_change', 'N/A')}<br>"
            f"<b>Variant:</b> {row.get('chrom', '')}:{row.get('pos', '')} "
            f"{row.get('ref', '')}>{row.get('alt', '')}<br>"
            f"<b>{score_column}:</b> {score_str}<br>"
            f"<b>AC_genomes:</b> {row.get('AC_genomes', 0)}<br>"
            f"<b>AC_joint:</b> {row.get('AC_joint', 0)}<br>"
            f"<b>nhomalt_genomes:</b> {row.get('nhomalt_genomes', 0)}<br>"
            f"<b>nhomalt_joint:</b> {row.get('nhomalt_joint', 0)}"
        )
        hover_text.append(text)
    return hover_text

@app.callback(
    Output('pathogenicity-plot', 'figure'),
    Output('info-display', 'children'),
    Input('transcript-dropdown', 'value'),
    Input('score-dropdown', 'value'),
    Input('threshold-field-dropdown', 'value'),
    Input('threshold-value', 'value')
)
def update_plot(transcript, score, threshold_field, threshold_value):
    if not transcript or not score:
        return go.Figure(), "Please select a transcript and pathogenicity score"
    
    # Create figure
    fig = go.Figure()
    info_parts = []
    
    # Trace 1: gnomAD data (filtered dynamically by user)
    try:
        # Convert threshold_value to numeric, default to 0 if None
        threshold_val = float(threshold_value) if threshold_value is not None else 0
        
        # Filter gnomAD data
        gnomad_filtered = GNOMAD_DATA[
            (GNOMAD_DATA['transcript'] == transcript) &
            (GNOMAD_DATA[score].notna())
        ].copy()
        
        # Apply threshold filter
        if not gnomad_filtered.empty and threshold_field in gnomad_filtered.columns:
            # Convert threshold field to numeric
            gnomad_filtered[threshold_field] = pd.to_numeric(gnomad_filtered[threshold_field], errors='coerce')
            gnomad_filtered = gnomad_filtered[gnomad_filtered[threshold_field] > threshold_val]
        
        if not gnomad_filtered.empty:
            gnomad_filtered = gnomad_filtered.sort_values('aa_position')
            
            # Create color scale based on threshold field
            color_values = gnomad_filtered[threshold_field].astype(float).replace(0, 1e-9)
            log_color = np.log10(color_values)
            
            fig.add_trace(go.Scatter(
                x=gnomad_filtered['aa_position'],
                y=gnomad_filtered[score],
                mode='markers',
                name='gnomAD Variants',
                marker=dict(
                    size=8,
                    color=log_color,
                    colorscale='GnBu',
                    showscale=True,
                    colorbar=dict(
                        title=f"log10({threshold_field})",
                        x=1.15
                    ),
                    line=dict(width=1, color='LightGray')
                ),
                text=create_hover_text(gnomad_filtered, score),
                hoverinfo='text',
                hovertemplate='%{text}<extra></extra>'
            ))
            fig.update_traces(
                marker_colorbar=dict(
                    tickvals=np.log10([1, 10, 100, 1000, 10000]),
                    ticktext=['1', '10', '100', '1000', '10000']
                            )
            )
            info_parts.append(f"{len(gnomad_filtered)} gnomAD variants (filtered: {threshold_field} > {threshold_val})")
        else:
            info_parts.append(f"No gnomAD variants with {threshold_field} > {threshold_val}")
    except Exception as e:
        print(f"Error processing gnomAD data: {e}")
        info_parts.append(f"Error loading gnomAD variants: {str(e)}")
    
    # Trace 2: ClinVar data (ALWAYS displayed when available for the transcript)
    if not CLINVAR_DATA.empty:
        try:
            clinvar_filtered = CLINVAR_DATA[
                (CLINVAR_DATA['transcript'] == transcript) &
                (CLINVAR_DATA[score].notna())
            ].copy()
            
            if not clinvar_filtered.empty:
                fig.add_trace(go.Scatter(
                    x=clinvar_filtered['aa_position'],
                    y=clinvar_filtered[score],
                    mode='markers',
                    name='ClinVar P/LP',
                    marker=dict(
                        color='darkred',
                        size=12,
                        symbol='diamond'
                    ),
                    text=create_hover_text(clinvar_filtered, score),
                    hoverinfo='text',
                    hovertemplate='%{text}<extra></extra>'
                ))
                info_parts.append(f"{len(clinvar_filtered)} ClinVar P/LP variants")
        except Exception as e:
            print(f"Error processing ClinVar data: {e}")
            info_parts.append(f"Error loading ClinVar variants: {str(e)}")
    
    # Trace 3: Custom Variant (ALWAYS displayed when available for the transcript)
    if not CUSTOM_VARIANT_DATA.empty:
        try:
            custom_filtered = CUSTOM_VARIANT_DATA[
                (CUSTOM_VARIANT_DATA['transcript'] == transcript) &
                (CUSTOM_VARIANT_DATA[score].notna())
            ].copy()
            
            if not custom_filtered.empty:
                fig.add_trace(go.Scatter(
                    x=custom_filtered['aa_position'],
                    y=custom_filtered[score],
                    mode='markers',
                    name='Custom Variant',
                    marker=dict(
                        color='purple',
                        size=16,
                        symbol='star'
                    ),
                    text=create_hover_text(custom_filtered, score),
                    hoverinfo='text',
                    hovertemplate='%{text}<extra></extra>'
                ))
                info_parts.append(f"1 custom variant")
        except Exception as e:
            print(f"Error processing custom variant data: {e}")
            info_parts.append(f"Error loading custom variant: {str(e)}")
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"<b>{TARGET_GENE} Missense Variants - Pathogenicity Comparison</b><br>"
                 f"<sub>Transcript: {transcript} | Score: {score}</sub>",
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            title="<b>Amino Acid Position</b>",
            gridcolor='lightgray'
        ),
        yaxis=dict(
            title=f"<b>{score} Score</b>",
            gridcolor='lightgray'
        ),
        hovermode='closest',
        plot_bgcolor='white',
        height=600,
        showlegend=False,
        font=dict(size=12),
        margin=dict(r=150)  # Extra margin for colorbar
    )
    
    fig.update_xaxes(showgrid=True, zeroline=False)
    fig.update_yaxes(showgrid=True, zeroline=False)
    
    # Create info text
    if info_parts:
        info_text = "Displaying: " + " | ".join(info_parts)
    else:
        info_text = "No variants to display with current filters."
    
    if len(fig.data) == 0:
        return fig, "⚠️ No variants found for this transcript with the selected score. Check if data is loaded correctly."
    
    return fig, info_text

if __name__ == '__main__':
    app.run(debug=True, port=8050)
