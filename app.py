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

VCF_FILE = "/path/to/gnomad_ms.fully_annotated.vcf.gz"
CLINVAR_VCF = "/path/to/clinvar_plp_ms.fully_annotated.vcf.gz"
GENE_COORD = "/path/to/gene_coord.csv.gz"

# first, the gene coordinates
gene_coord = pd.read_csv(GENE_COORD, compression='gzip')
conv_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']

# Now filter only those rows that have a non conventional chromosome
gene_coord = gene_coord[gene_coord['Chromosome/scaffold name'].isin(conv_chroms)]

import requests
import sys
import json

TARGET_GENE = "kdm2b".upper()
variant = "NM_032590.5:c.1638G>T"
server = "https://rest.ensembl.org"
# Added fields parameter to ensure we get gene_symbol for all transcripts
ext = "/vep/human/hgvs/" + variant + "?mane=1&canonical=1&REVEL=1&AlphaMissense=1&CADD=1&vcf_string=1&fields=transcript_consequences&dbNSFP=BayesDel_noAF_score,VARITY_R_LOO_score,gnomAD_exomes_nhomalt"

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

def parse_score_string(score_value):
    """
    Parses a score that might be a comma-separated string with multiple values.
    Returns the maximum numerical value, or the original value if it's not a string.
    
    Examples:
    - '.,0.083826095,.,.,.,.,.,.' -> 0.083826095
    - '0.5' -> 0.5
    - 0.5 -> 0.5
    """
    if not isinstance(score_value, str):
        return score_value
    
    # Check if it contains commas (multiple values)
    if ',' not in score_value:
        try:
            return float(score_value)
        except (ValueError, TypeError):
            return score_value
    
    # Parse comma-separated values and find the maximum
    values = score_value.split(',')
    numeric_values = []
    
    for val in values:
        val = val.strip()
        if val and val != '.':
            try:
                numeric_values.append(float(val))
            except (ValueError, TypeError):
                continue
    
    if numeric_values:
        return max(numeric_values)
    
    return None

def get_best_score(transcript_consequences, score_key, target_gene):
    """
    Finds the best score from a list of transcript consequences.
    
    Logic:
    1. Collect all available scores for the given score_key.
    2. Parse scores that might be comma-separated strings.
    3. If there's only one unique score value, return it.
    4. If there are multiple, return the one associated with the target_gene.
    5. If none are found, or none match the target gene in a conflict, return None.
    """
    scores_with_genes = []
    for tc in transcript_consequences:
        if score_key in tc:
            raw_score = tc[score_key]
            
            # Special handling for scores that might be comma-separated strings
            if score_key in ['varity_r_loo_score', 'bayesdel_noaf_score']:
                parsed_score = parse_score_string(raw_score)
                if parsed_score is not None:
                    scores_with_genes.append({
                        'gene': tc.get('gene_symbol'),
                        'score': parsed_score
                    })
            else:
                scores_with_genes.append({
                    'gene': tc.get('gene_symbol'),
                    'score': raw_score
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

bayesdel_score = get_best_score(all_tcs, 'bayesdel_noaf_score', TARGET_GENE)
if bayesdel_score is not None:
    info_parts.append(f"BayesDel_nsfp33a_noAF={bayesdel_score}")

varity_score = get_best_score(all_tcs, 'varity_r_loo_score', TARGET_GENE)
if varity_score is not None:
    info_parts.append(f"VARITY_R_LOO={varity_score}")

info_string = ";".join(info_parts)

# --- 3. Assemble Final VCF Row ---
vcf_row = "\t".join(map(str, [chrom, pos, rsid, ref, alt, qual, filter_col, info_string]))
vcf_string = data['vcf_string'] # define vcf-string, to be displayed in the legend

# get the date of the clinvar vcf file, as written in the header, to be displayed in the legend
import gzip
with gzip.open(CLINVAR_VCF, 'rt') as f:
    for line in f:
        if line.startswith('##fileDate='):
            clinvar_date = line.strip().split('=')[1]
            break

# Initialize the Dash app with Bootstrap theme for modal support
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], 
                suppress_callback_exceptions=True)

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
        variant_id = str(variant.ID) if variant.ID else None
    else:  # From our custom row parser
        info_dict = variant  # The dict itself is the info
        chrom, pos, ref, alt = variant['chrom'], variant['pos'], variant['ref'], variant['alt']
        variant_id = variant.get('id')
    
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
            'variant_id': variant_id,
            'gene': entry_gene, 'transcript': entry_transcript, 'biotype': biotype,
            'aa_position': aa_position, 'aa_change': aa_change,
            'AC_joint': info_dict.get('AC_joint', 0),
            'AC_genomes': info_dict.get('AC_genomes', 0),
            'nhomalt_joint': info_dict.get('nhomalt_joint', 0),
            'nhomalt_genomes': info_dict.get('nhomalt_genomes', 0),
            # Convert all scores to float
            'REVEL': to_float(info_dict.get('REVEL')),
            'am_pathogenicity': to_float(info_dict.get('am_pathogenicity')),
            'cadd_v1.7': to_float(info_dict.get('cadd_v1.7')),
            'MPC': to_float(info_dict.get('MPC')),
            'MISTIC_score': to_float(info_dict.get('MISTIC_score')),
            'MISTIC_pred': info_dict.get('MISTIC_pred'),  # Keep as string
            'popEVE': to_float(info_dict.get('popEVE')),
            'BayesDel_nsfp33a_noAF': to_float(info_dict.get('BayesDel_nsfp33a_noAF')),
            'VARITY_R_LOO': to_float(info_dict.get('VARITY_R_LOO')),
            # === CLINVAR SPECIFIC FIELDS ===
            'CLNDN': info_dict.get('CLNDN'),
            'CLNREVSTAT': info_dict.get('CLNREVSTAT'),
            'CLNSIG': info_dict.get('CLNSIG')
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
        
        chrom, pos, variant_id, ref, alt, _, _, info_str = fields[0:8]
        info_dict = parse_info_field(info_str)
        
        # Add required fields for the parser
        info_dict['chrom'] = chrom
        info_dict['pos'] = int(pos)
        info_dict['ref'] = ref
        info_dict['alt'] = alt
        info_dict['id'] = variant_id if variant_id != '.' else None
        
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
    GNOMAD_DATA['source'] = 'gnomAD'
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
    CLINVAR_DATA['source'] = 'ClinVar'
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
    CUSTOM_VARIANT_DATA['source'] = 'Custom'
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
    html.H1(f"missense-visual of {TARGET_GENE}", 
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
                    {'label': 'MPC2', 'value': 'MPC'},
                    {'label': 'MISTIC', 'value': 'MISTIC_score'},
                    {'label': 'popEVE', 'value': 'popEVE'},
                    # === NEW PREDICTORS ADDED TO DROPDOWN ===
                    {'label': 'BayesDel noAF', 'value': 'BayesDel_nsfp33a_noAF'},
                    {'label': 'VARITY R LOO', 'value': 'VARITY_R_LOO'}
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
    
    # Hidden div to store clicked point data
    html.Div(id='clicked-point-data', style={'display': 'none'}),
    
    # Hidden div to store the variant string for copying
    dcc.Store(id='variant-to-copy', storage_type='memory'),
    
    # Copy notification toast
    dbc.Toast(
        "Variant copied to clipboard!",
        id="copy-toast",
        header="Success",
        is_open=False,
        dismissable=True,
        duration=2000,
        icon="success",
        style={"position": "fixed", "top": 66, "right": 10, "width": 350, "zIndex": 9999},
    ),
    
    # Modal for variant details
    dbc.Modal([
        dbc.ModalHeader(dbc.ModalTitle(id='modal-title')),
        dbc.ModalBody(id='modal-body'),
        dbc.ModalFooter(
            dbc.Button("Close", id="close-modal", className="ms-auto", n_clicks=0)
        ),
    ],
    id="variant-modal",
    is_open=False,
    size="lg"
    ),
    
    html.Div([
        html.H3("Legend:", style={'marginTop': 5}),
        html.Ul([
            html.Li([html.Span("★", style={'color': 'gold', 'fontSize': '20px'}), 
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
            f"<b>Change:</b> {row.get('aa_change', 'N/A')}<br>"
            f"<b>{score_column}:</b> {score_str}<br>"
        )
        
        # Add MISTIC prediction if available and relevant
        if 'MISTIC_pred' in row and pd.notna(row.get('MISTIC_pred')) and score_column == 'MISTIC_score':
            text += f"<b>MISTIC_pred:</b> {row['MISTIC_pred']}<br>"
        
        # Add CLNREVSTAT for ClinVar variants
        if row.get('source') == 'ClinVar' and pd.notna(row.get('CLNREVSTAT')):
            text += f"<b>Review Status:</b> {row['CLNREVSTAT']}<br>"
        
        text += (
            f"<b>AC_genomes:</b> {row.get('AC_genomes', 0)}<br>"
            f"<b>AC_joint:</b> {row.get('AC_joint', 0)}<br>"
            f"<b>nhomalt_genomes:</b> {row.get('nhomalt_genomes', 0)}<br>"
            f"<b>nhomalt_joint:</b> {row.get('nhomalt_joint', 0)}<br>"
            "<b><i>Click for more details</i></b>"
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
    
    # Store all data for click events
    all_plot_data = pd.DataFrame()
    
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
            all_plot_data = pd.concat([all_plot_data, gnomad_filtered])
            
            # Create color scale based on threshold field
            color_values = gnomad_filtered[threshold_field].astype(float).replace(0, 1e-9)
            log_color = np.log10(color_values)
            
            fig.add_trace(go.Scatter(
                x=gnomad_filtered['aa_position'],
                y=gnomad_filtered[score],
                mode='markers',
                name='gnomAD Variants',
                customdata=gnomad_filtered.to_dict('records'),
                marker=dict(
                    size=10,
                    color=log_color,
                    colorscale='GnBu',
                    showscale=True,
                    colorbar=dict(
                        title=f"{threshold_field}",
                        x=1,
                        thickness=12
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
            
            if not clinvar_filtered .empty:
                all_plot_data = pd.concat([all_plot_data, clinvar_filtered])
                fig.add_trace(go.Scatter(
                    x=clinvar_filtered['aa_position'],
                    y=clinvar_filtered[score],
                    mode='markers',
                    name='ClinVar P/LP',
                    customdata=clinvar_filtered.to_dict('records'),
                    marker=dict(
                        color='darkred', 
                        size=10,
                        symbol='diamond',
                        opacity=0.7
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
                all_plot_data = pd.concat([all_plot_data, custom_filtered])
                fig.add_trace(go.Scatter(
                    x=custom_filtered['aa_position'],
                    y=custom_filtered[score],
                    mode='markers',
                    name='Custom Variant',
                    customdata=custom_filtered.to_dict('records'),
                    marker=dict(
                        color='gold',
                        size=16,
                        symbol='star',
                        line=dict(width=1, color='white')
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

@app.callback(
    Output("variant-modal", "is_open"),
    Output("modal-title", "children"),
    Output("modal-body", "children"),
    Output("variant-to-copy", "data"),
    Input("pathogenicity-plot", "clickData"),
    Input("close-modal", "n_clicks"),
    State("variant-modal", "is_open"),
    State('score-dropdown', 'value'),
)
def toggle_modal(clickData, n_close, is_open, current_score):
    ctx = callback_context
    if not ctx.triggered:
        return False, "", "", ""
    
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if trigger_id == "pathogenicity-plot" and clickData:
        # Get the clicked point's data
        point_data = clickData['points'][0].get('customdata', {})
        
        if not point_data:
            return False, "", "", ""
        
        # Build VCF string
        vcf_string = f"{point_data.get('chrom', '')}-{point_data.get('pos', '')}-{point_data.get('ref', '')}-{point_data.get('alt', '')}"
        
        # Build URL based on source
        url = ""
        if point_data.get('source') == 'ClinVar' and point_data.get('variant_id'):
            url = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{point_data['variant_id']}/"
        elif point_data.get('source') == 'gnomAD':
            url = f"https://gnomad.broadinstitute.org/variant/{vcf_string}?dataset=gnomad_r4"
        
        # Build modal title
        modal_title = f"Variant Details: {point_data.get('aa_change', 'N/A')}"
        
        # Create modal body with all information
        modal_content = [
            html.Div([
                html.H5("Variant Information"),
                html.Hr(),
                html.Div([
                    dbc.Alert([
                        dbc.Row([
                            dbc.Col([
                                html.Strong("hg38: "),
                                html.Span(vcf_string, style={'fontFamily': 'monospace'})
                            ], width="auto", className="align-self-center"),
                            dbc.Col([
                                dbc.Button(
                                    "📋 Copy Variant",
                                    id="copy-variant-btn",
                                    color="secondary",
                                    size="sm",
                                    className="ms-2",
                                    n_clicks=0,
                                    style={'backgroundColor': '#0D6EFD', 'borderColor': '#0D6EFD',}
                                )
                            ], width="auto", className="align-self-center")
                        ], className="g-2 align-items-center")  # g-2 for small gap
                    ], color="light"),
                ]),
                # Add hyperlink
                html.Div([              
                dbc.Button(
                    f"View in {point_data.get('source', 'Database')} →",
                    href=url,
                    target="_blank",
                    color="primary",
                    className="w-100"
                )
                ]),
                html.Br(),
                html.P([html.Strong("Gene: "), point_data.get('gene', 'N/A')]),
                html.P([html.Strong("Transcript: "), point_data.get('transcript', 'N/A')]),
                html.P([html.Strong("AA Change: "), point_data.get('aa_change', 'N/A')]),
                html.P([html.Strong("Source: "), point_data.get('source', 'N/A')]),
            ]),
            html.Hr(),
            html.H5("Pathogenicity Scores"),
            html.Hr(),
        ]
        # Add all available scores
        scores_to_show = ['REVEL', 'am_pathogenicity', 'cadd_v1.7', 'MPC', 
                         'MISTIC_score', 'popEVE', 'BayesDel_nsfp33a_noAF', 'VARITY_R_LOO']
        
        for score in scores_to_show:
            if score in point_data and pd.notna(point_data.get(score)):
                score_val = float(point_data[score])
                score_label = score
                if score == current_score:
                    modal_content.append(
                        html.P([
                            html.Strong(f"{score_label}: "), 
                            html.Span(f"{score_val:.3f}", 
                                    style={'backgroundColor': '#ffc107', 'padding': '2px 6px', 
                                          'borderRadius': '4px', 'fontWeight': 'bold'}),
                            html.Span(" (current)", style={'fontStyle': 'italic', 'marginLeft': '5px'})
                        ])
                    )
                else:
                    modal_content.append(
                        html.P([html.Strong(f"{score_label}: "), f"{score_val:.3f}"])
                    )
        
        # Add ClinVar specific information
        if point_data.get('source') == 'ClinVar':
            if pd.notna(point_data.get('CLNDN')):
                modal_content.append(html.Hr())
                modal_content.append(html.H5("ClinVar Information"))
                modal_content.append(html.Hr())
                modal_content.append(
                    html.P([html.Strong("Disease: "), point_data['CLNDN'].replace('_', ' ')])
                )
            if pd.notna(point_data.get('CLNREVSTAT')):
                modal_content.append(
                    html.P([html.Strong("Review Status: "), point_data['CLNREVSTAT'].replace('_', ' ')])
                )
            if pd.notna(point_data.get('CLNSIG')):
                modal_content.append(
                    html.P([html.Strong("Clinical Significance: "), point_data['CLNSIG'].replace('_', ' ')])
                )
        
        # Add gnomAD counts
        if point_data.get('source') == 'gnomAD' or point_data.get('source') == 'Custom':
            modal_content.append(html.Hr())
            modal_content.append(html.H5("Population Data"))
            modal_content.append(html.Hr())
            modal_content.append(
                html.P([html.Strong("AC Genomes: "), str(point_data.get('AC_genomes', 0))])
            )
            modal_content.append(
                html.P([html.Strong("AC Joint: "), str(point_data.get('AC_joint', 0))])
            )
            modal_content.append(
                html.P([html.Strong("Homozygous Genomes: "), str(point_data.get('nhomalt_genomes', 0))])
            )
            modal_content.append(
                html.P([html.Strong("Homozygous Joint: "), str(point_data.get('nhomalt_joint', 0))])
            )
        
        return True, modal_title, modal_content, vcf_string
    
    elif trigger_id == "close-modal":
        return False, "", "", ""
    
    return is_open, "", "", ""

# Callback to handle copying variant
@app.callback(
    Output("copy-toast", "is_open"),
    Input("copy-variant-btn", "n_clicks"),
    State("variant-to-copy", "data"),
    prevent_initial_call=True
)
def copy_variant(n_clicks, variant_string):
    if n_clicks and variant_string:
        # Note: Actually copying to clipboard requires JavaScript
        # This callback just shows the notification
        # The actual copying is handled by the clientside callback below
        return True
    return False

# Clientside callback to actually copy to clipboard
app.clientside_callback(
    """
    function(n_clicks, variant_string) {
        if (n_clicks > 0 && variant_string) {
            navigator.clipboard.writeText(variant_string);
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("copy-variant-btn", "n_clicks_timestamp"),  # Dummy output
    Input("copy-variant-btn", "n_clicks"),
    State("variant-to-copy", "data"),
    prevent_initial_call=True
)

if __name__ == '__main__':
    app.run(debug=True, port=8050)
