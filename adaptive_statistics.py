import os
from pathlib import Path

import pandas as pd
import plotly.express as px
from neurosnap.protein import Protein

### Prepare CSV:
## Read results the CSV
path = Path("results")
df = pd.read_csv(path / "result_summary.csv")

## Map expression to numeric values
expression_map = {"high": 3, "medium": 2, "low": 1, "none": 0}
# Fill in 'none' where missing and encode
df["expression_numeric"] = df["expression"].fillna("none").map(expression_map)

## add contact count
contacts = []
for name in df.name:
    pdb_path = path / "structures" / (name + ".pdb")
    print(pdb_path)
    if os.path.exists(pdb_path):
        prot = Protein(str(pdb_path))
        contacts.append(prot.calculate_contacts_interface("A", "B"))  # NOTE: Assumes everything is dimer
    else:
        print("skipped")
        contacts.append(-1)  # NOTE: Assumes everything is dimer

df["contacts"] = contacts

## Add destress columns
df_destress = pd.read_csv(path / "destress_binder_with_egfr.csv")

df_indexed = df.set_index('name')
df_destress_indexed = df_destress.set_index('design_name')

# Find columns in df_destress that are not already in df
new_columns = [col for col in df_destress.columns if col not in df.columns and col != 'design_name']

# Subset and align df_destress
df_destress_subset = df_destress_indexed[new_columns]

# Join the new columns to df
df_combined = df_indexed.join(df_destress_subset, how='left')

# Reset index if you want the original 'name' as a column
df_combined.reset_index(inplace=True)
df = df_combined



### Correlation matrix of the entire results
# Compute correlation matrix with encoded expression
corr_matrix = df.corr(method = "spearman", numeric_only=True)

# Plot the heatmap
fig = px.imshow(corr_matrix, text_auto=True, color_continuous_scale="RdBu_r", title="Correlation Matrix Heatmap")
fig.update_layout(xaxis_title="Features", yaxis_title="Features", coloraxis_colorbar=dict(title="Correlation"))
fig.show()

# Plot the heatmap (JUST KD)
fig = px.bar(abs(corr_matrix.kd).dropna(), text_auto=True, title="Correlation With KD")
fig.update_layout(xaxis_title="Features", yaxis_title="Absolute Corelation", coloraxis_colorbar=dict(title="Correlation"))
fig.show()


## Correlation matrix of the results segregated by method used
# import pandas as pd
# import plotly.express as px
# import ast

# # Load CSV
# df = pd.read_csv("results/result_summary.csv")

# # Safely parse the JSON-like strings to actual lists, skip NaNs
# def safe_eval(x):
#     if pd.isna(x):
#         return []
#     try:
#         return ast.literal_eval(x)
#     except:
#         return []

# df['design_models'] = df['design_models'].apply(safe_eval)

# # Explode the design_models into one model per row
# df_exploded = df.explode('design_models')

# # Drop rows with no models after exploding
# df_exploded = df_exploded[df_exploded['design_models'] != '']

# # Get unique model names
# unique_models = df_exploded['design_models'].dropna().unique()

# # Plot correlation matrix for each model
# for model in unique_models:
#     sub_df = df_exploded[df_exploded['design_models'] == model]

#     numeric_df = sub_df.select_dtypes(include='number')

#     if numeric_df.shape[0] < 2:
#         print(f"Skipping {model}: not enough data.")
#         continue

#     corr_matrix = numeric_df.corr()

#     fig = px.imshow(
#         corr_matrix,
#         text_auto=True,
#         color_continuous_scale='RdBu_r',
#         title=f"Correlation Matrix for {model}"
#     )

#     fig.update_layout(
#         xaxis_title="Features",
#         yaxis_title="Features",
#         coloraxis_colorbar=dict(title="Correlation")
#     )

#     fig.show()


# ## Show bar plots for which methods had the best results
# import pandas as pd
# import plotly.express as px
# import ast

# # Load the CSV
# df = pd.read_csv("results/result_summary.csv")

# # Convert design_models to actual lists, handle NaNs safely
# def safe_eval(x):
#     if pd.isna(x):
#         return []
#     try:
#         return ast.literal_eval(x)
#     except:
#         return []

# df['design_models'] = df['design_models'].apply(safe_eval)

# # Explode the lists into individual rows
# df_exploded = df.explode('design_models')
# df_exploded = df_exploded[df_exploded['design_models'] != '']  # Remove empty rows

# # Drop rows where binding_strength or expression is missing
# df_exploded = df_exploded.dropna(subset=['binding_strength', 'expression'])

# # --------------------
# # PLOT 1: Binding Strength per Model
# # --------------------
# bind_counts = df_exploded.groupby(['design_models', 'binding_strength']).size().reset_index(name='count')

# fig1 = px.bar(
#     bind_counts,
#     x='design_models',
#     y='count',
#     color='binding_strength',
#     barmode='group',
#     title="Binding Strength per Design Model"
# )
# fig1.update_layout(xaxis_title="Design Model", yaxis_title="Count", xaxis_tickangle=-45)
# fig1.show()

# # --------------------
# # PLOT 2: Expression Level per Model
# # --------------------
# expr_counts = df_exploded.groupby(['design_models', 'expression']).size().reset_index(name='count')

# fig2 = px.bar(
#     expr_counts,
#     x='design_models',
#     y='count',
#     color='expression',
#     barmode='group',
#     title="Expression Level per Design Model"
# )
# fig2.update_layout(xaxis_title="Design Model", yaxis_title="Count", xaxis_tickangle=-45)
# fig2.show()


## Show bar plots for which methods had the best results (NORMALIZED)
# import pandas as pd
# import plotly.express as px
# import ast

# # Load CSV
# df = pd.read_csv("results/result_summary.csv")

# # Safely parse design_models column
# def safe_eval(x):
#     if pd.isna(x):
#         return []
#     try:
#         return ast.literal_eval(x)
#     except:
#         return []

# df['design_models'] = df['design_models'].apply(safe_eval)

# # Explode design_models lists into rows
# df_exploded = df.explode('design_models')
# df_exploded = df_exploded[df_exploded['design_models'] != '']

# # Drop rows with missing target data
# df_exploded = df_exploded.dropna(subset=['binding_strength', 'expression'])

# # --------------------
# # PLOT 1: Stacked Binding Strength (%)
# # --------------------
# bind_counts = df_exploded.groupby(['design_models', 'binding_strength']).size().reset_index(name='count')
# bind_total = bind_counts.groupby('design_models')['count'].transform('sum')
# bind_counts['percent'] = (bind_counts['count'] / bind_total) * 100

# fig1 = px.bar(
#     bind_counts,
#     x='design_models',
#     y='percent',
#     color='binding_strength',
#     barmode='relative',  # Stack bars
#     title="Binding Strength Distribution per Design Model (Normalized)"
# )
# fig1.update_layout(
#     xaxis_title="Design Model",
#     yaxis_title="Percentage",
#     xaxis_tickangle=-45,
#     yaxis_ticksuffix="%",
# )
# fig1.show()

# # --------------------
# # PLOT 2: Stacked Expression Level (%)
# # --------------------
# expr_counts = df_exploded.groupby(['design_models', 'expression']).size().reset_index(name='count')
# expr_total = expr_counts.groupby('design_models')['count'].transform('sum')
# expr_counts['percent'] = (expr_counts['count'] / expr_total) * 100

# fig2 = px.bar(
#     expr_counts,
#     x='design_models',
#     y='percent',
#     color='expression',
#     barmode='relative',  # Stack bars
#     title="Expression Level Distribution per Design Model (Normalized)"
# )
# fig2.update_layout(
#     xaxis_title="Design Model",
#     yaxis_title="Percentage",
#     xaxis_tickangle=-45,
#     yaxis_ticksuffix="%",
# )
# fig2.show()
