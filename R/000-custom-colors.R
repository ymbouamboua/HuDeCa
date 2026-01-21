
# HuDeCA custom colors for A Single-Cell and Spatial Atlas of Early Human Olfactory Development paper

output <- "/Users/yvon.mbouamboua/Documents/projects/singlecell/hudeca/analysis/"

custom_colors <- list()

# Level 2 — Fine Cell Types (ann_level_2)

custom_colors$ann_level_2 <- c(
  
  # Epithelium
  "CycBasal" = "#C2A523",
  "Deut"     = "#3498db",
  "GBC"      = "#926B54",
  "INP"      = "#f05b43",
  "iOSN"     = "#33b8ff",
  "MCC"      = "#1f618d",
  "MV"       = "#efe13c",
  "OHBC"     = "#E41A1C",
  "RHBC"     = "#C82C73",
  "SUS"      = "#C09ACA",
  
  # Neuronal
  "GABA"     = "#800EF1",
  "GLUT"     = "#706fd3",
  "GnRH"     = "#2EECDB",
  "NP"       = "#6A0B78",
  "NOS1"     = "#95819F",
  
  # Neural Crest / Glia
  "Mel"      = "#eb10fd",
  "NCC"      = "#E3D9AC",
  "OEC"      = "#A90D55",
  "SC"       = "#95ccff",
  "SCP"      = "orangered",
  
  # Mesenchyme
  "Cartl"    = "#0B4B19",
  "Mes0"     = "#99D6A9",
  "Mes1"     = "#1B8F76",
  "Mes2"     = "#9DAF07",
  "MSC"      = "#27F557",
  "Osteo"    = "#4CAD4C",
  
  # Vascular
  "EC"       = "#E788C2",
  "Lymph"    = "#F78896",
  "Peri"     = "#BBD870",
  
  # Muscle
  "Sat"      = "#CB7647",
  "SM"       = "#E69F00",
  
  # Immune
  "MG"       = "#91BFB7"
)

# Lineage Subsets (clean + grouped) 
custom_colors$OE <- c(
  "CycBasal" = "#C2A523",
  "GBC"      = "#926B54",
  "INP"      = "#f05b43",
  "iOSN"     = "#33b8ff",
  "MV"       = "#efe13c",
  "OHBC"     = "#E41A1C",
  "SUS"      = "#C09ACA"
)

custom_colors$neuronal_lineage <- c(
  "CycBasal" = "#C2A523",
  "GBC"      = "#926B54",
  "INP"      = "#f05b43",
  "iOSN"     = "#33b8ff",
  "OHBC"     = "#E41A1C"
)

custom_colors$mv_lineage <- c(
  "CycBasal" = "#C2A523",
  "MV"       = "#efe13c",
  "OHBC"     = "#E41A1C"
)

custom_colors$sus_lineage <- c(
  "CycBasal" = "#C2A523",
  "OHBC"     = "#E41A1C",
  "SUS"      = "#C09ACA"
)


# Level 1 — Lineages (ann_level_1) 
custom_colors$ann_level_1 <- c(
  "Immune"               = "#91BFB7",
  "Mesenchyme"           = "#27F557",
  "Muscle"               = "#E69F00",
  "Neuronal"             = "#800EF1",
  "NeuralCrest_Glia"     = "#FF4500",
  "Olfactory Epithelium" = "#E41A1C",
  "Respiratory Epithelium" = "#3498DB",
  "Vascular"             = "#E788C2"
)


# Level 0 — Supergroups (ann_level_0) 
custom_colors$ann_level_0 <- c(
  "Epithelium"        = "#C82C73",
  "Immune"            = "#91BFB7",
  "Mesenchyme"        = "#1B8F76",
  "Muscle"            = "#E69F00",
  "NeuralCrest_Glia"  = "#A90D55",
  "Neuronal"          = "#706FD3",
  "Vascular"          = "#BBD870"
)

# Sample Colors (Developer Review Names)
custom_colors$sample <- c(
  "S1-PCW7"   = "#6baed6",
  "S2-PCW8"   = "#1f78b4",
  "S3-PCW10"  = "#b2df8a",
  "S4-PCW10.5"= "#ff7f00",
  "S5-PCW10"  = "#33a02c",
  "S6-PCW12"  = "#e31a1c",
  "S7-PCW12"  = "#fb9a99",
  "S8-PCW12"  = "#6a3d9a"
)

# Developmental Stage Colors (PCWs)   
custom_colors$PCW <- c(
  "PCW7"   = "#6baed6",
  "PCW8"   = "#1f78b4",
  "PCW10"  = "#b2df8a",
  "PCW10.5"= "#33a02c",
  "PCW12"  = "#e31a1c"
)

# Combined Stage Groups (PCW groups) 
custom_colors$group <- c(
  "PCW7_8" = "#1f78b4",
  "PCW10"  = "#33a02c",
  "PCW12"  = "#e31a1c"
)

# BMI-Related Groups
custom_colors$BMI <- c(
  "PCW7"      = "#6baed6",
  "PCW8"      = "#1f78b4",
  "PCW10"     = "#33a02c",
  "OB-PCW10.5"= "#ff7f0e",
  "OB-PCW12"  = "#e31a1c",
  "PCW12"     = "#e31a1c"
)

# Sex
custom_colors$sex <- c(
  "Female" = "#b15928",
  "Male" = "#6a3d9a"
)

# Cell Cycle / Phase
custom_colors$phase <- c(
  "G1"         = "#9CD4FA",
  "G2M"        = "#f1c232",
  "Non-Cycling"= "gray50",
  "PostM"      = "#6aa84f",
  "S"          = "#90301B"
)

# Cell State
custom_colors$cell_state <- c(
  "Cycling"     = "#E41A1C",
  "G1-like"     = "#FF7F00",
  "Postmitotic" = "#4DAF4A",
  "Quiescent"   = "#377EB8"
)

saveRDS(custom_colors, file.path(output, "custom_colors.rds"))
