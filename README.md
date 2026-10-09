# Yeast-maps
This repository contains the metabolic maps of [Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) that [Metabolic Atlas](https://metabolicatlas.org) serves, in four formats, and the files from which the first maps were made.

The maps match **Yeast-GEM 9.1.1**.

## Content

- `subsystem/`: 102 maps, one per subsystem of the model, 27 of them for the transport subsystems (`Transport [a, b]`, as `transport_<a>_<b>`); exchange reactions have no map.
- `compartment/`: 19 maps, one per compartment. The cytoplasm, the endoplasmic reticulum membrane and the mitochondrion are each divided over several maps, each holding whole subsystems.

Each folder holds the same maps in each format:

| Folder | Format | Content |
| --- | --- | --- |
| `svg/` | SVG | The map as Metabolic Atlas shows it. |
| `sbgn/` | [SBGN-ML](https://sbgn.github.io) 0.3, process description | Metabolites as simple chemicals, reactions as processes, genes as macromolecules catalysing them, at the positions of the SVG. Validated against the libsbgn schema. |
| `sbml/` | [SBML](https://sbml.org) Level 3 Version 1, with the layout and groups packages | The reactions on the map with their stoichiometry, reversibility and genes (as modifiers) from Yeast-GEM 9.1.1; the layout of every drawn metabolite, gene and edge; one group per subsystem. Checked with libsbml. |
| `escher/` | [Escher](https://escher.github.io) map (JSON, schema 1-0-0) | Metabolites and reactions where the SVG draws them, edges as Bézier curves through the SVG's bends, with stoichiometry, reversibility and gene rules from Yeast-GEM 9.1.1; the title and compartment headings as text labels (Escher has no compartment boxes or gene boxes). Open it in Escher with *Load map JSON*; loading the model as well (COBRA JSON) shows names and lets you overlay data. Checked against Escher's schema and its own consistency checks. |

Gene boxes show the gene name over its systematic (ORF) name.

## How the maps are made

The first maps were drawn for Yeast8 (Yeast-GEM 8.3) in CellDesigner; `SBMLfiles/` holds those drawings, with the whole map in `SBMLfiles/Yeast8.xml`, and `scripts/` the files used to split them by subsystem. Since Yeast-GEM 9.1.1 the maps are fitted to each model release by script, keeping the drawing: identifiers follow the model, removed reactions and metabolites are taken off, reactions of a map's subsystem or compartment that were not drawn are added, maps for new subsystems and for the compartments are drawn, and parts that share a metabolite are connected. The script and its rules are in [MetabolicAtlas/data-generation](https://github.com/MetabolicAtlas/data-generation) (`maps`, rules in `maps/RULES.md`). The SVG maps that Metabolic Atlas serves are kept in [MetabolicAtlas/data-files](https://github.com/MetabolicAtlas/data-files) (`svg/Yeast-GEM`), and the SBGN, SBML and Escher files are written from them with `maps/publish_maps.py`. After each model update in data-files, a workflow there opens a pull request here with the new maps.

The custom map (chloroalkane and chloroalkene degradation) is only in data-files.

This repository is administered by Angelo Limeta ([@angelolimeta](https://github.com/angelolimeta)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.
