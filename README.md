# EBI Causality Course 2026

## Starting a notebook

Install the environment (first time only):

```bash
pixi install
```

Launch JupyterLab:

```bash
pixi run notebook
```

JupyterLab starts with this project's environment as the kernel. No separate
kernel installation is needed.

## Reusing the signaling plot

`utils.py` contains the reusable theme, CORNETO styling processor,
and a convenience function. From a notebook in this directory:

```python
from utils import plot_signaling

plot_signaling(
    graph,
    vertex_values,
    edge_values=edge_values,
    receptors=receptor_names,
    regulators=tf_names,
    measured=measured_names,
    theme={"edge_width": 0.9},  # optional overrides
)
```

Values must follow the graph's vertex/edge order and represent one sample.
`edge_values` selects nonzero solution edges; alternatively pass `edge_indexes`.
Without either, all edges are shown. Annotations are supplied explicitly; this
module does not download data. Regulators take precedence over receptor
candidates when annotations overlap. Copy the module alongside other notebooks
to reuse it in another project.

Node fills represent inferred activity. Arrowheads retain the prior interaction
sign; dashed grey edges indicate reduced upstream activity, including release
of inhibition. All displayed solution edges remain selected. The notebook
contains the full legend. Graphviz/WASM preserve the biological shapes and bars;
CORNETO's NetworkX renderer may not preserve all these attributes.

Colors, shapes, widths and the activity threshold are in `BIOLOGY_THEME`.
Layout can be changed with `graph_attr={"rankdir": "LR"}`; individual elements
can be overridden with CORNETO's `custom_vertex_attr` and `custom_edge_attr`.

For future CORNETO integration, `biology_signaling_style(graph, data, theme)`
already follows its public processor API. It can be used directly with
`graph.plot(processor=biology_signaling_style, theme=BIOLOGY_THEME, data=...)`.
A future built-in preset would combine this processor, theme and layout defaults;
the theme alone cannot encode release-of-inhibition semantics. Biological class
annotations should remain explicit and independent of input/output roles.
