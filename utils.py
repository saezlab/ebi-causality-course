"""Reusable CORNETO dataframe conversion and biological graph styling.

The processor uses CORNETO's public ``processor(graph, data, theme)`` API.
It does not fetch annotations or depend on CARNIVAL/OmniPath/Decoupler.
Node classes must be supplied by the caller; regulator annotations win overlaps.
Arrowheads encode the prior interaction, dashes encode reduced source activity,
and fills encode inferred node activity. Widths do not encode confidence.
"""

import numpy as np
import pandas as pd

BIOLOGY_THEME = {
    "positive_color": "#F6D8B5",
    "negative_color": "#CCE2F2",
    "zero_color": "#F2F2F2",
    "node_border": "#555555",
    "active_edge": "#555555",
    "reduced_edge": "#888888",
    "edge_width": 1.1,
    "arrow_size": 0.65,
    "node_border_width": 0.8,
    "measured_border_width": 1.8,
    "receptor_shape": "box",
    "regulator_shape": "hexagon",
    "protein_shape": "ellipse",
    "activity_threshold": 0.5,
}


def _activity_vector(values, size, name):
    """Accept one sample in graph order, including CARNIVAL's column vectors."""
    values = np.asarray(values, dtype=float)
    if values.shape == (size, 1):
        values = values[:, 0]
    if values.shape != (size,) or not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain {size} finite values for one sample in graph order.")
    return values


def biology_signaling_style(graph, plot_data, theme):
    theme = {**BIOLOGY_THEME, **theme}
    threshold = float(theme["activity_threshold"])
    values = _activity_vector(plot_data["vertex_values"], len(graph.V), "vertex_values")
    activity = dict(zip(graph.V, values, strict=True))
    receptors = set(plot_data.get("receptors", ()))
    regulators = set(plot_data.get("regulators", ()))
    measured = set(plot_data.get("measured", ()))
    edge_attrs, vertex_attrs = {}, {}

    for vertex in graph.V:
        # Regulator annotations take precedence over the broad LigRecExtra
        # target list, which can also contain intracellular proteins (e.g. SMAD3).
        if vertex in regulators:
            shape, style = theme["regulator_shape"], "filled"
            node_class = "TF / transcriptional regulator"
        elif vertex in receptors:
            shape, style = theme["receptor_shape"], "rounded,filled"
            node_class = "Receptor"
        else:
            shape, style = theme["protein_shape"], "filled"
            node_class = "Other signaling protein"
        value = activity[vertex]
        fill = (theme["positive_color"] if value > threshold else
                theme["negative_color"] if value < -threshold else theme["zero_color"])
        vertex_attrs[vertex] = {
            "shape": shape, "style": style, "fillcolor": fill,
            "color": theme["node_border"],
            "penwidth": str(theme["measured_border_width"] if vertex in measured else theme["node_border_width"]),
            "tooltip": f"{vertex}: {node_class}; inferred activity {value:+g}",
        }

    for index, (sources, targets) in graph.edges():
        # CARNIVAL adds boundary edges with an empty source or target.
        # These are hidden by orphan_edges=False below.
        if len(sources) != 1 or len(targets) != 1:
            continue
        source = next(iter(sources))
        reduced = activity[source] < -threshold
        inhibitory = graph.get_attr_edge(index).get("interaction", 0) < 0
        mechanism = "inhibition" if inhibitory else "activation"
        effect = (f"Reduced {mechanism}" if reduced else f"Increased {mechanism}")
        if reduced and inhibitory:
            effect = "Release of inhibition"
        if abs(activity[source]) <= threshold:
            effect = f"No inferred change in upstream activity ({mechanism})"
        edge_attrs[index] = {
            "arrowhead": "tee" if inhibitory else "normal",
            "style": "dashed" if reduced else "solid",
            "color": theme["reduced_edge"] if reduced else theme["active_edge"],
            "penwidth": str(theme["edge_width"]), "arrowsize": str(theme["arrow_size"]), "tooltip": effect,
        }
    return edge_attrs, vertex_attrs



def plot_signaling(
    graph,
    vertex_values,
    *,
    edge_values=None,
    edge_indexes=None,
    receptors=(),
    regulators=(),
    measured=(),
    theme=None,
    renderer="auto",
    graph_attr=None,
    node_attr=None,
    edge_attr=None,
    **kwargs,
):
    """Plot a signed signaling solution using the biological theme.

    Values must follow ``graph.V`` / ``graph.E`` order. Supply a single sample
    as a vector or column vector; select a column yourself for multiple samples.
    If ``edge_indexes`` is omitted, ``edge_values`` selects edges whose absolute
    value exceeds ``activity_threshold``. With neither argument, show all edges.
    Boundary edges are hidden by default. Colors, shapes and widths can be
    overridden with a partial ``theme`` dictionary. Graphviz attributes and
    sparse ``custom_vertex_attr`` / ``custom_edge_attr`` overrides are forwarded.
    The return value is CORNETO's plot object (Graphviz by default).
    """
    resolved_theme = {**BIOLOGY_THEME, **(theme or {})}
    threshold = float(resolved_theme["activity_threshold"])
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("activity_threshold must be finite and nonnegative.")
    values = _activity_vector(vertex_values, len(graph.V), "vertex_values")
    if edge_indexes is None and edge_values is not None:
        edge_values = _activity_vector(edge_values, len(graph.E), "edge_values")
        edge_indexes = np.flatnonzero(np.abs(edge_values) > threshold)
    kwargs.setdefault("orphan_edges", False)
    return graph.plot(
        processor=biology_signaling_style,
        theme=resolved_theme,
        data={"vertex_values": values, "receptors": receptors,
              "regulators": regulators, "measured": measured},
        edge_indexes=edge_indexes,
        renderer=renderer,
        graph_attr={"rankdir": "TB", "bgcolor": "white", "nodesep": "0.3",
                    "ranksep": "0.5", **(graph_attr or {})},
        node_attr={"fixedsize": "false", "fontname": "Helvetica", "fontsize": "10",
                   "fontcolor": "#222222", "margin": "0.10,0.06", "height": "0.35",
                   **(node_attr or {})},
        edge_attr={"fontname": "Helvetica", "fontsize": "9", **(edge_attr or {})},
        **kwargs,
    )


def plot_dag_solution(
    method,
    *,
    decimals=2,
    renderer="auto",
    graph_attr=None,
    node_attr=None,
):
    """Plot a solved LinearDAGDiscovery model using the course notebook style.

    Accepts a fitted method, or a built method whose problem has been solved.
    Only selected edges are shown. Labels contain signed coefficients in the
    original measurement units; green edges are positive and orange edges
    are negative. The method's public solution API validates assigned values.

    Returns CORNETO's plot object for display in a notebook. ``graph_attr``
    and ``node_attr`` override the default layout and rounded-box styling.
    """
    if not isinstance(decimals, int) or isinstance(decimals, bool) or decimals < 0:
        raise ValueError("decimals must be a nonnegative integer.")
    graph = method.get_solution_graph()
    edge_styles = {}
    for index in range(graph.num_edges):
        coefficient = float(graph.get_attr_edge(index)["coefficient"])
        if not np.isfinite(coefficient):
            raise ValueError("Selected edges must have finite fitted coefficients.")
        edge_styles[index] = {
            "label": f"{coefficient:+.{decimals}f}",
            "color": "#16856b" if coefficient > 0 else "#dc7448" if coefficient < 0 else "#43566b",
            "arrowhead": "normal",
        }
    return graph.plot(
        renderer=renderer,
        graph_attr={
            "rankdir": "LR", "ranksep": "0.65", "nodesep": "0.35",
            "splines": "spline", "pad": "0.2", "bgcolor": "transparent",
            **(graph_attr or {}),
        },
        node_attr={
            "shape": "box", "style": "rounded,filled", "fillcolor": "#edf3f8",
            "color": "#43566b", "fontcolor": "#213247", "fontname": "Helvetica",
            "fontsize": "10", "margin": "0.16,0.09", **(node_attr or {}),
        },
        edge_attr={"fontname": "Helvetica", "fontsize": "10", "fontcolor": "#344054",
                   "penwidth": "2.4"},
        custom_edge_attr=edge_styles,
    )


def to_df(data_obj):
    data = data_obj.to_dict()
    return pd.DataFrame.from_dict(
        {
            obs: {f["id"]: f["value"] for f in item["features"]}
            for obs, item in data.items()
        },
        orient="index"
    )


def from_df(df, *, sample_id_column=None, feature_columns=None):
    """Convert a samples-by-features dataframe to CORNETO ``Data``.

    Use the index as sample IDs, or supply ``sample_id_column`` for a CSV table.
    By default, exclude the course CSV metadata columns: sample_id, regime,
    target, level, intervention_value, and intervention_type. Explicit
    ``feature_columns`` can select the measurements and their order.

    When target/intervention_type columns are present, preserve hard
    interventions on the target feature and use regime as intervention_group.
    This helper supports observational and hard-intervention data only.
    ``to_df`` exports values alone, so that round trip cannot restore metadata
    unless it is supplied separately in the input dataframe.
    """
    from corneto.data import Data

    if not df.columns.is_unique:
        raise ValueError("Dataframe columns must be unique.")
    table = df.set_index(sample_id_column) if sample_id_column is not None else df
    if not table.index.is_unique or table.index.hasnans:
        raise ValueError("Sample IDs must be unique and non-missing.")
    metadata_columns = {
        "sample_id", "regime", "target", "level", "intervention_value", "intervention_type",
    }
    columns = list(feature_columns) if feature_columns is not None else [
        column for column in table.columns if column not in metadata_columns
    ]
    if not columns or len(columns) != len(set(columns)):
        raise ValueError("Select at least one measurement column, without duplicates.")
    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise ValueError(f"Unknown measurement columns: {missing_columns!r}.")
    has_target = "target" in table.columns
    has_intervention = "intervention_type" in table.columns
    if has_target != has_intervention:
        raise ValueError("Provide both target and intervention_type columns, or neither.")

    samples = {}
    for sample_id, row in zip(table.index, table.to_dict(orient="records"), strict=True):
        intervention = row.get("intervention_type", "none")
        if pd.isna(intervention) or intervention == "":
            intervention = "none"
        if intervention not in {"none", "hard"}:
            raise ValueError(f"Unsupported intervention type {intervention!r} for sample {sample_id!r}.")
        target = row.get("target")
        if intervention == "hard" and target not in columns:
            raise ValueError(f"Hard-intervention target {target!r} is not a measurement column.")
        features = {}
        for column in columns:
            attributes = {"mapping": "vertex", "value": float(row[column])}
            if intervention == "hard" and column == target:
                attributes["intervention"] = "hard"
                group = row.get("regime")
                if group is not None and not pd.isna(group) and group != "":
                    attributes["intervention_group"] = group
            features[column] = attributes
        samples[sample_id] = features
    return Data.from_cdict(samples)



def evaluate_dag(method, data, *, groups=None):
    """Diagnose one solved DAG on complete measurements, without refitting.

    Groups are a Series/mapping indexed by sample ID, or are inferred from
    intervention_group labels. All measured variables must occur in the DAG.
    Hard-clamped targets are excluded. Returns summary, by_variable, by_regime,
    responses, scales and sample groups. Normalized errors use evaluation-data standard
    deviations excluding clamped targets (a constant response uses scale 1).
    Mean-response errors give every available group/variable equal weight.

    Forward predictions propagate assigned interventions; local residuals use
    observed parents. Group means estimate expectations, not individual outcomes.
    Callers must choose held-out samples if they want a generalization score.
    """
    observed = to_df(data).astype(float)
    if observed.empty or not np.isfinite(observed.to_numpy()).all():
        raise ValueError("Evaluation data must contain complete, finite measurements.")
    if groups is None:
        labels = {}
        for sample_id, sample in data.samples.items():
            group_values = {
                feature.data["intervention_group"] for feature in sample.features
                if feature.data.get("intervention_group") is not None
            }
            if len(group_values) > 1:
                raise ValueError("A sample has multiple intervention groups; supply groups explicitly.")
            is_intervened = any(
                feature.data.get("intervention", "none") != "none"
                or feature.data.get("intervened", False) for feature in sample.features
            )
            if is_intervened and not group_values:
                raise ValueError("Intervened samples need intervention_group labels or explicit groups.")
            labels[sample_id] = next(iter(group_values), "observation")
        groups = pd.Series(labels)
    groups = pd.Series(groups).reindex(observed.index)
    if groups.isna().any():
        raise ValueError("Provide a group label for every evaluation sample ID.")

    evaluation = method.evaluate(data)
    predicted = to_df(evaluation["predictions"])
    if set(predicted.columns) != set(observed.columns):
        raise ValueError("The DAG must include every measured variable, including isolated vertices.")
    predicted = predicted.reindex(index=observed.index, columns=observed.columns).astype(float)
    clamped = pd.DataFrame.from_dict({
        sample_id: {feature.id: bool(feature.data.get("clamped", False))
                    for feature in sample.features}
        for sample_id, sample in evaluation["predictions"].samples.items()
    }, orient="index").reindex(index=observed.index, columns=observed.columns)
    valid = ~clamped.astype(bool)
    if not np.isfinite(predicted.to_numpy()).all():
        raise ValueError("The fitted method must produce finite predictions for all variables.")
    if not valid.to_numpy().any():
        raise ValueError("No non-clamped measurements are available for scoring.")
    scales = observed.where(valid).std(ddof=0)
    scales = scales.where(np.isfinite(scales) & scales.gt(np.finfo(float).eps), 1.0)
    forward_error = (observed - predicted).where(valid)
    local_error = to_df(evaluation["residuals"]).reindex(
        index=observed.index, columns=observed.columns,
    ).astype(float).where(valid)
    if not np.isfinite(local_error.to_numpy()[valid.to_numpy()]).all():
        raise ValueError("The method lacks local residuals for some scored measurements.")

    responses = []
    for variable in observed.columns:
        frame = pd.DataFrame({
            "observed": observed[variable].where(valid[variable]),
            "predicted": predicted[variable].where(valid[variable]),
            "regime": groups,
        }).dropna(subset=["observed", "predicted"])
        grouped = frame.groupby("regime", sort=False).agg(
            observed_mean=("observed", "mean"),
            predicted_mean=("predicted", "mean"),
            observed_sem=("observed", "sem"),
            n=("observed", "size"),
        ).reset_index()
        grouped["variable"] = variable
        grouped["mean_response_nmae"] = (
            grouped["predicted_mean"] - grouped["observed_mean"]
        ).abs() / scales[variable]
        responses.append(grouped)
    responses = pd.concat(responses, ignore_index=True)
    scored = pd.DataFrame({
        "variable": np.tile(observed.columns, len(observed)),
        "regime": np.repeat(groups.to_numpy(), len(observed.columns)),
        "forward_nmae": (forward_error.abs() / scales).to_numpy().reshape(-1),
        "local_nmae": (local_error.abs() / scales).to_numpy().reshape(-1),
    }).dropna(subset=["forward_nmae", "local_nmae"])
    by_variable = scored.groupby("variable", sort=False)[["forward_nmae", "local_nmae"]].mean()
    by_variable.insert(0, "mean_response_nmae", responses.groupby("variable")["mean_response_nmae"].mean())
    by_regime = scored.groupby("regime", sort=False)[["forward_nmae", "local_nmae"]].mean()
    by_regime.insert(0, "mean_response_nmae", responses.groupby("regime")["mean_response_nmae"].mean())
    return {
        "summary": pd.Series({
            "mean_response_nmae": responses["mean_response_nmae"].mean(),
            "forward_nmae": scored["forward_nmae"].mean(),
            "local_nmae": scored["local_nmae"].mean(),
            "selected_edges": len(method.get_selected_edge_indices()),
            "scored_values": len(scored),
        }),
        "by_variable": by_variable,
        "by_regime": by_regime,
        "responses": responses,
        "scales": scales,
        "groups": groups,
    }


def plot_dag_errors(method, data, *, groups=None, metric="mean_response_nmae"):
    """Plot one solved method's errors by variable and intervention group.

    Accepts mean_response_nmae (default), forward_nmae, or local_nmae.
    No fitting is performed; all directly clamped targets are excluded.
    """
    import matplotlib.pyplot as plt

    titles = {
        "mean_response_nmae": "Normalized error in group means",
        "forward_nmae": "Normalized forward absolute error",
        "local_nmae": "Normalized local absolute error",
    }
    if metric not in titles:
        raise ValueError(f"Choose a metric from {list(titles)!r}.")
    report = evaluate_dag(method, data, groups=groups)
    tables = [report["by_variable"], report["by_regime"]]
    fig, axes = plt.subplots(1, 2, figsize=(11, max(4, .4 * max(len(t) for t in tables))),
                             constrained_layout=True)
    maximum = max(float(table[metric].max()) for table in tables)
    for ax, table, title in zip(axes, tables, ["Error by variable", "Error by intervention group"], strict=True):
        values = table[metric].sort_values()
        ax.barh(range(len(values)), values, color="#16856b", height=.65)
        ax.set_yticks(range(len(values)), labels=values.index)
        ax.set(title=title, xlabel=titles[metric], xlim=(0, max(maximum * 1.2, .01)))
        for row, value in enumerate(values):
            ax.text(value + max(maximum * .02, .0001), row, f"{value:.3f}", va="center", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("DAG prediction errors (lower is better)")
    return fig


def plot_dag_predictions(method, data, *, groups=None, variables=None):
    """Show observed means, predicted means, and signed errors as heatmaps.

    Rows are experimental groups and columns are variables. All panels share
    a zero-centered color scale, so zero responses do not obscure the effects
    elsewhere. Gray cells mark excluded hard-intervention targets. Errors are
    predicted minus observed means, in original measurement units. Group means
    remain empirical estimates; their sampling errors are not plotted here.
    """
    import matplotlib.pyplot as plt

    report = evaluate_dag(method, data, groups=groups)
    responses = report["responses"]
    variables = list(variables) if variables is not None else list(responses["variable"].unique())
    if not variables or not set(variables).issubset(set(responses["variable"])):
        raise ValueError("Choose variables with non-clamped evaluation responses.")
    regimes = list(report["groups"].drop_duplicates())
    observed = responses.pivot(index="regime", columns="variable", values="observed_mean").reindex(
        index=regimes, columns=variables,
    )
    predicted = responses.pivot(index="regime", columns="variable", values="predicted_mean").reindex(
        index=regimes, columns=variables,
    )
    tables = [observed, predicted, predicted - observed]
    limit = max(float(np.nanmax(np.abs(table.to_numpy()))) for table in tables)
    limit = max(limit, .1)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e4e7ec")
    fig, axes = plt.subplots(
        1, 3, figsize=(max(10.5, 3 * (.4 * len(variables) + 1.1)),
                       max(3.5, .34 * len(regimes) + 1.8)),
        sharey=True, constrained_layout=True,
    )
    titles = ["Observed group means", "Predicted group means", "Error: predicted − observed"]
    for ax, table, title in zip(axes, tables, titles, strict=True):
        im = ax.imshow(np.ma.masked_invalid(table.to_numpy()), cmap=cmap,
                       vmin=-limit, vmax=limit, aspect="auto")
        ax.set_xticks(range(len(variables)), labels=variables, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(regimes)), labels=regimes, fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.set_xticks(np.arange(-.5, len(variables), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(regimes), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=.5)
        ax.tick_params(which="minor", bottom=False, left=False)
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].set_ylabel("Experimental group")
    fig.colorbar(im, ax=axes, label="Mean score / difference in mean score", shrink=.8)
    fig.suptitle("Intervention responses · gray = directly clamped target (excluded)", fontsize=12)
    return fig
